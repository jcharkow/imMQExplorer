import copy 
import numpy as np
import pandas as pd
from ncls import NCLS
import sqlite3
from imMQExplorer.rangeTable import *
import imMQExplorer.timsdata as td




''' Class to store MS2 windows'''
class Experiment():
    def __init__(self, expr_path):
        exp = td.TimsData(expr_path)
        conn = exp.conn
        self.ms1 = Ms1Frames(conn=conn)
        self.ms2 = Ms2Frames(conn=conn)

''' stores the MS1 frames of an experiment (start/stop time) '''
class Ms1Frames():
    def __init__(self, conn, frameWidth=0.1):
        frameTimeQuery = conn.execute("select Frames.Time as startTime, Frames.Time + ? as endTime, Id FROM Frames WHERE MsMsType == 0 order by Id", [frameWidth])  
        self.timeTable = FrameRange(frameTimeQuery)

'''stores properties on MS2 frame conn can be a sql query connection or a pandas dataframe.pkl with information on each frame'''    
class Ms2Frames():
    def __init__(self, conn, frameWidth=0.1):
        if type(conn) is sqlite3.Connection:

            frameTimeQuery = conn.execute("select Frames.Time as startTime, Frames.Time + ? as endTime, Frame, ScanNumBegin from PasefFrameMsMsInfo, Frames WHERE PasefFrameMsMsInfo.Frame == Frames.Id order by Frame, ScanNumBegin", [frameWidth]) #take frame and scanNum begin because this is the primary key, want to make sure that all data get is in the same order 
            
            frameIMQuery = conn.execute("select ScanNumBegin, ScanNumEnd, Frame from PasefFrameMsMsInfo order by Frame, ScanNumBegin")

            frameMZQuery = conn.execute("select IsolationMz - (IsolationWidth/2) as startMz, IsolationMz + (IsolationWidth/2) as endMz, ScanNumBegin, Frame from PasefFrameMsMsInfo order by Frame, ScanNumBegin") #select ScanNumBegin and Frame because this is the primary key
            
            self.timeTable = FrameRange(frameTimeQuery)
            self.imTable = FrameRange(frameIMQuery)
            self.mzTable = FrameRange(frameMZQuery)

        else: #initiate using pandas dataframe
            self.timeTable = RangeArray(conn['rtStart'], conn['rtEnd'])
            self.mzTable = RangeArray(conn['mzStart'], conn['mzEnd'])
            
            if 'scanNumBegin' in conn.columns.values: #means IM values set
                self.imTable = RangeArray(conn['scanNumBegin'], conn['scanNumEnd'])
            else:
                self.imTable = None

    '''saves the object as a pandas dataframe'''
    def save_pandas(self, file_path):
        df = pd.DataFrame(np.column_stack([self.timeTable.idx, self.timeTable[:,1:], self.mzTable[:,1:], self.imTable[:,1:]]), columns=['idx','rtStart','rtEnd', 'mzStart', 'mzEnd', 'scanNumBegin', 'scanNumEnd']) #have to exclude index from all the tables so don't have repeat index
        df.to_pickle(file_path)


''' class for storing peptide features '''
class Precursor():

        def __init__(self, filePath, rawFile = None, ionFWHM = False, retentionFWHM = False, **kwargs):
                '''
                This function initiates the precursor object
                filePath = the path where the precursors.txt file from maxQuant was generated
                rawFile = if set, precursors list will be filtered to only include raw files
                ionFWHM = if TRUE will base boundries off of Ion Mobility Length (FWHM) instead of Ion Mobility Length
                retentionFWHM = if TRUE will base boundries off of Retention Length (FWHM) instead of Retention Length
                **kwargs for pd.read_csv
                '''
                self.precData = pd.read_csv(filePath, '\t', **kwargs) #read the allpeptides.txt file

                #filter out data to only include precursors in raw file (if raw file specified)
                if not rawFile is None:
                        # before filtering make sure raw file exists 
                        self.precData = self.precData[self.precData['Raw file'] == rawFile]
                        try:
                            assert(len(self.precData) > 0)
                        except AssertionError:
                            raise AssertionError("Raw file is not in the list of features")

                if retentionFWHM:
                    self.retentionTable = RangeTable(self.precData['Retention time'], self.precData['Retention length (FWHM)'])
                else:
                    self.retentionTable = RangeTable(self.precData['Retention time'], self.precData['Retention length'])

                if ionFWHM:
                    self.ionMobilityTable = RangeTable(self.precData['Ion mobility index'], self.precData['Ion mobility index length (FWHM)'])
                else:
                    self.ionMobilityTable = RangeTable(self.precData['Ion mobility index'], self.precData['Ion mobility index length'])

                self.mzTable = MzTable(self.precData, idx = self.precData.index)

                self.chargeTable = np.column_stack((np.arange(0, len(self.precData)), self.precData['Charge'].values))
                self.__vecMzImOverlap = np.vectorize(self.__mzImOverlap, excluded=['idx2_data']) #vectorized mz_im_overlap function
                print(self.precData)

        def __len__(self):
            return len(self.retentionTable)

        '''helper to the overlap function. From two indices provided determine if there is overlap in the mz and im dimensions. Other is the data which the indicies of idx2 correspond to '''
        def __mzImOverlap(self, idx1, idx2, idx2_data = None):
            if idx2_data is not None:
                mz_overlap = self.mzTable.idx_overlap(idx1, idx2, yData = idx2_data.mzTable)
                
                if mz_overlap: #if overlap in mz look for overlap in im
                    return self.ionMobilityTable.idx_overlap(idx1, idx2, yData = idx2_data.imTable)
                else:
                    return False #if have not returned then no overlap (means no mz_overlap) else: #idx2 data not set so don't pass anything 
            else:
                mz_overlap = self.mzTable.idx_overlap(idx1, idx2) 
                
                if mz_overlap: #if overlap in mz look for overlap in im
                    return self.ionMobilityTable.idx_overlap(idx1, idx2) 
                else:
                    return False #if have not returned then no overlap (means no mz_overlap)
        
        ''' 
        this method maps the "fake" index back to the real one
        rslt = pairwise overlap both using the fake index
        true_index the true index
        '''
        @staticmethod
        def unindex(rslt, true_index):
            for fake_idx, true_idx  in enumerate(true_index):
                rslt[rslt == fake_idx] = true_idx
            return rslt


        ''' determine overlap with im and without im '''
        def overlap(self, im = True, decimal_places = 5):
            #if index not linear then filter and use this hidden index
            use_hidden = False #if true that means reindex done for overlap (with have to unindex before return results)
            if not np.all(self.retentionTable.idx == np.arange(0,len(self.retentionTable.start))):
                hidden_idx = np.arange(0,len(self.retentionTable.start))
                use_hidden = True
            else:
                hidden_idx = self.retentionTable.idx

            #decimal_places = accuracy of the retentionTime axis (how many decimal places should take into acoount)
            
            # convert the retention times to integers so that it is compatible with ncls
            ret_int_start = (self.retentionTable.start * (10**decimal_places)).astype(int) 
            ret_int_end = (self.retentionTable.end * (10**decimal_places)).astype(int)

            #create the ncls object
            ncls = NCLS(ret_int_start, ret_int_end, hidden_idx)

            #find all pairwise retention time overlaps, store in a vertical 2XN numpy nd array
            ret_idx = np.column_stack(ncls.all_overlaps_both(ret_int_start,ret_int_end, hidden_idx)) #column stack puts the two lists vertically (easier iteration for np.vectorize)

            
            #filter out pairs where x=y, although these overlap not interested in them
            ret_idx = ret_idx[ret_idx[:,0] != ret_idx[:,1]]

            if ret_idx.size > 0: #only look for overlap if there is overlap in retention time
                #if im flag on, then have to check for overlap in both mz and im dimensions
                if im:
                    rslt = ret_idx[self.__vecMzImOverlap(ret_idx[:,0], ret_idx[:,1])]
                else:
                    rslt = ret_idx[self.mzTable.vec_idx_overlap(ret_idx[:,0], ret_idx[:,1])]
            else:
                rslt = np.array([])
            
            #unindex if need to 
            if use_hidden:
                return Precursor.unindex(rslt, self.retentionTable.idx)
            else:
                return rslt
        ''' calculates the pariwise overlaps for a subset of the precursors, possible use: calculate precursor density of a window '''
        def overlap_subset(self, idx, im = True, decimal_places = 5):
           #create a new precursor object with only the filtered idx (retain the idx of the original precursor object
           filtered = self.filter_precursor(idx)
           return filtered.overlap(im = im)

            

        '''overlap the given precursor against different ranges not specified here: e.g. double charge against all others, query should be Precursor object, e.g. overlap against DDA windows --> DDA windows is other'''
        def overlap_against_other(self, query, im = True, decimal_places = 5):
            #decimal_places = accuracy of the retentionTime axis (how many decimal places should take into acoount)

            # convert the retention times to integers so that it is compatible with ncls
            ret_int_start = (self.retentionTable.start * (10**decimal_places)).astype(int) 
            ret_int_end = (self.retentionTable.end * (10**decimal_places)).astype(int)

            query_ret_int_start = (query.retentionTable.start * (10**decimal_places)).astype(int) 
            query_ret_int_end = (query.retentionTable.end * (10**decimal_places)).astype(int) 
            
            #create the ncls object
            ncls = NCLS(ret_int_start, ret_int_end, self.retentionTable.idx)

            #find all pairwise retention time overlaps, store in a vertical 2XN numpy nd array
            ret_idx = np.column_stack(ncls.all_overlaps_both(query_ret_int_start,query_ret_int_end, query.retentionTable.idx)) #column stack puts the two lists vertically (easier iteration for np.vectorize)            
            print(ret_idx)
            #reverse columns so have the retention time idx first then the frame idx
            ret_idx = np.flipr(ret_idx)
            print(ret_idx)
            
            #can't filter out because the index are not the same 
                #ret_idx = ret_idx[ret_idx[:,0] != ret_idx[:,1]]
            
            #if im flag on, then have to check for overlap in both mz and im dimensions
            if im:
                return ret_idx[self.__vecMzImOverlap(ret_idx[:,0], ret_idx[:,1])]
            else:
                return ret_idx[self.mzTable.vec_idx_overlap(ret_idx[:,0], ret_idx[:,1])]
        
        '''overlap of precursors against a window object 
            win - Experiment object'''
        def overlap_against_window(self, exp, im=True, decimal_places = 5):
            print("starting overlap against window")
            ms2 = exp.ms2
            #decimal_places = accuracy of the retentionTime axis (how many decimal places should take into acoount)

            # convert the retention times to integers so that it is compatible with ncls
            ret_int_start = (self.retentionTable.start * (10**decimal_places)).astype(int) 
            ret_int_end = (self.retentionTable.end * (10**decimal_places)).astype(int)

            ms2_time_int_start = (ms2.timeTable.start * (10**decimal_places)).astype(int) 
            ms2_time_int_end = (ms2.timeTable.end * (10**decimal_places)).astype(int) 

            print(ret_int_start)
            print(ret_int_end)
            
            #create the ncls object
            ncls = NCLS(ret_int_start, ret_int_end, self.retentionTable.idx)

            #find all pairwise retention time overlaps, store in a vertical 2XN numpy nd array
            ret_idx = np.column_stack(ncls.all_overlaps_both(ms2_time_int_start,ms2_time_int_end, ms2.timeTable.idx)) #column stack puts the two lists vertically (easier iteration for np.vectorize)            
            ret_idx = np.fliplr(ret_idx)
            
            #if im flag on, then have to check for overlap in both mz and im dimensions
            if im:
                return ret_idx[self.vec_mz_im_overlap(ret_idx[:,0], ret_idx[:,1], idx2_data=ms2)]
            else:
                #print(ms2.mzTable)
                return ret_idx[self.mzTable.vec_idx_overlap(ret_idx[:,0], ret_idx[:,1], yData=ms2.mzTable)]
 
        def filter_precursor(self, idx):
            '''filter based on idx, return new precursor object '''

            #make a copy of precursor object
            new = copy.deepcopy(self)
            
            new.retentionTable = self.retentionTable.filter(idx)
            new.ionMobilityTable = self.ionMobilityTable.filter(idx)
            new.mzTable = self.mzTable.filter(idx)
            new.chargeTable = new.chargeTable[np.isin(new.chargeTable[:,0],idx)]
            return new
