# This contains classes for dealing with ranges

import numpy as np
from ncls import NCLS


'''This class is a numpy array with start, end index and a table with both start and end '''
class RangeTable():
    ''' default way to start and end is using create is using center and length '''
    def __init__(self, center, length, idx = None):
        '''
                center = pandas series of center value of the measurement
                length = pandas series of length value of the measurment 
                idx = Index
                **center and length must be in the same unit**
        '''
        print("welcome to rangeTable!")
        self.start = np.array(center - (length / 2))
        self.end = np.array(center + (length / 2))
        self.idx = None #initialized below
        self.table = None #initialized below
        self.initialize(idx)

    ''' initiate supplimentary variables (idx, table and vectorized overlap) '''
    def initialize(self, idx = None):
        '''
        idx = pandas series of specific idx values for the RangeTable if none provided will use (1:len_array)
        table = table with start and end in the same object
        vec_idx_overlap = vectorized version of the idx_overlap method
        '''
        if idx is None:
            self.idx = np.arange(0, len(self.start))
        else:
            self.idx = idx.astype(int).values

        ''' create a table as a np ndarray '''
        self.table = np.column_stack((self.idx, self.start, self.end))

        #vectorize the idx_overlap function so can run with vectors
        self.vec_idx_overlap = np.vectorize(self.idx_overlap, excluded =['yData'])


    def __str__(self):
        return str(self.table)

    #allows for easy access of table
    def __repr__(self):
        return repr(self.table)

    def __getitem__(self, key):
        return self.table[key]

    def __len__(self):
        return len(self.table)

    ''' given two indecies, determine if the intervals at these indices overlap. If y data == None x and y from same dataset however if y data provided then will fetch y interval from y data '''
    def idx_overlap(self, x, y, yData = None):
        '''determine overlap between two list element'''
        if yData == None:
            return self.start[x] <= self.end[y] and self.start[y] <= self.end[x]
        else:
            return self.start[x] <= yData.end[y] and yData.start[y] <= self.end[x]
    
    ''' using ncls determine the overlap across a subset or the entire axes '''
    def ncls_overlap(self, decimal_places = 5, start_idx = 0, end_idx = None):
        # if end_idx is none set as the end of the list
        if end_idx is None:
            end_idx == len(self.start)
        #decimal_places = accuracy of the retentionTime axis (how many decimal places should take into acoount) (since NCLS works optimally with integers)

        # convert the retention times to integers so that it is compatible with ncls
        int_start = (self.start * (10**decimal_places)).astype(int) 
        int_end = (self.end * (10**decimal_places)).astype(int)

        #create the ncls object
        ncls = NCLS(int_start, int_end, self.idx)

        #find all pairwise retention time overlaps, store in a vertical 2XN numpy nd array

        return np.column_stack(ncls.all_overlaps_both(int_start[start_idx:end_idx], int_end[start_idx:end_idx], self.idx[start_idx:end_idx])) #column stack puts the two lists vertically (easier iteration for np.vectorize)
    
    ''' return a filtered rangeTable to the correponding indices, element returned is a deep copy of the original object '''
    def filter(self, idx):
        idx = np.array(idx) #convert idx to numpy array


        # deepcopy the rangeTable
        new = copy.deepcopy(self)
        
        #filter the rangeTable based on index
        new.table = new.table[np.isin(new.table[:,0], idx)]
        
        #set other variables 
        new.start = new.table[:,1]
        new.end = new.table[:,2]
        new.idx = new.table[:,0].astype(int) #make sure integer for ncls not to throw error

        return new #return the new RangeTable object
        





########## Child Classes ##############

''' a type of range table that can be initialized from array like elements specfying the start and stop values '''
class RangeArray(RangeTable):
    def __init__(self, start, end):
        self.start = np.array(start)
        self.end = np.array(end)
        self.idx = None #initialized in super
        self.table = None #initialized in super
        super().initialize()

''' range table specific for Frames (initiated from sql query) '''
class FrameRange(RangeTable):
    def __init__(self, query):
        rslt = np.array(query.fetchall())
        self.start = rslt[:,0]
        self.end = rslt[:,1]
        self.idx = None #initialized in super
        self.table = None #initialized in super
        super().initialize()
     
''' Class which holds the minimum and maximum Mz values '''
class MzTable(RangeTable):
        def __init__(self, precFile, **kwargs):
            # where store all information in form (minMz, maxMz, idx) 
            self.start = precFile['m/z'].values #the start values are just the provided m/z
            
            self.idx = np.arange(0,len(precFile))
            
            #initialize the maximum m/z values
            self.end = np.zeros(len(precFile))
            
            # fill in the maximum m/z values
            
            # m/z and Number of isotopic peaks column name is invalid for itertuples() rename them (note Charge column is valid and does not have to be renamed)
            precFile = precFile.rename(columns={"m/z":"mz", "Number of isotopic peaks":"numPeaks"})

            for idx, row in enumerate(precFile.itertuples()):
                self.end[idx] = getattr(row, 'mz') + ((getattr(row, 'numPeaks'))/getattr(row, 'Charge')) #note: should not have -1 because monoisotopic peak is not included in numPeaks (ran tests Mar. 3 2020)                
           
            super().initialize(**kwargs)


