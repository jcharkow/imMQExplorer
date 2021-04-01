# This contains functions for composing theoretical isolation windows based on a list of features (MQ allPeptides.txt output)

# This script a creates pandas dataframe with DIA windows. 
# DIA windows correspond with each precursor however are not centered around the precursor.


# Ret time --> Retention time of precursor apex with 0.1 second width (0.01 on either side) --> Is centered around the precursor, assume the cycling speed is going to be fast enough that this may be possible


# Mz time --> 25 m/z not centered around the precursor 

# IM time --> 0.3, 0.6 or 0.9 1/k0 not centered around the precursor
    #NOTE: 0.1 1/k0 ~= 1 scan. Therefore 0.3 ~= 30 scans

from imMQExplorer.featureExperiment import *
import numpy as np
import argparse
import pandas as pd


class TheoreticalIsolationWindows():
        ''' 
        class for creating Theoretical isolation windows from a list of features
        pep = path to allpeptides.txt file
        mode = method in which windows shold be created (DDA/DIA)

        '''
        def __init__(self, pep, rawFile=None, mode='DIA'):
            self.pep = Precursor(pep, rawFile=rawFile, ionFWHM=True, retentionFWHM=True )
            self.windows = None

            self.imRange = (0, 918) # for 100ms ramp time
            self.rtWidth = 0.1
            ### default DIA parameters
            if mode.upper() == "DIA":
                self.mzWidth = 25
                self.imWidth = 250
            elif mode.upper() == 'DDA':
                self.mzWith = 2
                self.imWith = 25 
            print(self.pep.precData)

        # return dictionary of parameters 
        def getParams(self):
            return {'mzWidth':self.mzWidth, 'imRange':self.imRange, 'imWidth':self.imWidth, 'rtWidth':self.rtWidth}

        def setParams(self, param, value):
            if param == 'mzWidth':
                self.mzWidth = value
            elif param == 'imRange':
                if isinstance(value, tuple):
                    self.imRange = value
                else:
                    raise ValueError("imRange must be a tuple")
            elif param == 'imWidth':
                self.imWidth = value
            
            elif param == 'rtWidth':
                self.rtWidth = value
            else:
                raise ValueError("Invalid Parameter")


        def createTheoreticalWindows(self, mode):

            # template to store windows
            windows = np.empty((len(self.pep), 6)) #columns are rtStart, rtEnd, mzStart, mzEnd, scanNumBegin, scanNumEnd

            # retention time windows create same way for both DDA and DIA
            rtWindows = np.column_stack([self.pep.precData['Retention time'] - (self.rtWidth / 2), (self.pep.precData['Retention time'] + (self.rtWidth / 2))])
            windows[:,0:2] = rtWindows #set the rtStart and rtEnd

            if mode.upper() == "DIA":
                
                # create a set of possible "windows" that can be assigned to the feature
                mzBins = TheoreticalIsolationWindows.__createBins(self.pep.mzTable.start, self.pep.mzTable.end, self.mzWidth)
                imBins = TheoreticalIsolationWindows.__createBins(self.pep.ionMobilityTable.start, self.pep.ionMobilityTable.end, self.imWidth)

                ### correct windows if values exceed maximum ##
                imBins  = TheoreticalIsolationWindows.__correctBins(imBins, self.imRange[0], self.imRange[1]) 


                # for each precursor find the correct bin
                for idx, row in self.pep.precData[['m/z', 'Ion mobility index']].iterrows():
                    # set the mzStart and end
                    windows[idx,2:4] = TheoreticalIsolationWindows.__findBin(mzBins, row['m/z'])


                    #set the im start and end 
                    windows[idx, 4:6] = TheoreticalIsolationWindows.__findBin(imBins, row['Ion mobility index'])

            elif mode.upper() == 'DDA':
                windows[:, 2:4] = np.column_stack([(self.pep.precData['m/z'] - (self.mzWidth / 2)).values, (self.pep.precData['m/z'] + (self.mzWidth / 2))])
                windows[:, 4:6] = np.column_stack([np.floor(self.pep.precData['Ion mobility index'] - (self.imWidth / 2)), np.floor(self.pep.precData['Ion mobility index'] + (self.imWidth / 2))])


            self.windows =  pd.DataFrame(windows, columns=['rtStart', 'rtEnd', 'mzStart', 'mzEnd', 'scanNumBegin', 'scanNumEnd'])
            return self.windows

        ''' create numpy array of start and end for each window'''
        @staticmethod
        def __createBins(start, end, width):
            #note endWin is not inclusive 
            startWin = (min(start) // width) * width
            endWin = (max(end) // width) * width + width #have to add width so includes the maximum values (because // cuts off the remainder and need to technically round up for this measurment)
            diaWindows = np.column_stack([np.arange(startWin, endWin, width), np.arange(startWin + width, endWin + width, width)]).astype(int) #windows for the dia width of width and range from min to max of the m/z in the pep retention table

            return diaWindows

        @staticmethod
        def __correctBins(window, minVal, maxVal):
            ''' for an array of window start and ends, "correct" for start/end values outside of a certain range '''
            starts = window[:,0]
            ends = window[:,1]

            starts[ (starts < minVal) ] = minVal
            ends[ (ends > maxVal) ] = maxVal

            return np.column_stack([starts, ends])


        # given a list of ranges find which one the value belongs to
        @staticmethod
        def __findBin(bins, val):
            rslt =  bins[(val >= bins[:,0]) & (val < bins[:,1])] #filter is inclusive at the start and not inclusve at the end 
            assert(rslt.size == 2)
            return rslt

