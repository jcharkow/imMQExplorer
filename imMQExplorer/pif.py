# This script combines all the functions used in dda_pif in order to make it easier to run from a different script
# This october update uses the IM from MQ, requries additional data in the table (MQImStart, MQImEnd) corresponding with MQ start and end.  
# Note using IM index length (not the FWHM length)

# This has been updated from the may one for compatibility with new timspy package (my extension of it)
#from myTimsDDA import myTimsDDA
import imMQExplorer.timsdata as td
import numpy as np
import pandas as pd
import sqlite3
''' This class contains all the important information in calculaing MS1 overlap
    Input:
        exp_dir = the .d folder
        meta = precomuted metaData, needs to have columns Feature | MS1Frame | ScanNumBegin | ScanNumEnd | Mz | Charge | NumberOfIsotopicPeaks
'''
class PIFOverlapCalculations:
    def __init__(self, exp_dir, meta, lim = None):
        self.exp = td.TimsData(exp_dir)
        self.meta = meta
        self.pifCalc = False #stores if pifCalculations were done
        self.lim = lim # used for testing purposes trim metadata for faster runtime 
        
    
    ''' calculates the mz range to a given ppm
    Input: mz = target mz
          ppm = ppm value
    '''
    # collapse the meta so that there is one row per MS1-Feature Combination, also sorts the metaData by MS1Frame to prevent double loading later (more efficient analysis)
    def collapseAndSortMeta(self):
        return self.meta[['Feature','MS1Frame', 'ScanNumBegin', 'ScanNumEnd', 'Mz', 'Charge', 'NumberOfIsotopicPeaks', 'IsolationMzBegin', 'IsolationMzEnd', 'MQImStart', 'MQImEnd']].drop_duplicates().sort_values(by='MS1Frame')[:self.lim]
    
    # given a collapsed form of the metadata, uncollpase it with the added columns
    def uncollapse(self, collapsed):
        self.meta = pd.merge(collapsed, uncollapsed)
        
    
    def savePIFMeta(self, save_path):
        self.meta.to_pickle(save_path) 

    # updates the meta table to include PIF results based on collasedPIF table
    def __updateMeta(self, collapsedPIF):
        self.meta = collapsedPIF.merge(self.meta)
        self.pifCalc = True


    ''' MetaFunction for calculating the PIF.
        noIM = Whether to calculate PIF with no IM (across all scans)
        im = whether to calculate PIF with IM (across certain scans)
    '''
    def calcPIF(self, noIm = True, im = True, target = True):

        # 1. collapse and sort data
        collapsedMeta = self.collapseAndSortMeta()
        
        # 0. Create numpy arrays to store results
        noImPIFArray = np.empty(len(collapsedMeta))
        imPIFArray = np.empty(len(collapsedMeta))
        targetImIntensArray = np.empty(len(collapsedMeta))
        targetNoImIntensArray = np.empty(len(collapsedMeta))
        totalImIntensArray = np.empty(len(collapsedMeta))
        totalNoImIntensArray = np.empty(len(collapsedMeta))


        # 2. Load the scan if needed
        previousSpectrum = None
        previousFeature = None
        for idx, row in enumerate(collapsedMeta.itertuples()):

            if idx % 10000 == 0:
                print(idx)

            #2. Define the feature 
            feature = Feature(row)

            #3. If the feature MS1 is equal to the loaded MS1 then don't need to load
            if (previousSpectrum is None) or (not feature.ms1Frame == previousSpectrum.frame): #need to load
                spectrum = Spectrum(feature.ms1Frame, feature, self.exp)
            else: # don't need to load (however link to new feature so can filter)
                spectrum = Spectrum(previousSpectrum.data, feature)


            #4. Compute values needed for PIF
            targetImIntensArray[idx] = spectrum.getTargetIntensity(feature,im = True)
            targetNoImIntensArray[idx] = spectrum.getTargetIntensity(feature, im = False)
            totalImIntensArray[idx] = spectrum.getTotIntens(im = True)
            totalNoImIntensArray[idx] = spectrum.getTotIntens(im = False)

        
        #using the array of values calculate the PIF values 
        collapsedMeta['ImPIF'] = targetImIntensArray / totalImIntensArray
        collapsedMeta['noImPIF'] = targetNoImIntensArray / totalNoImIntensArray

        #include extra information 
        collapsedMeta['TargetIntensityNoIm'] = targetNoImIntensArray
        collapsedMeta['TargetIntensityIm'] = targetImIntensArray
        collapsedMeta['totIntensIm'] = totalImIntensArray
        collapsedMeta['totIntensNoIm'] = totalNoImIntensArray

        #update metaTable 
        self.__updateMeta(collapsedMeta)

    ''' computes the PIF of the X number of top scans (only done with IM)
    input: NumScans: How many scans to compute PIF from (e.g. 13 is 50% of the scans when have 25 scans)'''
    def calcPIFTopScans(self, numScans):

        # 1. collapse and sort data
        collapsedMeta = self.collapseAndSortMeta()
        
        # 0. Create numpy arrays to store results
        imPIFArray = np.empty(len(collapsedMeta))
        targetIntensArray = np.empty(len(collapsedMeta))
        totalIntensArray = np.empty(len(collapsedMeta))
        topScansIdx = np.empty(len(collapsedMeta), dtype=object)
        
        # 2. Load the scan if needed
        previousSpectrum = None
        previousFeature = None
        for idx, row in enumerate(collapsedMeta.itertuples()):
            if idx % 10000 == 0:
                print(idx)

            #2. Define the feature 
            feature = Feature(row)

            #3. If the feature MS1 is equal to the loaded MS1 then don't need to load
            if (previousSpectrum is None) or (not feature.ms1Frame == previousSpectrum.frame): #need to load
                spectrum = Spectrum(feature.ms1Frame, feature, exp = self.exp)
            else: # don't need to load (however link to new feature so can filter)
                spectrum = Spectrum(previousSpectrum.data, feature)


            # 5. Compute the PIF per scan and get best quality scans 
            topScans = spectrum.getBestQualityScans(numScans, feature)
            
            #6. Compute values needed for PIF with only the topScans

            targetIntensArray[idx] = spectrum.getTargetIntensity(feature, im=True, scans=topScans)
            totalIntensArray[idx] = spectrum.getPartialTotIntens(topScans)
            topScansIdx[idx] = topScans - feature.scanNumBegin # subtract scanNumBegin so that all topScans in common axes 


        # 6. using the array of values calculate the PIF values 
        collapsedMeta['ImPIFTop{}Scans'.format(numScans)] = targetIntensArray / totalIntensArray

        #include extra information 
        collapsedMeta['targetIntensity'] = targetIntensArray
        collapsedMeta['totIntensIm'] = totalIntensArray
        collapsedMeta['topScansIdx'] = topScansIdx.tolist()

        #update metaTable 
        self.__updateMeta(collapsedMeta)



# a helper class this is holds information on the target feature
class Feature:
    #takes a pd dataframe row and gets information on the expected mz Ranges
    def __init__(self,pre, ppm=20): 
        self.isotopicMzs = self.calcIsotopicMzs(pre, ppm)
        self.isolationMzBegin = pre.IsolationMzBegin
        self.isolationMzEnd = pre.IsolationMzEnd
        self.scanNumBegin = int(pre.ScanNumBegin)
        self.scanNumEnd = int(pre.ScanNumEnd)
        self.ms1Frame = int(pre.MS1Frame)
        self.mqImStart = int(pre.MQImStart)
        self.mqImEnd = int(pre.MQImEnd)

    #takes a tuple with information from pd dataframe and returns the istopic mass ranges (in interval index)
    def calcIsotopicMzs(self, pre, ppm):
        # Calculate the target m/z signals (theoretical signals based on charge, number of isotopic peaks from maxQuant data)
        tar_mzs = [] #will contain range of the isotopic peaks for the precursor
        for isotope in range(0, int(pre.NumberOfIsotopicPeaks) + 1):
            tar_mzs.append(Feature.__ppmToRange(pre.Mz + (isotope/pre.Charge), ppm))

        return pd.IntervalIndex.from_tuples(tar_mzs, closed='neither')
            
    ''' calculates the mz range to a given ppm
    Input: mz = target mz
          ppm = ppm value
    '''
    @staticmethod
    def __ppmToRange(mz, ppm):
        err_da = mz / (1000000 / ppm)
        return (mz - err_da, mz + err_da)

# class for the loaded spectrum data, contains functions for PIF computations
class Spectrum:
    def __init__(self, data, feature, exp = None):
        # data can be the loaded data or the data want to load find out which one
        if isinstance(data, int):
            # load the data 
            if exp == None:
                raise Exception("Must set the experiment to load the data")
            self.data = self.__loadData(exp, feature)

        elif isinstance(data, pd.DataFrame):
            self.data = data
        else:
            raise Exception("Invalid Dataype")

        # save filtered data here so don't have to recompute 
        self.imData = None # all signals from the IM isolation window 
        self.targetData = None # all target signals based off the MQ start and end (target signals based on MQ not filtered by IM isolation window)

        # compute filtered data here
        self.__getFilteredData(feature)
 
    '''loads MS1 data as a pandas df only loads data in the target m/z window inclusive '''
    def __loadData(self, exp, feature):
         data = exp.readPandas(feature.ms1Frame)
         mz = data['mz'].values
         return data[np.logical_and(mz >= feature.isolationMzBegin, mz <= feature.isolationMzEnd)] #logical_and 3X faster than pandas between

    # filter the data to that stemming from the target precursor 
    def __getFilteredData(self, feature):
        # filter to the isolation window #inclusive on lower boundary but not on the upper boundary
        scans = self.data['scan'].values
        self.targetData = self.data[np.logical_and(scans >= feature.mqImStart, scans <= feature.mqImEnd)] #filter to IM scans where target feature elutes (MQ range)

        # filter to the target mzs
        tarMzBoolMask = [ np.any(feature.isotopicMzs.contains(i)) for i in self.targetData['mz'].values ]
        self.targetData = self.targetData[tarMzBoolMask]

        self.imData = self.data[np.logical_and(scans >= feature.scanNumBegin, scans <= feature.scanNumEnd)] # get the IM data (all peaks within IM isolation window)
    # gets a list of the best X quality scans as computed by PIF 
    # scan is considered "best" if in the top X PIF (when PIF calculated per scan)
    # return a minimum of numBest scans, if there are more scans that can be used (e.g. all 1) then use more scans
    def getBestQualityScans(self, numBest, feature):
        # compute the target Intensity 

        # tot intens by scan only of scans in QP isolation window
        scan = self.imData['scan'].values
        isolationWindowTot = self.imData[np.logical_and(scan >= feature.scanNumBegin, scan <= feature.scanNumEnd)] #filter target signals to only be in the scans actually isolated (If MQ scan is smaller than isolation window this step does not filter out any)
        totIntensByScan = isolationWindowTot[['scan', 'i']].groupby('scan').sum()

        # target intens by scan 
        scan = self.targetData['scan'].values
        isolationWindowTarget = self.targetData[np.logical_and(scan >= feature.scanNumBegin, scan <= feature.scanNumEnd)] #filter target signals to only be in the scans actually isolated (If MQ scan is smaller than isolation window this step does not filter out any)
        tarIntensByScan = isolationWindowTarget[['scan', 'i']].groupby('scan').sum()

        # PIF by scan
        pif = (tarIntensByScan / totIntensByScan).sort_values(by= 'i', ascending=False).reset_index().values
        #if < 13 scans add more scans

        # return top scans list 
        # check how many scans are 100% purity
        pureScans = pif[pif[:,1] == 1]
        if pureScans.shape[0] > numBest:
            return pureScans[:,0]
        else:
            return pif[:numBest, 0]
            

    # get the total target intensity 
    # Scans: (array of scans) if provided compute target intensity only across certain scans
    def getTargetIntensity(self, feature, scans=None, im = True):
        if scans is not None:
            data = self.targetData[np.isin(self.targetData['scan'].values, scans)]
        elif im:
            scan = self.targetData['scan'].values
            data  = self.targetData[np.logical_and(scan >= feature.scanNumBegin, scan <= feature.scanNumEnd)] #filter target signals to only be in the scans actually isolated in an IM scan
        else:
            data = self.targetData # all target signals whether isolated or not from an IM scan.

        return data['i'].values.sum()

    # get the total intensity from some scans
    def getPartialTotIntens(self, scans):
        return self.imData[np.isin(self.imData['scan'].values, scans)]['i'].values.sum()

    # get the total intensity (with or without IM filter)
    def getTotIntens(self, im = True):
        if im:
            return self.imData['i'].values.sum()
        else:
            return self.data['i'].values.sum()
