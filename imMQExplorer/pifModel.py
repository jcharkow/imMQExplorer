# This python script creates a class for the categorical model which predicts the number of peptides identified with PIF. 

import numpy as np
import pandas as pd


class pifModel:
    # pif - pif of the non theoretical peptides (e.g. ones with IM), outputted from pifComputations (in pd dataframe)
    # pep - allpeptides.txt output from maxquant (in pd dataframe)
    # evi - evidence.txt output from maxquant (in pd dataframe)
    def __init__(self, pif, pifTopScans, pep, evi):
        self.pep = pep
        self.evi = evi
        self.identified = self.getIdentified().index.drop_duplicates() #only 1 pif per precursor
        self.pif = self.mergePif(pif, pifTopScans)

        self.mergePif(pif, pifTopScans)
        self.model = None

        self.collapsePif()
        self.createModel()

    def mergePif(self, pif, pifTopScans):
        # combine top scans and regular pif into one table
        pifTopScans  = pifTopScans.drop(columns=['targetIntensity','totIntensIm','topScansIdx'])
        pif = pif.drop(columns=['TargetIntensityNoIm', 'TargetIntensityIm', 'totIntensIm', 'totIntensNoIm'])
        pif = pd.merge(pifTopScans, pif)

        return pif

    # filters the allPeptides table to those identified peptides (by comparing Ms/MS scans with identified
    def getIdentified(self):
        return self.pep[self.pep['MS/MS scan number'].isin(self.evi['MS/MS scan number'].values)]

    def collapsePif(self):
        # to be the most "conservative" sort by noIm PIF first 
        self.pif =  (self.pif[['Feature', 'ImPIF', 'noImPIF', 'ImPIFTop13Scans']]
                                 .drop_duplicates()
                                 .sort_values(by=['noImPIF', 'ImPIF'], ascending=False)
                                 .groupby('Feature')
                                 .head(1))
        self.pif = self.pif.sort_values(by='Feature')

    '''apply the model to this column name (in the PIF data)'''
    def applyModel(self, colName):
        # Expected Number identified without IM
        appliedModel = pifModel.bin_pif(self.pif[colName], name='lenTar')
        appliedModel['lenIdent'] = appliedModel['lenTar'] * self.model['percentIdent']
        return appliedModel


    ''' function to get number identified '''
    def getNumIdentified(self, what='im'):
        # gets the number
        if what == 'im':
            return self.model['lenIdent'].sum()
        elif what == 'noIm':
            return self.applyModel('noImPIF')['lenIdent'].sum()
        elif what == 'imTop':
            return self.applyModel('ImPIFTop13Scans')['lenIdent'].sum()
        else:
            raise ValueError("'what' parameter must be ['im', 'noIm', 'imTop']")
        
    @staticmethod
    def bin_pif(series, name='counts'):
        # Input:
            # SeriesName - columnName in pif want to bin
            # name - name of binned column
        binned = pd.DataFrame(pd.cut(series, bins =np.arange(0,1.1,0.1), include_lowest = True))
        binned[name] = binned.index
        binned = binned.groupby(series.name).count()
        binned = binned.reset_index().rename(columns={series.name:'interval'})
        binned['binCenter'] = [ round(i.mid,2) for (_, i) in binned['interval'].iteritems() ]
        return binned

    def createModel(self, colName='ImPIF'):
       # Percent (and Number) identified with IM
        self.model = pd.merge(pifModel.bin_pif(self.pif[colName], name='lenTar'), 
                             pifModel.bin_pif(self.pif[self.pif['Feature'].isin(self.identified)][colName], name='lenIdent'))
        self.model['percentIdent'] = self.model['lenIdent'] / self.model['lenTar']

        #modify model so that don't observe "dip at high pif"
        self.model.loc[9:, 'percentIdent'] = self.model.loc[8, 'percentIdent']

'''
       # Percent (and Number) identified with IM
        self.model = pd.merge(bin_pif(self.pif['ImPIF'], name='lenTar'), 
                             bin_pif(self.pif[self.pif['Pre'].isin(identified)]['ImPIF'], name='lenIdent'))
        self.model['percentIdent'] = pifImTarCol['lenIdent'] / pifImTarCol['lenTar']


        #modify model so that don't observe "dip at high pif"
        self.model.loc[9:, 'percentIdent'] = self.model.loc[8, 'percentIdent']

        # Expected Number identified without IM
        pifNoImTarCol = bin_pif(pif_targeted_collapse['noImPIF'], name='lenTar')
        pifNoImTarCol['lenIdent'] = pifNoImTarCol['lenTar'] * pifImTarCol['percentIdent']

        # Expected Number identified with optimized IM
        pifOpImTarCol = bin_pif(pif_targeted_collapse['ImPIFTop13Scans'], name='lenTar')
        pifOpImTarCol['lenIdent'] = pifOpImTarCol['lenTar'] * pifImTarCol['percentIdent']

''' 
