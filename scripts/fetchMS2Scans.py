# This python script fetches the MS2 scans from a PASEF experiment 
import argparse
import pandas as pd
from pyopenms import *
import imMQExplorer.timsdata as td

#convert string argument to bool (using code from stack overflow) https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# parse args
parser = argparse.ArgumentParser(description='fetch all MS2 windows of precursor in MzML format')
parser.add_argument('fIn', help='.d file in')
parser.add_argument('fOut', help='mzml File out')
parser.add_argument('--verbose', help='whether output is verbose', default="True")

args = parser.parse_args()

verbose = str2bool(args.verbose)

# open experiment 
exp = td.TimsData(args.fIn)



############## Get MetaData ####################
# get Precursor table
q = exp.conn.execute("select Id, LargestPeakMz from Precursors")
rslt = q.fetchall()
pre = pd.DataFrame(rslt, columns=['Id','mz'])

#get Pasef table
q = exp.conn.execute("select Frame, ScanNumBegin, ScanNumEnd, Precursor from PasefFrameMsMsInfo")
rslt = q.fetchall()
pasef = pd.DataFrame(rslt, columns=['Frame', 'ScanNumBegin', 'ScanNumEnd', 'Precursor'])
#merge tables
pasef = pasef.sort_values(by='Precursor')
pasef_pre = pd.merge(pasef, pre, left_on='Precursor', right_on='Id').drop(columns='Id')


############## Fetch MS2 Frames ###################
curPrecursor = (0,None)
mzmlExp = MSExperiment()
for frame in pasef_pre.itertuples():
        if verbose and frame.Frame % 10000 == 0:
            print(frame.Frame)
        # If not the same as previous precursor 
        if curPrecursor[0] != frame.Precursor: 
            p = Precursor()
            p.setMZ(frame.mz)
            p.setMetaValue("Id", frame.Precursor)
            curPrecursor = (frame.Precursor, p)

        # get the raw data using timspy module
        raw = exp.readPandas(frame.Frame, scanBegin=frame.ScanNumBegin, scanEnd=frame.ScanNumEnd)

        # set mz and i
        s = MSSpectrum()
        s.set_peaks([raw['mz'].values, raw['i'].values])
        s.setPrecursors([curPrecursor[1]]) # add the precursor

        #set the IM scan
        im = IntegerDataArray()
        im.set_data(raw['scan'].values.astype('int32'))
        im.setName('IonMobilityScan')
        s.setIntegerDataArrays([im])

        ## NOTE: if the data array is empty, upon saving the data array will be erased. 

        # assert that placed integer data array in 
        assert(s.getIntegerDataArrays() != [])

        # add metadata 
        s.setMetaValue('Frame', frame.Frame)
        s.setMetaValue('ScanNumBegin', frame.ScanNumBegin)
        s.setMetaValue('ScanNumEnd', frame.ScanNumEnd)
        s.setMetaValue('Precursor', frame.Precursor)
        
        # add spectrum to experiment 
        mzmlExp.addSpectrum(s)

############### Save File ######################
MzMLFile().store(args.fOut, mzmlExp)
