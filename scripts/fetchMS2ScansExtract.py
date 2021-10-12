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
parser.add_argument('pif', help='pif .tsv file for determining top scans (.pkl protocol was not working)')
parser.add_argument('topScansCutoff', type=int, help='cutoff for filtering scans, if do not have this many "good quality scans" then dont filter')
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


def string_to_array(s):
    if isinstance(s, float):
        return []
    s_array = s.split(" ")
    tmp_arr = []
    for str_split in s_array:
        tmp_str = ""
        for char in str_split:
            if char.isdigit():
                tmp_str += char
        if tmp_str != "":
            tmp_str = int(tmp_str)
            tmp_arr.append(tmp_str)
    return tmp_arr


# NOTE: not all precursors (targeted by instrument) were detected by MQ in the MS1 spectra, those precursors not excluded but not optimized (Only 261367 / 306385 precursors detected)
#get PIF results
pif = pd.read_csv(args.pif, '\t') # pickle protocol was not working so changed to csv

#PIF processing. Some "Pre" selected by the instrument have more than one feature (Precursor). select the "best" one.
pif_processed = pif.sort_values(by='ImPIFTop13Scans', ascending=False).groupby("Precursor").head(1)

#merge all tables 
pif_pasef_pre = pd.merge(pif_processed[['topScansIdx', 'Precursor', 'topScansLen']], pasef_pre, left_on=['Precursor'], right_on=['Precursor'], how='right')

### new process topScansIdx to be array ###
pif_pasef_pre['topScansIdx'] = pif_pasef_pre['topScansIdx'].apply(string_to_array)


############## Fetch MS2 Frames ###################
curPrecursor = (0,None)
mzmlExp = MSExperiment()
for frame in pif_pasef_pre.itertuples():
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

        # filter the MS2 data based on the topScans idx (if possible)
        if frame.topScansLen >= args.topScansCutoff: # only filter scans that have a reported top 13 scans, if less than top 13 reported do not filter. 

            # turn numbers from relative to absolute scan numbers
            absScans = np.array(frame.topScansIdx) + frame.ScanNumBegin

            # filter "raw" to only those entries in target IM scans
            raw = raw[raw['scan'].isin(absScans)]


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
