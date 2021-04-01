# provide mapping between isolation windows from a TIMS experiment and peptide features

import argparse
import pandas as pd
import imMQExplorer.timsdata as td
import numpy as np

parser = argparse.ArgumentParser(description="provide mapping between isolation windows from a TIMS experiment and peptide features, useful for computing PIF.")
parser.add_argument('fIn', help='allPeptides.txt output from MQ')
parser.add_argument('exp', help='.d experiment file')
parser.add_argument('fOut', help='Name of metaData file')
parser.add_argument('--mz', help='MS2 mz window', nargs = '?', default='2', const = '2', type=int)
parser.add_argument('--im', help='MS2 im window', nargs = '?', default='26', const = '26', type=int)
parser.add_argument('--rawFile', help='rawFile name as specified in allPeptides.txt (for the corresponding experiment)')

args = parser.parse_args()
mz_width = args.mz
im_width = args.im

pep = pd.read_csv(args.fIn, '\t')
pep = pep[pep['Raw file'] == args.rawFile]

#filter to untargeted pep
pep = pep[pep['Pasef MS/MS IDs'].isna()]

ms1_scans = np.empty(len(pep))

exp = td.TimsData(args.exp)

##### Compute pre-to-frames ####
for row_num, (_,row) in enumerate(pep.iterrows()):
    if row_num % 1000 == 0:
        print("{} of {} complete".format(row_num, len(pep)))
    query = exp.conn.execute("select Id, min(abs(? - Frames.Time)) from Frames where MsMsType == 0", [row['Retention time']])
    rslt = query.fetchall()
    #only one window should be fetched
    assert(len(rslt) == 1)

    #from the SQL output record the target MS1 frame
    ms1_scans[row_num] = rslt[0][0]

#store precursor and ms1scan information
pre_to_frame = pd.DataFrame(np.column_stack([pep.index, ms1_scans]), columns=['pre','Ms1Scan'])

# convert to int (if they were float)
pre_to_frame['pre'] = pre_to_frame['pre'].astype(int)
pre_to_frame['Ms1Scan'] = pre_to_frame['Ms1Scan'].astype(int)

print("add columns ...")
####add columns to on pseudo scan num start and end####
pep['ScanNumBegin'] = np.floor(pep['Ion mobility index'] - (im_width/2)).astype(int) #these can only be integer values, overestimate the range by rounding down
pep['ScanNumEnd'] = np.ceil(pep['Ion mobility index'] + (im_width/2)).astype(int) #ceil the end to overestimate the range of IM (cannot have half a scan)

##### add columns on MQ Im start and end ####
pep['MQImStart'] = pep['Ion mobility index'] - (pep['Ion mobility index length (FWHM)'] / 2)
pep['MQImEnd'] = pep['Ion mobility index'] + (pep['Ion mobility index length (FWHM)'] / 2)

#add columns on theoretical mz window
pep['IsolationMzBegin'] = pep['m/z'] - (mz_width / 2)
pep['IsolationMzEnd'] = pep['m/z'] + (mz_width/2)
#merge the important data from both the tables
meta_data = pd.merge(pep[['m/z', 'Charge','Number of isotopic peaks', 'ScanNumBegin', 'ScanNumEnd', 'IsolationMzBegin', 'IsolationMzEnd', 'MQImStart', 'MQImEnd']], pre_to_frame, left_on=pep.index, right_on='pre')

meta_data['Feature'] = meta_data['pre'].astype(int)
meta_data = meta_data.drop(columns='pre')

#sort the values by ms1 scan
meta_data = meta_data.sort_values(by="Ms1Scan")
meta_data = meta_data.rename(columns={"Ms1Scan":"MS1Frame"})

#rename columns for access with itertuples
meta_data = meta_data.rename(columns={"m/z":"Mz", "Number of isotopic peaks":"NumberOfIsotopicPeaks"})

#reindex
meta_data.index = np.arange(0, len(meta_data))

#save
meta_data.to_pickle(args.fOut)
