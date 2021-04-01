# provide mapping between isolation windows from a TIMS experiment and peptide features

import argparse
import pandas as pd
import imMQExplorer.timsdata as td


parser = argparse.ArgumentParser(description="provide mapping between isolation windows from a TIMS experiment and peptide features, useful for computing PIF.")
parser.add_argument('fIn', help='allPeptides.txt output from MQ')
parser.add_argument('--rawFile', help='rawFile name as specified in allPeptides.txt (for the corresponding experiment)')
parser.add_argument('exp', help='.d experiment file')
parser.add_argument('pasef', help='pasefMSMSScans.txt file from MQ output')
parser.add_argument("fOut", help="name of .pkl dataframe of outputted mapping")

args = parser.parse_args()


pep = pd.read_csv(args.fIn, '\t')
pep = pep[pep['Raw file'] == args.rawFile]

pasef = pd.read_csv(args.pasef, '\t')
pasef = pasef[pasef['Raw file'] == args.rawFile]

pasefMSMSIDs = pep['Pasef MS/MS IDs'].dropna()
pre_to_frame = []
for pre, frames in pasefMSMSIDs.iteritems():
    frame_lst = frames.split(';')
    for i in frame_lst:
        pre_to_frame.append((pre, int(i)))
pre_to_frame = pd.DataFrame(pre_to_frame, columns = ['Feature', 'Index']) #this is the peptide feature to pasefMSMSIndex


#get mapping of precursor to MS1 frame
exp = td.TimsData(args.exp)
query = exp.conn.execute("select id as precursor, parent from Precursors")
precursorToMs1 = pd.DataFrame(query.fetchall(), columns=['Precursor', 'MS1Frame'])

metaData = (pre_to_frame
            .merge(pasef)
            .merge(pep[['m/z','Charge', 'Number of isotopic peaks', 'Ion mobility index length', 'Ion mobility index' ]], left_on='Feature', right_index=True) #use the IM index (not FWHM to be extra conservative
            .merge(precursorToMs1))

#add the isolation window start and end
metaData['IsolationMzBegin'] = metaData['IsolationMz'] - (metaData['IsolationWidth'] / 2)
metaData['IsolationMzEnd'] = metaData['IsolationMz'] + (metaData['IsolationWidth'] / 2)
metaData['MQImStart'] = metaData['Ion mobility index'] - (metaData['Ion mobility index length'] / 2)
metaData['MQImEnd'] = metaData['Ion mobility index'] + (metaData['Ion mobility index length'] / 2)


#rename column entries 
metaDataFinal = metaData.rename(columns={'Number of isotopic peaks':'NumberOfIsotopicPeaks', 'm/z':'Mz', 'Frame':'MS2Frame'}).drop(['CollisionEnergy', 'Raw file', 'IsolationMz', 'IsolationWidth', 'Ion mobility index length', 'Ion mobility index'], axis = 1)


metaDataFinal.to_pickle(args.fOut)

