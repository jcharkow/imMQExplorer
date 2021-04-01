# This script invokes the may-2020 DDA PIF script to compute the top 50% of scans.

from imMQExplorer.pif import PIFOverlapCalculations
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='calculates optimized PIF (PIF of the top X scans)')

parser.add_argument('exp', help='.d experiment file')
parser.add_argument('metaData', help='.pkl metaData file')
parser.add_argument('fOut', help='name of resulting .pkl file to produce')
parser.add_argument('--numScans', type=int, help='number of top scans to extract (may not always be this number if IM signal does not span whole window or there is a tie for the number of top scans) default = 13')


args = parser.parse_args()

if not 'args.numScans' in locals():
    numScans = 13 
else:
    numScans = args.numScans


metaData = pd.read_pickle(args.metaData)
pifObj = PIFOverlapCalculations(args.exp, metaData)

#do the pif calculations 
pifObj.calcPIFTopScans(numScans)

pifObj.savePIFMeta(args.fOut)

