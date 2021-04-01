# This script converts a list of features from MQ output allpeptides.txt to theoretical isolation windows. 


import argparse
from imMQExplorer.makeTheoreticalIsolationWindows import TheoreticalIsolationWindows
import pandas as pd

parser = argparse.ArgumentParser(description="computes theoretical isolation windows based on a set of features. isolation windows are supplied as a range in mz, rt and im.") 

parser.add_argument('fIn', help='allPeptides.txt input from MQ')
parser.add_argument('--rawFile', help='if set only consider features from this experiment', default=None)
parser.add_argument('acquisition', help='acquisition mode, can be dia or dda')
parser.add_argument("fOut", help="name of .pkl dataframe of outputted theoretical isolation windows")
parser.add_argument("--mzWidth", help='width for mz window', type=float)
parser.add_argument("--imWidth", help='width for im window', type=float)
parser.add_argument("--rtWidth", help='width in retention time', type=float)
parser.add_argument("--imRangeStart", help='starting scan for imRange', type=int)
parser.add_argument("--imRangeEnd", help='ending scan for imRange', type=int)


args = parser.parse_args()

win = TheoreticalIsolationWindows(args.fIn, rawFile=args.rawFile)

## set params
if 'args.mzWidth' in locals():
    win.setParams(self, mzWidth, args.mzWidth)
if 'args.imWidth' in locals():
    win.setParams(self, imWidth, args.imWidth)
if 'args.rtWidth' in locals():
    win.setParams(self, rtWidth, args.rtWidth)
if 'args.imRangeStart' in locals() and 'args.imRangeEnd' in locals():
    win.setParams(self, imWidth, (args.imWidthStart, args.imWidthEnd))


# compute theoretical windows
theoWindows = win.createTheoreticalWindows(args.acquisition)


theoWindows.to_pickle(args.fOut)

