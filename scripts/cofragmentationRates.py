# This script computes all of the pairwise overlaps with and without IM separation.


import argparse
from  imMQExplorer.featureExperiment import *
from imMQExplorer.commonFunctions import *
import pandas as pd

parser = argparse.ArgumentParser(description="computes the pariwise overlaps between peptide features and a given set of frames. Frames can either be supplied by a .d experiment or a pandas df with columns <rtStart, rtEnd, mzStart, mzEnd, scanNumBegin (optional), scanNumEnd (optional)> Indices outputted correspond to row # of allPeptides.txt features")

parser.add_argument('fIn', help='allPeptides.txt input from MQ')
parser.add_argument('windows', help='frames supplied as a pandas df (with columsn mzStart, mzEnd, rtStart, rtEnd, scanNumBegin (optional), scanNumEnd (optional)) or .d experiment')
parser.add_argument('--rawFile', help='if set only consider features from this experiment', default=None)
parser.add_argument('--ionFWHM', help='IM feature boundaries set by FWHM', default=True)
parser.add_argument('--rtFWHM', help='RT feature boundaries set by FWHM', default=True)
parser.add_argument("im", help='Should IM dimension be included when computting overlaps, type should be bool')
parser.add_argument("fOut", help="name of .pkl dataframe of outputted pairwise overlaps")

args = parser.parse_args()


# convert IM to bool
compute_im = str2bool(args.im)


# load allpeptides.txt
exp = Experiment(args.windows)
pep = Precursor(args.fIn, rawFile=args.rawFile, ionFWHM= args.ionFWHM, retentionFWHM = args.rtFWHM)

# compute overlaps
rslt = pep.overlap_against_window(exp, im=compute_im)

#save as pandas DF
pd.to_pickle(pd.DataFrame(rslt, columns=['Feature', 'FrameIndex']), args.fOut)


