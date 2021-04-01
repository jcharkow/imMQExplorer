# This script computes all of the pairwise overlaps with and without IM separation.


import argparse
from  imMQExplorer.featureExperiment import *
from imMQExplorer.commonFunctions import *
import pandas as pd

parser = argparse.ArgumentParser(description="computes the pariwise overlaps of peptide features. Indices outputted correspond to row # of allPeptides.txt features")

parser.add_argument('fIn', help='allPeptides.txt input from MQ')
parser.add_argument('--rawFile', help='if set only consider features from this experiment', default=None)
parser.add_argument('--ionFWHM', help='IM feature boundaries set by FWHM', default=True)
parser.add_argument('--rtFWHM', help='RT feature boundaries set by FWHM', default=True)
parser.add_argument("im", help='Should IM dimension be included when computting overlaps, type should be bool')
parser.add_argument("fOut", help="name of .pkl pandas dataframe of outputted pairwise overlaps")

args = parser.parse_args()


# convert IM to bool
compute_im = str2bool(args.im)


# load allpeptides.txt
pep = Precursor(args.fIn, rawFile=args.rawFile, ionFWHM= args.ionFWHM, retentionFWHM = args.rtFWHM)

# compute overlaps
rslt = pep.overlap(im=compute_im)

#save 
(pd.DataFrame(rslt, columns=['feature1','feature2'])).to_pickle(args.fOut)


