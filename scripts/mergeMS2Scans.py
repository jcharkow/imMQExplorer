# This script merges Spectra of the same precursor

#NOTE: the peak picking merger does not does peak pick on the IM axes (IM info outputted is garbage)


from pyopenms import *
import pandas as pd
import imMQExplorer.timsdata as td
import argparse


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

# merge peaks that are close to one another which are isolating the same precursor to boost S/N 
# NOTE: not exactly the same as the bruker peak picker however simmilar results. (e.g. m/z reported on combining is different
def mergePeaks(df, ppm = 20):
    df = df.sort_values(by='mz')
    
    ppm_arr = df['mz'] / ((1/ppm) * 1000000) #array of ppm values

    groups = (df['mz'].diff() >= ppm_arr).cumsum()
    wmean = lambda x: np.average(x, weights=df.loc[x.index, "i"])
    return df.groupby(groups).agg({'mz':wmean, 'i':sum})
    
    

# parse args
parser = argparse.ArgumentParser(description='fetch all MS2 windows of precursor in MzML format')
parser.add_argument('fInd', help='.d file in')
parser.add_argument('fInMzML', help='.mlML file in')
parser.add_argument('fOut', help='mzml File out')
parser.add_argument('--limit', help='number of precursor to include (useful in testing)', default=0, type=int)
parser.add_argument('--verbose', help='whether output is verbose', default="True")

args = parser.parse_args()
verbose = str2bool(args.verbose)

e = td.TimsData(args.fInd)

#### get table which tells us how many spectra there are for precursor take the precursor mass as the largest peak mz
q = e.conn.execute('select Precursor, count(Precursor), LargestPeakMz, IsolationWidth from PasefFrameMsMsInfo, Precursors where id = precursor group by Precursor order by Precursor')
rslt = pd.DataFrame(q.fetchall(), columns=['Precursor', 'NumPre', 'mz', 'IsolationWidth'])


# load the ondisc experiment (this takes a while)
exp = OnDiscMSExperiment()
exp.openFile(args.fInMzML)

if verbose:
    print("Experiment Loaded")


# experiment to store data
exp_out = MSExperiment()


ptr = 0 #index of the next spectrum to be loaded (0 means no spectrum has been loaded)  (0 indexed)

if args.limit == 0: #if 0 then this value is not set, sequence all precursors
    limit = len(rslt)
else:
    limit = args.limit


print("starting merge...")
for pre in rslt[:limit].itertuples():
    
    # load X frames (where X is the number of frames associated with the current precursor
    spec = [ exp.getSpectrum(i) for i in range(ptr, ptr + pre.NumPre) ]

    if len(spec) == 0:
        raise Exception("Length is 0")

    ptr += pre.NumPre # increment ptr

    if len(spec) == 1: #if only 1 frame just take that spectrum (don't have to merge anything)
        s = spec[0]

    else: #length greater then 1 have to merge

        # merge the array of spectrum into one spectrum
        

        # assert that all the spectrums precursors are equal to the expected precursor
        pres = np.array([ i.getMetaValue("Precursor") for i in spec ])
        try:
            assert(np.all(pres == pre.Precursor))
        except AssertionError:
            print(spec)
            print(pres)
            print(pre.Precursor)
            print(np.all(pres) == pre.Precursor)
            raise AssertionError()

        pres = np.array([ i.getPrecursors()[0].getMetaValue('Id') for i in spec ])
        assert(np.all(pres == pre.Precursor))



        # fetch the mz data and concatenate 
        s = MSSpectrum()
        peaks = [ i.get_peaks() for i in spec ]
        s.set_peaks(tuple(np.concatenate(peaks,  axis = 1)))

        # fetch and set IM values
        im_data = []
        for i in spec:
            try:
                im_data.append(i.getIntegerDataArrays()[0].get_data())
            except IndexError:
                assert(i.get_peaks()[0].size  == 0) #assert no output for getpeaks as well
                print("WARNING: Precursor {} has empty IM frame".format(pre.Precursor))
        im = IntegerDataArray()
        im.set_data(np.concatenate(im_data))
        im.setName('IonMobilityScan')
        s.setIntegerDataArrays([im])


    ## peak picking 
    peaks = pd.DataFrame(s.get_peaks()).T
    peaks = peaks.rename(columns={0:'mz', 1:'i'})
    peaks = mergePeaks(peaks)
    s.set_peaks(((peaks['mz'].values, peaks['i'].values)))


    # set the spectrum as having been centroided 
    s.setType(1)

    #### once have the scan data add the metadata
    s.setMetaValue('Precursor', pre.Precursor)

    # add precursor
    p = Precursor()
    p.setMZ(pre.mz) #set precursor mz as the largest peak mz detected by the instrument. 

    p.setIsolationWindowLowerOffset(pre.IsolationWidth / 2)
    p.setIsolationWindowUpperOffset(pre.IsolationWidth / 2)


    p.setMetaValue('Id', pre.Precursor)
    s.setPrecursors([p])


    # set to ms2 scan (not ms1 scan)
    s.setMSLevel(2)

    # add spectrum to experiment 
    exp_out.addSpectrum(s)


# save the experiment 
############### Save File ######################
MzMLFile().store(args.fOut, exp_out)
