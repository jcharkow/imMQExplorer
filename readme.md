The following contains source code is from the manuscript "Trapped Ion Mobility Spectrometry Reduces Spectral Complexity in Mass Spectrometry Based Workflow". 



Some possible analyses are included below:

1. Compute pariwise overlaps across peptide features
	- With IM
	`python featureDensity.py allPeptides.txt True mzrtim_ovlps.npy --rawFile 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768`
	- No IM
	`python featureDensity.py allPeptides.txt True mzrtim_ovlps.npy --rawFile 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768`

2. Compute overlaps between isolation windows and peptide features
	- With IM
		`python cofragmentationRates.py allPeptides.txt 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d True cofragIm.npy --rawFile 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768`
	
	- No IM
		`python cofragmentationRates.py allPeptides.txt 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d False cofragNoIm.npy --rawFile 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768`


3. Compute the Precursor Ion Fraction (PIF) of targeted peptides

	1. Fetch Isolation Window Feature Mapping 
		`python fetchIsolationWindowFeatureMapping.py allPeptides.txt 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d pasefMsmsScans.txt featureWindowMapping.pkl --rawFile 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768`

	2. Compute PIF 
		`python computePIF.py 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d featureWindowMapping.pkl pifComputations.pkl`
	
	3. Compute PIF Top Scans
		`python computePIFTopScans.py 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d featureWindowMapping.pkl pifComputationsTopScans.pkl`

4. Compute the Precursor Ion Fraction of untargeted peptides 
	1. Fetch Isolation Window Feature Mapping
		`python fetchUntargetedIsolationWindowFeatureMapping.py allPeptides.txt 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d/  untargetedFeatureWindowMapping.pkl --rawFile 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768`
	2. Compute PIF
		`python computePIF.py 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d untargetedFeatureWindowMapping.pkl untargetedPifComputations.pkl`


4. Predict Identification Rates using PIF
	`python pifToIdentifications.py pifComputations.pkl pifComputationsTopScans.pkl allPeptides.txt evidence.txt 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768` 

5. Predict Degree of cofragmentation in theoretical windows 
	1. Create Theoretical windows
		`python fetchUntargetedIsolationWindowFeatureMapping.py allPeptides.txt 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d untargetedFeatureWindowMapping.pkl --rawFile 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768`

	2. Compute PIF
		`python computePIF.py 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d UntargetedFeatureWindowMapping.pkl UntargetedPifComputations.pkl`


6. Perform Post-Acqusition Extraction of high quality TOF pushes
	- We found that this does not have a huge impact on identifications. The scripts are here to see our methodology. 
	1. Fetch MS2 scans from the .d file and merge (Control File)
		`python fetchMS2Scans.py 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d MS2-Frames.mzML`
		`python mergeMS2Scans.py 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d MS2-Frames.mzML MS2-Frames-Merged.mzML`
	2. Fetch (Based on post acquisition extraction) and merge
		`python fetchMS2ScansExtract.py 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d MS2-Frames-Extract.mzML pifComputationsTopScans.tsv 13`
		`python mergeMS2Scans.py 20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d MS2-Frames.mzML MS2-Frames-Merged-Extract.mzML`


Furthermore python scripts can be called directly in the imMQExplorer python functions



If used please cite the following manuscript:
`Charkow and Rost, “Trapped Ion Mobility Spectrometry Reduces Spectral Complexity in Mass Spectrometry Based Workflow.”`
