spectrum_max = 90000 
spectrum_min = 40000
vanadium_raw_file = "C:/Users/spu92482/Desktop/vanadium/WISH00018589.RAW"
sample_raw_file = "C:/Users/spu92482/Desktop/vanadium/WISH00017986.RAW"

#load workspaces
LoadRaw(Filename=vanadium_raw_file, OutputWorkspace="vanadium",SpectrumMin=spectrum_min,SpectrumMax=spectrum_max);
LoadRaw(Filename=sample_raw_file, OutputWorkspace="sample",SpectrumMin=spectrum_min,SpectrumMax=spectrum_max);
#smooth neighbours on both input workspaces to get rid of noise across spectra
SmoothNeighbours(InputWorkspace="vanadium", OutputWorkspace="smoothed_vanadium",ProvideRadius=False,WeightedSum=False,NumberOfNeighbours=35)
SmoothNeighbours(InputWorkspace="sample", OutputWorkspace="smoothed_sample",ProvideRadius=False,WeightedSum=False,NumberOfNeighbours=35)
#clean-up original workspaces
mantid.deleteWorkspace('vanadium')
mantid.deleteWorkspace('sample')
#smooth data to get rid of noise accross bins in vanadium workspace.
SmoothData(InputWorkspace="smoothed_vanadium", OutputWorkspace="smoothed_vanadium",NPoints=50);
#divide to normalise by vanadium
Divide(LHSWorkspace="smoothed_sample", RHSWorkspace="smoothed_vanadium",OutputWorkspace="normalised_sample");
#clear out everything apart from the normalised by vanadium sample workspace
mantid.deleteWorkspace('smoothed_vanadium')
mantid.deleteWorkspace('smoothed_sample')
#Replace NaN and infinite values, also threshold big numbers
ReplaceSpecialValues(InputWorkspace='normalised_sample',OutputWorkspace='normalised_sample',InfinityValue='800',BigNumberThreshold='1000',BigNumberValue='1000')
#convert to md workspace

#---------
#Note that with your Peak finding algorithm you would be able to run the following
#1 FindPeaks -- Find your peaks using your detector space peak finding routine
#2 CalculateUMatrix 
#3 CopySample -- Copy the sample off the peaks workspace and onto your original MatrixWorkspace
#4 ConvertToDiffractionMDWorkspace -- using Q HKL
#----------

#Split into and split threshold parameters have been set-up to favour peak finding in Q-space.
ConvertToDiffractionMDWorkspace(InputWorkspace='normalised_sample',SplitInto=2,SplitThreshold=100,OutputWorkspace='md', Extents='-4,4,-4,4,-4,4')





