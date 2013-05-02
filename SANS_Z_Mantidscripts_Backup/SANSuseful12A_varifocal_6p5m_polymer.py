path="c:/mantidinstall/data/"
# ===========================the NEW scripting commands =================================

#Make the reduction module available
from ISISCommandInterface import *
from math import *

# Set a default data path where we look for raw data
DataPath("Y:/cycle_11_3/")
DataPath("Z:/processed/SANS2d/")
# Set a user path where we look for the mask file
UserPath("Z:/MASKS/")


# Set instrument to LOQ
#LOQ()
#Set reduction to 1D (note that if this is left out, 1D is the default)
#Set1D()

# Read a mask file
#MaskFile('MASKLOQ_102I.txt')
# Assign run numbers (.nxs for nexus)
#AssignSample('59890.nxs')
#TransmissionSample('59891.nxs', '59889.nxs')
#AssignCan('59888.nxs')
#TransmissionCan('59889.nxs', '59889.nxs')
# WavRangeReduction runs the reduction for the specfied wavelength range
#wav1 = "2.2"
#wav2 = "10.0"
#reduced = WavRangeReduction(wav1, wav2, DefaultTrans)
#wav1 = "3.0"
#wav2 = "4.0"
#reduced = WavRangeReduction(wav1, wav2, DefaultTrans,no_clean=True)
#wav1 = "4.0"
#wav2 = "5.0"
#reduced = WavRangeReduction(wav1, wav2, DefaultTrans,no_clean=True)


# Set instrument to SANS2D, FOCUS loops 
SANS2D()
#Set reduction to 1D (note that if this is left out, 1D is the default)
Set1D()
# Read a mask file.   NOTE same direct beam file needs to be mentioned below !
maskfilename="z:/Masks/masksans2d_113cc_12m_M1_RKH_v13.txt"
# pull in the direct beam that we need to modify MAKE SURE IT IS THE SAME AS IN THE MASK FILE !
LoadRKH("z:/Masks/DIRECTM1_10302_12m_14Oct11_v13.dat","directbeam")
rsam="10302-add"
rcan="10300-add"
tsam="10293"
tdb="10291"
tcan="10291"
wlist=[1.75,2.125,2.5,2.875,3.25,3.625,4.0,4.5,5.0,6.0,8.0,9.0,10.0,11.0,12.5]
#wlist=[1.75,1.975,2.125,2.375,2.625,2.875,3.125,3.625,4.125,4.625,5.125,5.625,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.5]

MaskFile(maskfilename)
#     WON'T WORK   as string ???             qlimits="0.002,1.0,0.08"
LimitsQ(0.002,.15,0.08,'LOG')
W1=1.75
Q1A=.06
Q1B=.2
W2=12.5
Q2A=0.004
Q2B=0.02
# Change default limits if we want to
#LimitsR(50.,170.)
#   Here the arguments are min, max, step, step type
#LimitsWav(4.,8.,0.125, 'LIN')
#LimitsQXY(0, 0.1, 0.002, 'LIN')
# command to switch banks in from Python is Detector('bank-name'), i.e. Detector('front-detector').
#Detector('front-detector')
# Gravity(True) or Gravity(False) that will switch accounting for gravity on and off 
# in the next call to WavRangeReduction.
#
AssignSample(rsam+'.nxs')
AssignCan(rcan+'.nxs')
# this loads the workspaces but does not do the calc, so you have to do one data reduction to get it calculated, then
# modify it and do the data reduction again
TransFit('Log',lambdamin=1.75,lambdamax=12.5)
TransmissionSample(tsam+'.nxs', tdb+'.nxs')
TransmissionCan(tcan+'.nxs', tdb+'.nxs')
# Update the centre coordinates
# Arguments are rmin,rmax, niterations
#FindBeamCentre(50., 170., 2)
# WavRangeReduction runs the reduction for the specfied wavelength range
# The final argument can either be DefaultTrans or CalcTrans where
#  (a) DefaultTrans calculates the transmission for the whole range specified by L/WAV 
#      and then crops it to the current range and
#  (b) CalcTrans calculates the transmission for the specifed range
#                             [ DefaultTrans seems only a logical, which is True ]
# The returned value is the name of the fully reduced workspace
nloops=len(wlist)-1
reducedAll = WavRangeReduction(wlist[0], wlist[nloops], False)
#Integration(reducedAll,"integral",Q1,Q2)
#norm = mtd["integral"].readY(0)[0]
#normE = mtd["integral"].readE(0)[0]
#print norm,normE
# BEWARE to get the gui to give correct results for SANS2d after running from python here you may need
# to switch to LOQ, load a LOQ mask, switch back to SANS2d, then load a SANS2d mask!
# 12/04/11 RKH checked that gui & python give same results
#
plt = plotSpectrum(reducedAll,0)
#
# reduce data
outlist=[]
scaleX=[]
scaleY=[]
scaleE=[]
#  must go to unity at extremes of ranges, these get modified again below
scaleX.append(0.5)
scaleY.append(1.0*(1.125-0.5))
scaleE.append(0.0)
scaleX.append(1.125)
scaleY.append(1.0*(wlist[0]-1.125))
scaleE.append(0.0)
#
for i in range(nloops):
	SANS2D()
#Set reduction to 1D (note that if this is left out, 1D is the default)
	Set1D()
	MaskFile(maskfilename)
	LimitsQ(0.002,0.15,0.08,'LOG')
# Assign run numbers (.nxs for nexus)
	AssignSample(rsam+'.nxs')
	AssignCan(rcan+'.nxs')
	wav1 = wlist[i]
	wav2 = wlist[i+1]
	TransFit('Log',lambdamin=1.75,lambdamax=12.5)
	TransmissionSample(tsam+'.nxs', tdb+'.nxs')
	TransmissionCan(tcan+'.nxs', tdb+'.nxs')
#  reduced = WavRangeReduction(wav1, wav2, DefaultTrans, no_clean=True)
	reduced = WavRangeReduction(wav1, wav2, False)
	outlist.append(reduced)
	mergePlots(plt,plotSpectrum(reduced,0))
# set up varying Q range to integrate over
	QA=Q2A+(1./wav2-1./W2)*(Q1A-Q2A)/(1./W1-1./W2)
	QB=Q2B+(1./wav1-1./W2)*(Q1B-Q2B)/(1./W1-1./W2)
#	print W1,Q1A,Q1B,W2,Q2A,Q2B
#	print QA,QB
# integrate full wavelength data over chosen Q range
	Integration(reducedAll,"integral",QA,QB)
	norm = mtd["integral"].readY(0)[0]
	normE = mtd["integral"].readE(0)[0]
#  now integrate this loop 
	Integration(reduced,"integral",QA,QB)
	scaleX.append(wav1)
	int=mtd["integral"].readY(0)[0]
	print wav1,wav2,QA,QB,norm,normE,int
	scaleY.append((int/norm)*(wav2-wav1))
	scaleE.append(sqrt( ((mtd["integral"].readE(0)[0])/int)**2 + (normE/norm)**2 )*(wav2-wav1))
scaleX.append(wav2)
#  must go to unity or same as last value at extremes of ranges 
scaleX.append(19.0)
scaleY.append((int/norm)*(19.0-wav2))
scaleE.append(sqrt( ((mtd["integral"].readE(0)[0])/int)**2 + (normE/norm)**2 ))
scaleX.append(20.0)
scaleY.append((int/norm)*(20.0-19.0))
scaleE.append(sqrt( ((mtd["integral"].readE(0)[0])/int)**2 + (normE/norm)**2 ))
# now change the two dummy points at the start 
scaleY[0]=scaleY[2]*(1.125-0.5)/(wlist[1]-wlist[0])
scaleY[1]=scaleY[2]*(wlist[0]-1.125)/(wlist[1]-wlist[0])
# ======================================================================================
# need sort some errors for this
CreateWorkspace("scale",scaleX,scaleY,scaleE,UnitX="Wavelength")
# small problem, LoadRKH brings in a distribution,, InterpolatingRebin needs a histogram
# ConvertToHistogram leaves Y values unchanged, just splits the x values, adding one extra x 
ConvertToHistogram("directbeam","directbeam_hist")
# interpolating rebin expects a proper histogram wtih counts and time,
# it seems to need a big over-run at the ends of the range!
InterpolatingRebin("scale","scale2","0.8,0.1,18.3")
# ConvertToDistribution only divides by bin width, leave original X values !
ConvertToDistribution("scale")
ConvertToDistribution("scale2")
#CropWorkspace("directbeam_hist","directbeam_hist_crop",XMin=1.0,XMax=15.0)
Multiply("directbeam_hist","scale2","directbeam_new_hist")
SaveRKH("directbeam_new_hist","z:/Masks/DIRECTM1_10302_12m_14Oct11_v14.dat")
#====================================================================================================
#
#  add results of loop to original plot
for i in range(nloops):
    mergePlots(plt, plotSpectrum(outlist[i],0))

# ===========================END of the NEW scripting commands =================================

# thus far SANS2d transmission files are same as for sans runs, but will change before too long
#SANS2d_trans_inst_def_file="c:/mantidinstall/instrument/SANS2d_definition.xml"
#path="//isis/inst$/NDXSANS2D/Instrument/data/cycle_09_2/"
#path = '/home/dmn58364/Desktop/SANS2D_Data/'
# get whole raw file once only
#instr_dir = '/home/dmn58364/Mantid/trunk/Test/Instrument/'
#instr_dir = '/home/dmn58364/Mantid/trunk/Test/Instrument/'
#
#MaskStringRear = 'h0,h190>h191,v0,v191,v64>v127+h64>h191,v128>v191+h0>h127'
#MaskStringRear = 'h0,h190>h191,v0,v191'
#MaskStringRear = 'h0,h190>h191,v0,v191,h131>h140'
#detbank = "rear-detector"
#detbank = "front-detector"
#correction_type = '2D'
#bank='_mask'
#
# pick spectrum for centre, enter number below
wksp = mantid.getMatrixWorkspace("3316_sans_nxs")
instrument = wksp.getInstrument()
det = wksp.getDetector(11000)
# then check the coords for the shift 
out=str(det.getPos())
# NEED to read X&Y coords from string "out" here 
print "det is  "+out 
#
#
# D_H or D_V  sums, need get these into a gui
#
wksp = '3324_sans_nxs'
wksp = '6116'
CropWorkspace(wksp, OutputWorkspace = 'rear', StartWorkspaceIndex= 8, EndWorkspaceIndex=(192*192 + 7))
CropWorkspace(wksp, OutputWorkspace = 'front', StartWorkspaceIndex= (192*192+8), EndWorkspaceIndex=(192*192*2 + 7))
# help(SumRowColumn)
# SumRowColumn(InputWorkspace, OutputWorkspace, Orientation, XMin=-9999, XMax=-9999, HVMin=-9999, HVMax=-9999)
# Orientation is "D_H" or "D_V" at present, Xmin/max is for TIMES, HVmin/max for rows or cols 
#
# was SumRowColumn('rear', 'rear_H', 'D_H', XMin=-9999, XMax=-9999, HVMin=0, HVMax=9999)
SumRowColumn('rear', 'rear_H', 'D_H')
SumRowColumn('rear', 'rear_V', 'D_V')
SumRowColumn('front', 'front_H', 'D_H')
SumRowColumn('front', 'front_V', 'D_V')
#
# series of D_V or D_H to help find issues around the beam stop
for i in range (50,65,1):
	SumRowColumn('rear', 'rear_H'+str(i), 'D_H', HVMin=i, HVMax=i)
for i in range (60,74,1):
	SumRowColumn('rear', 'rear_V'+str(i), 'D_V', HVMin=i, HVMax=i)
#
SumRowColumn('rear', 'rear_V130', 'D_V', HVMin=130, HVMax=130)
SumRowColumn('rear', 'rear_H60', 'D_H', HVMin=60, HVMax=60)
SumRowColumn('front', 'front_H61', 'D_H', HVMin=61, HVMax=61)
#
# ================================================================================================================================
# 13/12/11  compare D_T over same radius range for different runs.   
import SANSUtility_RKH
# for SANS2d run data reduction FIRST to get detector in correct place!  then rename nxs wksp
# the starting workspace is preserved here, delete the _Q ones as they are LARGE.
wksp="9672_sans_nxs"
# select radius range in metres
rmin=0.041
rmax=100.0
CropWorkspace(wksp, OutputWorkspace = wksp+"_Q", StartWorkspaceIndex=(8+1-1), EndWorkspaceIndex=str(8+192*192-1))
NormaliseByCurrent(wksp+"_Q",wksp+"_Q")
print 'rmin = '+str(rmin)
print 'rmax = '+str(rmax)
if rmin > 0.0: 
		SANSUtility_RKH.MaskInsideCylinder(wksp+"_Q", rmin)
if rmax > 0.0 and rmax < 10.0 :
		SANSUtility_RKH.MaskOutsideCylinder(wksp+"_Q", rmax)
# 
ConvertUnits(wksp+"_Q",wksp+"_Q","MomentumTransfer",Emode="Elastic",AlignBins=0)
# must rebin here before final sum
Rebin(wksp+"_Q",wksp+"_Q","0.001,-0.08,2.0")
#ConvertUnits(wksp,wksp,"Wavelength",AlignBins=1)
SumSpectra(wksp+"_Q",wksp+"_Q_sum")
#
# ================================================================================================================================
#
# D_T  sums, need get these into a gui
#
wksp="SANS2D00001674"

CropWorkspace(wksp, OutputWorkspace = "D_T", StartWorkspaceIndex=(8+1-1), EndWorkspaceIndex=str(8+192*192-1))


path="c:/mantidinstall/data/"
run="10000816"
wksp=run+"_sans_nexus"
LoadNexus(path+"SANS2d"+run+".nxs",wksp)
wksp="10850_sans_nxs"
CropWorkspace(wksp, OutputWorkspace = "D_T", StartWorkspaceIndex=(8+1-1), EndWorkspaceIndex=str(8+192*192-1))
# sum all into a single spectrum  this does the COLETTE, DISPLAY TIME or D T
SumSpectra("D_T","D_T")
plot1=plotSpectrum("D_T",0)
Integration("D_T","D_T_integral",10,100000)	
#
wksp="59890_sans_nxs"
CropWorkspace(wksp, OutputWorkspace = "D_Tf", StartWorkspaceIndex=(8+192*192-1), EndWorkspaceIndex=str(8+192*192*2-1))
# sum all into a single spectrum  this does the COLETTE, DISPLAY TIME or D T
SumSpectra("D_Tf","D_Tf")
Integration("D_Tf","D_Tf_integral",10,100000)	
#  LOQ
CropWorkspace(wksp, OutputWorkspace = "D_T_LOQ", StartWorkspaceIndex=(3+1-1), EndWorkspaceIndex=str(3+128*128-1))
# sum all into a single spectrum  this does the COLETTE, DISPLAY TIME or D T
SumSpectra("D_T_LOQ","D_T_LOQ")
Integration("D_T_LOQ","D_LOQ_integral",3500,43500)	
#
#
import sys
sys.path.append("C:/Mantidinstall/scripts/SANS")
print sys.path
import SANSUtility
# D_T for nine blocks when dae playing up  in cycle 09/1 or 09/2
wksp="SANS2D00001674"
wksp="5976_sans_nxs"
dimdet=192
spec1=8+1+192*192
tag="frt"
#spec1=8+1
#tag="rear"
#print speclist
i=0
namelist=["A TL","D TC","G TR","B ML","E MC","H MR","C BL","F BC","I BR"]
masklist=["h0>h127,v64>v191","h0>h127,v0>v63,v128>v191","h0>h127,v0>v127","h0>h63,h128>h191,v64>v191","h0>h63,h128>h191,v0>v63,v128>v191","h0>h63,h128>h191,v0>v127","h63>h191,v64>v191","h63>h191,v0>v63,v128>v191","h63>h191,v0>v127"]
for name in namelist:
	print "i= "+str(i)
	# guess "orientation"
	list = SANSUtility.ConvertToSpecList(masklist[i], spec1, dimdet,'0')
	CropWorkspace(wksp, OutputWorkspace = "D_T "+name+tag, StartWorkspaceIndex=(spec1-1), EndWorkspaceIndex=str(spec1+192*192-2))
	SANSUtility.MaskBySpecNumber("D_T "+name+tag, list)
	SumSpectra("D_T "+name+tag,"D_T "+name+tag)
	i=i+1
print "done"
for wext in namelist:
	wsp="D_T "+wext
	if wext==namelist[0]:
		plot2=plotSpectrum(wsp+tag,0)
		layer=plot2.activeLayer()
		layer.setTitle("D_T modules")
	else:
		mergePlots(plot2,plotSpectrum(wsp+tag,0))
	layer.setCurveTitle(i,name)
#
#
#  ================================       VERY USEFUL =============================================================	
#
#17/12/09  USING NEW  CALL TO WavRangeReduction 
# Looping over other ranges
#  YOU MUST do the full range in gui first  (unless use scripting to set everything up first )
#first= WavRangeReduction(wlist[0], wlist[len(wlist)-1], DefaultTrans)
#  ( BEWARE - bug at present requires each limit to be INTEGERS, else get a divide error with transmissions )
#wlist=[2.0,4.0,6.0,8.0,10.0,11.5,12.75]
wlist=[2.2,3.0,4.0,6.0,8.0,10.0]
wlist=[2.0,4.0,6.0,8.0,10.0,12.0,14.0]
wlist=[2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0]
wlist=[4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0]
#wlist=[4.6,6.1,7.1,8.1,9.1,10.1,11.6,12.85]
wlist=[5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0]
# reduce data
outlist=[]
for i in range(len(wlist)-1):
    wav1 = wlist[i]
    wav2 = wlist[i+1]
    reduced = WavRangeReduction(wav1, wav2, DefaultTrans, no_clean=True)
#    reduced = WavRangeReduction(wav1, wav2, DefaultTrans)
    outlist.append(reduced)
#
#  PLOT the data
#  .split will generate a tuple  e,g, from  10002694front_11.5_12.75  get ['10002694front', '11.5', '12.75']
file=reduced.split('_')
# reassmble name of full range workspace as run from gui  
# BUG here if we have not actually run the full range of WLIST
first = file[0]+'_'+str(wlist[0])+'_'+str(wlist[len(wlist)-1])
# need to see if "first" exists, else then use outlist[0] ?
print first
plt = plotSpectrum(first, 0)
layer=plt.activeLayer()
layer.setTitle(file[0])
print i
#print outlist[6]
#
for i in range(len(outlist)):
    wksp= outlist[i] 
    mergePlots(plt, plotSpectrum(wksp,0))
#
# SAVEthe workspaces to file   (generate workspace name "first" above !)
filename="c:/mantidinstall/genie/"+file[0]+'.Q'
SaveRKH(first,filename,Append='0')
# now save all the partial ranges to file
for i in range(len(outlist)):
	SaveRKH(outlist[i],filename,Append='1')
print "Written to "+filename
#
#
# ADDING .raw files to make processed nexus file, NOTE add 1000 to front of run number
# NOTE until we can read the logs from the nexus file you also need to copy and rename a .log file for each, into the same directory.
# have seen issues of not subtracting bkg when SMALL changes in encoder positions happen, in  which case for now copy the can log file over the sample one...
#
#  the 999 version does three figure run numbers (need a string routine to put corrrect number of zeroes on the front !)
import SANSadd999
path="//isis/inst$/NDXSANS2D/Instrument/data/cycle_09_2/"
path="N:/cycle_09_1/"
SANSadd.add_runs(path, [541,567])
SaveNexusProcessed("added","i:/SANS2d10000567.nxs")

# rotating SDS 6m data
SANSadd999.add_runs(path, [769,776,782])
SaveNexusProcessed("added","c:/mantidinstall/data/SANS2d10000782.nxs")


# 17/03/10 new version of adding files, writes to U:/user/processed which should be in search list, see notes in SANSadd2.py
# to find SANSadd2 may need 
import sys
sys.path.append("U:/Mantidscripts")
print sys.path
import SANSadd2

path="O:/cycle_09_5/"
path="I:/"
# 6m data
SANSadd2.add_runs(path, [3333,3339])
#SaveNexusProcessed("added","i:/SANS2d00003362-add.nxs")


# new version, simpler, writes to z:/processed/sans2d/
# beware, not passing these params through ????
path="X:/"
import SANSadd2
SANSadd2.add_runs(path,[11971,11980,11984])


# MUST DO THIS BEFORE CREATING     D_T     ELSE SUMSPECTRA FLAGS AN ERROR!  
from ISISCommandInterface import *
from math import *

wksp="12364_sans_nxs"
CropWorkspace(wksp, OutputWorkspace = "D_T", StartWorkspaceIndex=(8+1-1), EndWorkspaceIndex=str(8+192*192-1))
# sum all into a single spectrum  this does the COLETTE, DISPLAY TIME or D T
SumSpectra("D_T","D_T")
Integration("D_T","D_T_integral",10,100000)	
#
CropWorkspace(wksp, OutputWorkspace = "D_Tf", StartWorkspaceIndex=(8+192*192-1), EndWorkspaceIndex=str(8+192*192*2-1))
# sum all into a single spectrum  this does the COLETTE, DISPLAY TIME or D T
SumSpectra("D_Tf","D_Tf")
Integration("D_Tf","D_Tf_integral",10,100000)	


print beamcoords
print XBEAM_CENTRE
print YBEAM_CENTRE
print RESCALE
print MONITORSPECTRUM
print NUMERATOR
print DENOMINATOR
