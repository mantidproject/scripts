#
#  25/10/11 takes slices out of SANS2d event mode file and turn to histogram
#
import time as tim
import calendar as cal
#from datetime import timedelta
from ISISCommandInterface import *
# example time handling from Rob Dalgleish (thanks !)   1319280436.54
#t1="2011-Oct-22 10:47:16.539999961"
#t3=t1.split('.')
#print t3
#t4=cal.timegm(tim.strptime(t3[0],"%Y-%b-%d %H:%M:%S"))
#print t4+float('0.'+t3[1])
#
#  ================================================ INPUT DATA HERE =================================     
DataPath("Z:/cycle_11_3/")
samplefilename='SANS2d00010196'
# list of times to slice between,
slicelist=[0,300,900,1200]
#  extensions .n001, .n002, .n003 etc, give first extension number here
savext=1
# ichan is the index for "time zero" in the per pulse proton log.  NOTE presume the proton log is for previous 0.1 sec
# so if you want 0.1 sec resolution you may need to adjust the slicelist numbers very carefully !
ichan=0
# ===================================================================================================
sampleext='.nxs'
# note LoadEventNexus can also do the additional filtering.  It does not by default load the monitors!
# (with LoadMonitors=True a separate workspace appears called SampleEvent_monitors with their 8 histogram spectra)
LoadEventNexus(samplefilename+sampleext,'SampleEvents',LoadMonitors=False)
# We noted that LoadEventNexus does not load all the neccessary matadata, hence use LoadNexus to get the metadata and monitors.
LoadNexus(samplefilename,'SampleMonHist',SpectrumMin=1,SpectrumMax=8)
#CropWorkspace('SampleMonHist','SampleMonHist',StartWorkspaceIndex=0,EndWorkspaceIndex=7)
#
ws=mtd["SampleEvents"]
r= ws.getRun()
time_array = r.getLogData("proton_charge").times
value_array = r.getLogData("proton_charge").value
#
unit_array = r.getLogData("proton_charge").units
#print dir(r.getLogData("proton_charge"))
#print dir(time_array)
#
tarr=[]
# converting individual elements of the time_arry with str(time_array[i]) results in format change and loss of the fractions of sec
# however Rob discovered (by accident) that str( ) of the whole array at once does work, then need to split it up again ...
tarr2str=str(time_array)
tarr2str=tarr2str[1:len(tarr2str)-1].split(',')
print " some debugging info ..."
print time_array[0],tarr2str[0]
print time_array[1],tarr2str[1]
print time_array[2],tarr2str[2]
print time_array[3],tarr2str[3]
print time_array[4],tarr2str[4]
#
# choose start time index
ichan=0
t1=tarr2str[ichan]
t2=t1.split('.')
t3=cal.timegm(tim.strptime(t2[0],"%Y-%b-%d %H:%M:%S"))
print t1,t2,t3
tstart=float(t3)+float('0.'+t2[1])
npts=len(tarr2str)
print "start index=",ichan,"  ",t1,"   mintime=",tstart,"  npts=",npts
#
print " more debugging info ..."
for i in range(npts):
	t1=tarr2str[i]
	t2=t1.split('.')
	# aagh is falling over when there are no fractions of a second, so no '.' and no t2[1]
	t2.append("0")
	t3=cal.timegm(tim.strptime(t2[0],"%Y-%b-%d %H:%M:%S"))
	
	tarr.append(float(t3)+float('0.'+t2[1])-tstart)
	if i<20:
		print t1,t3,tarr[i],value_array[i]
tmax=tarr[npts-1]
print "first 20 times =",tarr[0:20]
print "last 20 times =",tarr[npts-21:npts-1]
#
def findCharge(tarr,varr,tmin,tmax):
	imin=len(tarr)
	imax=0
	for i in range(len(tarr)):
		imin=i
		if tarr[i]>=tmin:
			break
	for i in range(len(tarr)):
		imax=i
		if tarr[i]>=tmax:
			break
	tCharge=0.0
	for i in range(imin,imax):
		tCharge=tCharge+varr[i]
	return tCharge

totaluamphr=findCharge(tarr,value_array,0.0,tmax)
print "total proton charge=",totaluamphr
#
nslices=len(slicelist)-1
for i in range(nslices):
	nzeros=3-len(str(savext))
	fpad=".n"
	for ii in range(nzeros):
		fpad+="0"
	newfilename=samplefilename+fpad+str(savext)
	savext=savext+1
  	# cut out a range of data collection times
	FilterByTime('SampleEvents','SampleEventsCut',StartTime=slicelist[i],Stoptime=slicelist[i+1])
	#FilterByTime('SampleEvents','SampleEventsCut',AbsoluteStartTime='2011-10-09T22:27:00',AbsoluteStoptime='2011-10-09T22:43:00')
	# impose time channels, identical to ones used in histogram part of file (blame Kelvin for the odd choices here)
	#params='5.5,4000,100005'
	params='5.5,44.5,50,50,1000,500,1500,750,99750,255,100005'
	Rebin('SampleEventsCut','SampleHist',params,PreserveEvents=False)
	#Rebin('SampleEvents','SampleHist',params,PreserveEvents=False)
	# use call to crop to make a copy of the monitor workspace to put at front of new file
	# will later rescale the new monitors to the correct run time fraction
	# was CropWorkspace('SampleEvents_monitors','SampleHistNew',StartWorkspaceIndex=0,EndWorkspaceIndex=7)
	uamphr=findCharge(tarr,value_array,slicelist[i],slicelist[i+1])
	print "slice ",i+1,"  from",slicelist[i],"  to ",slicelist[i+1],"  uamphr=",uamphr,"  of total=",totaluamphr
	rescale=uamphr/totaluamphr
	Scale('SampleMonHist','SampleHistNew',rescale)	
	# 8+ ((192^2)*2) -1 = 73735 for two detectors
	# 8 + 192^2 -1 =36871  use only the spectra we need, 
	# BEWARE due to Kelvin's wring tables Spectrum 9 is in ID=4
	CropWorkspace('SampleHist','SampleHistCropped',StartWorkspaceIndex=4,EndWorkspaceIndex=73731)
	#CropWorkspace('SampleHist','SampleHistCropped',StartWorkspaceIndex=4,EndWorkspaceIndex=36867)
	# conjoin expands the first workspace by tagging the second to the end of it, time channels I think have to
	# be identical and spectrum number ranges must not overlap.
	ConjoinWorkspaces('SampleHistNew','SampleHistCropped')
	#
	# SampleHist seems to have main detector counts OK, but the monitor spectra which were in hist mode to start with are ZERO
	SaveNexusProcessed('SampleHistNew',newfilename,PreserveEvents=False)
	print "saved file: ",newfilename
print 'DELETE all the temp workspaces else you may run out of memory!!'
print 'BEWARE Mantid reduction does not delete _n001 etc workspaces, and you need to rename the output ones before reducing another file!'
