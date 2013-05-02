# 22/1/13 look for average and peak count rates on LOQ
# version 1 crashed Mantid after 1499 entries to the .csv file - WHY ???
# version 2, don't store the D_T array, commented out all diagnostic print statements
from mantidsimple import *
path="I:/data/cycle_10_2/"
#path="U:\\processed\\"
pfix="LOQ"
sfix=".nxs"
# Beware this script may fall over with periods files or other unusual set ups, the loop below only check number of spectra and microamphrs
runlist1=58661
runlist2=61160
outname="survey_102b"
# LOQ spectrum range
spectrum_start=3
spectrum_end=3+128*128-1
pfix="LOQ"
sfix=".nxs"
wksp="LOQdata"
# number of digits in run numbers
len_runnum=5
#
# uamphr is stored in raw file,  NormaliseByCurrent sorts it.
# INTEGRATION ranges below will need changing to suit detector distance
#
sumlist1=[]
sumlist2=[]
sumlist1=[]
errsumlist1=[]
errsumlist2=[]
#
#  OPEN OUTPUT CSV FILE ================================================================================
#
filepath="z:\\user\\"
# mode 'a' appends; 'w' writes; 'r' read etc 
output_file = open(filepath+outname+'.csv','w')
# write column headers
output_file.write('Run_number,avg kHz,instant kHz\n')
#
#
# LOOP OVER RUNLIST ===================================================================================
#	Load the monitors only and divide by uamphr
for runnum in range(runlist1,runlist2+1):
        run=str(runnum)
	# get the first file in list,   5  digits for LOQ, 8 for SANS2d
	nzeros=len_runnum-len(run)
	fpad=""
	for ii in range(nzeros):
		fpad+="0"
	#
	filename=path+pfix+fpad+run+sfix
	#print "reading file:   "+filename
	#print runnum
        check=LoadISISNexus(Filename=filename,OutputWorkspace=wksp,SpectrumMin=1,SpectrumMax=1)
        nspectra = mtd[wksp].getRun().getLogData("nSpectra").value
        good_frames = mtd[wksp].getRun().getLogData("goodfrm").value 
	ntc = mtd[wksp].getRun().getLogData("nchannels").value
	#print "nspectra, good_frames, ntc =",nspectra,good_frames,ntc
	if(nspectra>=spectrum_end)and (good_frames>0):
		#print "uampHr = ", mtd[wksp].getRun().getProtonCharge()
		check=LoadISISNexus(Filename=filename,OutputWorkspace=wksp,SpectrumMin=spectrum_start,SpectrumMax=spectrum_end)
		SumSpectra(wksp,wksp)
		# wksp here contains actual counts so can sum them easily
		Integration(wksp,"integral")
		average=mtd["integral"].readY(0)[0]
		# need convert to cts/microsec here, else may not get the proper peak rate
		ConvertToDistribution(wksp)
		#yy=[]
		#xx=[]
		instant=0
		yy=mtd[wksp].readY(0)[0:ntc] 
		xx=mtd[wksp].readX(0)[0:ntc] 
		# then get an array, not a python list, so index(yy) does not work, so do from scratch here
		# find the peak rate in counts per micrisecond, as XX is in microsec
		for i in range(ntc):
			#yy1=mtd[wksp].readY(0)[i]
			#yy.append(yy1)
			#xx1=mtd[wksp].readX(0)[i]
			#xx.append(xx1)
			if(yy[i]>instant):
				instant=yy[i]
				ii=i
		#print ii, xx[ii],yy[ii],average,instant
		# average & max count rates in kHz
		average=1000.0*average/(good_frames*(xx[ntc-1]-xx[0]))
		instant=1000.0*instant/good_frames
		print runnum," avg = ",average,"  inst = ",instant,"  kHz"
	# note use the integer not string version of run number here
		#outlist.append(runnum)
		#sumlist1.append(average)
		#sumlist2.append(instant)
		#errsumlist1.append(0.0)
		#errsumlist2.append(0.0)
# WRITE TO A FILE ===========================================================================================
		output_file.write(run + ',' + str(average) + ',' + str(instant) + '\n') 

output_file.close()
print "written csv file:   "+filepath+outname+'.csv'
#	
#  CREATE RESULTS WORKSPACES ===============================================================================
#
#CreateWorkspace("average",outlist,sumlist1,errsumlist1)
#CreateWorkspace("instant",outlist,sumlist2,errsumlist2)
#
