# SANS2d    
#  1/10/9 change StartSpectrum etc to integer StartWorkspaceIndex not string to follow new convention
#  17/03/11 extract integrals automatically, write to workspaces and to a csv file
# 22/03/11 convert to wavelength, also integrals for M2 as well as M3
# 24/03/11 split M1 integral into two wavelength bands as the specturm shape keeps changing
path="\\\isis\inst$\NDXSANS2D\Instrument\data\cycle_10_3\\"
#path="U:\\processed\\"
pfix="sans2d"
sfix=".nxs"
#
#runlist=[5809,5810,5811,5812,5813,5814,5815,5816,5817,5818,5819,5820,5821,5822,5823,5824,5825,5826,5827,5828,5829]
#runlist=[5830,5831,5832,5833,5834,5835,5836,5837,5838,5839,5840]
#runlist=[5841,5842,5843,5844,5845,5846,5847,5848,5849,5850,5851]
#runlist=[5852,5853,5854,5855,5856,5857,5858,5859,5860,5861,5862]
#runlist=[5863,5864,5865,5866,5867,5868,5869,5870,5871,5872,5873]
#runlist=[5874,5875,5876,5877,5878,5879,5880,5881,5882,5883,5884]
#runlist=[5887,5888,5889,5890,5891,5892,5893,5894,5895,5896,5897,5898]
#runlist=[5899,5900,5901,5902,5903,5904,5905,5906,5907,5908,5909]
#runlist=[5910,5911,5912,5913,5914,5915,5916,5917,5918,5919,5920]
#runlist=[5921,5922,5923,5924,5925,5926,5927]
#runlist=[5928,5929,5930,5931,5932,5933,5934,5935,5936,5937,5938]
#runlist=[5939,5940,5941,5942,5943,5944,5945]
#runlist=[5946,5947,5948,5949,5950]
#runlist=[5951,5952,5953,5954,5955,5956,5957,5958]
#runlist=[5959,5960,5961,5962,5963,5964,5965,5966,5967,5968]
#runlist=[5979,5980,5981,5982,5983,5984,5985,5986,5987,5988,5989,5990,5991,5992,5993,5994,5995,5996,5997,5998,5999]
#runlist=[6000,6001,6002,6003,6004,6005,6006,6007,6008,6009,6010,6011,6012,6013,6014,6015,6016,6017,6018,6019,6020]
#runlist=[6095,6103,6111]
runlist=[6705,6706]
#
# uamphr is stored in raw file,  NormaliseByCurrent sorts it.
# INTEGRATION ranges below will need changing to suit detector distance
#
numfiles=0
monitor1=0
monitor2=1
monitor3=2
outlist=[]
sumlist1=[]
errsumlist1=[]
sumlist2=[]
errsumlist2=[]
sumlist3=[]
errsumlist3=[]
sumlist4=[]
errsumlist4=[]
sumlist5=[]
errsumlist5=[]
#
#  OPEN OUTPUT CSV FILE ================================================================================
#
filepath="u:\\user\\"
# mode 'a' appends; 'w' writes; 'r' read etc 
output_file = open(filepath+'monitors2.csv','w')
# write column headers
output_file.write('Run_number,M1,err,M2/M1_short,err,M2/M1_long,err,M3/M1_short,err,M3/M1_long,err\n')
#
#
# LOOP OVER RUNLIST ===================================================================================
#	Load the monitors only and divide by uamphr
for runnum in runlist:
        run=str(runnum)
        nzeros=8-len(run)
        fpad=""
        for ii in range(nzeros):
	    fpad+="0"
#
        filename=path+pfix+fpad+run+sfix
        print "reading file:   "+filename
        check=LoadNexus(Filename=filename,OutputWorkspace=run+"_norm",SpectrumMin=1,SpectrumMax=1).workspace()
	print "uamphr = " + str(check.getSampleDetails().getProtonCharge())
        m1=LoadNexus(Filename=filename,OutputWorkspace=run+"_norm",SpectrumMin=1,SpectrumMax=4)
	numfiles=numfiles+1

	NormaliseByCurrent(run+"_norm",run+"_norm")
#
	CropWorkspace(run+"_norm", run+"_m1", StartWorkspaceIndex=monitor1, EndWorkspaceIndex=monitor1)
	ConvertUnits(run+"_m1",run+"_m1","Wavelength")
	CropWorkspace(run+"_norm", run+"_m2", StartWorkspaceIndex=monitor2, EndWorkspaceIndex=monitor2)
	ConvertUnits(run+"_m2",run+"_m2","Wavelength")
	CropWorkspace(run+"_norm", run+"_m3", StartWorkspaceIndex=monitor3, EndWorkspaceIndex=monitor3)
	ConvertUnits(run+"_m3",run+"_m3","Wavelength")
	InterpolatingRebin(run+"_m1",run+"_m1_reb","2,0.25,12")
	InterpolatingRebin(run+"_m2",run+"_m2_div","2,0.25,12")
	InterpolatingRebin(run+"_m3",run+"_m3_div","2,0.25,12")
	Divide(run+"_m2_div" ,run+"_m1_reb",run+"_m2_div" )
	Divide(run+"_m3_div" ,run+"_m1_reb",run+"_m3_div" )
	#mantid.deleteWorkspace(run+"_norm")
	# sum all into a single spectrum  this does the COLETTE, DISPLAY TIME or D T
#	SumSpectra(run,run+"_sum")
# FIRST INTEGRATION ========================================================================================
	Integration(run+"_m1","integral",str(3.0),str(12.0))
	# THIS EXTRACTS THE SINGLE VALUE OF THE NORMALISED INTEGRAL INTO A PYTHON VARIABLE
	try:
				cr1 = mtd["integral"].readY(0)[0]
				err_cr1 = mtd["integral"].readE(0)[0]
	except AttributeError, exc2:
				if (str(exc2).find('attribute') > 0) and (exc2_value == 1.0):
					print 'WARNING: error reading integral run '+run
					exc2_value = 0.0
					cr1 = -1.0
					err_cr1=0.0
					exc2_value = 1.0
# note use the integer not string version of run number here
	outlist.append(runnum)
	sumlist1.append(cr1)
	errsumlist1.append(err_cr1)
#
# SECOND INTEGRATION ========================================================================================
	Integration(run+"_m2_div","integral",str(3.0),str(5.0))
#
	try:
				cr2 = mtd["integral"].readY(0)[0]
				err_cr2 = mtd["integral"].readE(0)[0]
	except AttributeError, exc2:
				if (str(exc2).find('attribute') > 0) and (exc2_value == 1.0):
					print 'WARNING: error reading integral run '+run
					exc2_value = 0.0
					cr2 = -1.0
					err_cr2=0.0
					exc2_value = 1.0
# note use the integer not string version of run number here
	sumlist2.append(cr2)
	errsumlist2.append(err_cr2)
#	
# THIRD INTEGRATION ========================================================================================
	Integration(run+"_m2_div","integral",str(7.0),str(12.0))
	try:
				cr3 = mtd["integral"].readY(0)[0]
				err_cr3 = mtd["integral"].readE(0)[0]
	except AttributeError, exc2:
				if (str(exc2).find('attribute') > 0) and (exc2_value == 1.0):
					print 'WARNING: error reading integral run '+run
					exc2_value = 0.0
					cr3 = -1.0
					err_cr3=0.0
					exc2_value = 1.0
# note use the integer not string version of run number here
	sumlist3.append(cr3)
	errsumlist3.append(err_cr3)
#
# FOURTH INTEGRATION ========================================================================================
	Integration(run+"_m3_div","integral",str(3.0),str(5.0))
#
	try:
				cr4 = mtd["integral"].readY(0)[0]
				err_cr4 = mtd["integral"].readE(0)[0]
	except AttributeError, exc2:
				if (str(exc2).find('attribute') > 0) and (exc2_value == 1.0):
					print 'WARNING: error reading integral run '+run
					exc2_value = 0.0
					cr4 = -1.0
					err_cr4=0.0
					exc2_value = 1.0
# note use the integer not string version of run number here
	sumlist4.append(cr4)
	errsumlist4.append(err_cr4)
#	
# FITH INTEGRATION ========================================================================================
	Integration(run+"_m3_div","integral",str(7.0),str(12.0))
	try:
				cr5 = mtd["integral"].readY(0)[0]
				err_cr5 = mtd["integral"].readE(0)[0]
	except AttributeError, exc2:
				if (str(exc2).find('attribute') > 0) and (exc2_value == 1.0):
					print 'WARNING: error reading integral run '+run
					exc2_value = 0.0
					cr5 = -1.0
					err_cr5=0.0
					exc2_value = 1.0
# note use the integer not string version of run number here
	sumlist5.append(cr5)
	errsumlist5.append(err_cr5)
#
# WRITE TO A FILE ===========================================================================================
	output_file.write(run + ',' + str(cr1) + ',' + str(err_cr1) + ',' + str(cr2) + \
	',' + str(err_cr2) + ',' + str(cr3) + ',' + str(err_cr3) + ',' + str(cr4) + \
	',' + str(err_cr4) + ',' + str(cr5) + ',' + str(err_cr5)  + '\n') 

output_file.close()
print "written csv file:   "+filepath+'monitors.csv'
#	
#  CREATE RESULTS WORKSPACES ===============================================================================
#
CreateWorkspace("sums_M1",outlist,sumlist1,errsumlist1)
CreateWorkspace("sums_M2div_short",outlist,sumlist2,errsumlist2)
CreateWorkspace("sums_M2div_long",outlist,sumlist3,errsumlist3)
CreateWorkspace("sums_M3div_short",outlist,sumlist4,errsumlist4)
CreateWorkspace("sums_M3div_long",outlist,sumlist5,errsumlist5)
#
# PLOTS,  =============== INCLUDE the i=0 if you execute a plot section again ! ============================================
i=0
for runnum in runlist:
	run=str(runnum)
	if runnum==runlist[0]:
		plot2=plotSpectrum(run+"_m1",0)
		layer=plot2.activeLayer()
		layer.setTitle("monitor 1")
	else:
		mergePlots(plot2,plotSpectrum(run+"_m1",0))
	layer.setCurveTitle(i,run)
	i=i+1
#
i=0
for runnum in runlist:
	run=str(runnum)
	if runnum==runlist[0]:
		plot3=plotSpectrum(run+"_m2",0)
		layer=plot3.activeLayer()
		layer.setTitle("monitor 2")
	else:
		mergePlots(plot3,plotSpectrum(run+"_m2",0))
	layer.setCurveTitle(i,run)
	i=i+1
#
i=0
for runnum in runlist:
	run=str(runnum)
	if runnum==runlist[0]:
		plot4=plotSpectrum(run+"_m3",0)
		layer=plot4.activeLayer()
		layer.setTitle("monitor 3")
	else:
		mergePlots(plot4,plotSpectrum(run+"_m3",0))
	layer.setCurveTitle(i,run)
	i=i+1
