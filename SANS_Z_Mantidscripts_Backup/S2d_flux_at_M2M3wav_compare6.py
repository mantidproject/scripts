# SANS2d    
#  1/10/9 change StartSpectrum etc to integer StartWorkspaceIndex not string to follow new convention
#  17/03/11 extract integrals automatically, write to workspaces and to a csv file
# 22/03/11 convert to wavelength, also integrals for M2 as well as M3
path="\\\isis\inst$\NDXSANS2D\Instrument\data\cycle_10_1\\"
#path="U:\\processed\\"
pfix="sans2d"
sfix=".nxs"
#
runlist=[5544,5545,5546,5332,5333,5300,5302,4734,4733,4635,4634,4373,4374]
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
output_file = open(filepath+'monitors.csv','w')
# write column headers
output_file.write('Run_number,M1,err,M2_short,err,M2_long,err,M3_short,err,M3_long,err\n')
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
	Integration(run+"_m2","integral",str(3.0),str(5.0))
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
	Integration(run+"_m2","integral",str(7.0),str(12.0))
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
	Integration(run+"_m3","integral",str(3.0),str(5.0))
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
	Integration(run+"_m3","integral",str(7.0),str(12.0))
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
CreateWorkspace("sums_M2_short",outlist,sumlist2,errsumlist2)
CreateWorkspace("sums_M2_long",outlist,sumlist3,errsumlist3)
CreateWorkspace("sums_M3_short",outlist,sumlist4,errsumlist4)
CreateWorkspace("sums_M3_long",outlist,sumlist5,errsumlist5)
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
