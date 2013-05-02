# SANS2d    1/10/9 change StartSpectrum etc to integer StartWorkspaceIndex not string to follow new convention
# compare SANS2d monitors in wavelength and time (so can phase the choppers)
datapath="\\\isis\inst$\NDXSANS2D\Instrument\data\cycle_09_2\SANS2D"
# only files in runlistnew are read in and summed, those in runlist are put on plots
# this speeds things up if you have already read in a lot of files and don't want to reload them.
# after ~5 raw files you will need to start deleting run_norm workspaces, or uncomment mantid.deleteWorkspace(run+"_norm") in second pass loop below
#runlist=[1660,1472,1355,1674,1026]
runlist=[1660,1472,1355,1674,1026,1017]
runlistnew=[1017]
# uamphr is stored in raw file !  NormaliseByCurrent sorts it.
# INTEGRATION ranges below will need changing to suit detector distance
#
# this really is spectrum number, starting at one
# chop off 4 rows at sides(?) of detector
first=9+2*192
last=192*192+8-2*192
#
numfiles=0
monitor1=0
monitor2=1
#	Load the whole raw file on first pass, and divide by uamphr
for runnum in runlistnew:
	run=str(runnum)
	data_file = datapath + "0000"+run+".raw"
	LoadRaw(Filename = data_file, OutputWorkspace=run+"_norm")
	numfiles=numfiles+1
	NormaliseByCurrent(run+"_norm",run+"_norm")
#
# second pass, cut out the monitors and rear det, then do a D_T sum of rear detector COULD LIMIT RADIUS HERE ?
#for runnum in runlistnew:
#	run=str(runnum)
	CropWorkspace(run+"_norm", run, StartWorkspaceIndex=(first - 1), EndWorkspaceIndex=(last - 1))
	CropWorkspace(run+"_norm", run+"_m1", StartWorkspaceIndex=monitor1, EndWorkspaceIndex=monitor1)
	CropWorkspace(run+"_norm", run+"_m2", StartWorkspaceIndex=monitor2, EndWorkspaceIndex=monitor2)
	#mantid.deleteWorkspace(run+"_norm")
	# sum all into a single spectrum  this does the COLETTE, DISPLAY TIME or D T
	SumSpectra(run,run+"_sum")
	Integration(run+"_sum",run+"_integral",str(15000),str(35000))
	Integration(run+"_sum",run+"_integral_2",str(55000),str(85000))
	#
#
# Now do the plots,  INCLUDE the i=0 if you execute a plot section again !
i=0
for runnum in runlist:
	run=str(runnum)
	if runnum==runlist[0]:
		plot1=plotSpectrum(run+"_sum",0)
		layer=plot1.activeLayer()
		layer.setTitle("D_T")
	else:
		mergePlots(plot1,plotSpectrum(run+"_sum",0))
	layer.setCurveTitle(i,run)
	i=i+1
#
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
