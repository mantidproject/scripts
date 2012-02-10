from qtiGenie import *
iliad_setup('let')

# program to crunch down event mode from LET and produce output SPE files. Program can automatically find all incident energies in rep rate mode and write out spe files in the 
# form of LET'run no: +ei'.spe

#############################################
# this is the user input section
wb=4211    # enter whitebeam run number here
run_no=[4962]   # enter event mode run numbers here or use next line for a continous sequence of runs i.e range(first run, last run +1)
#run_no=range(4181,4187)
ei=[5]           # incident energies you want analysed, or leave as ei=[]  if you want all incident energies analysed
ebin=[-0.2,0.002,0.8]    #binning of the energy for the spe file. The numbers are as a fraction of ei [from ,step, to ]
#mapping='rings_103'  # rings mapping file for powders, liquout=iliad("wb_wksp","w1reb",energy,ebinstring,mapping,bleed=False,norm_method='current',det_cal_file='det_corrected7.dat',detector_van_range=[0.5,200],bkgd_range=[int(t_elastic),int(tmax)])
mapping='one2one_103' # one to one mapping of detectors for crystals
file = '/Users/jon/Work-computing/mantid_test_data/let/hard.msk'    # standard hard mask file  for LET
#file = '/media/sf_D_DRIVE/MantidInstall/LET_maps/magnet_hard.msk'    # standard hard mask file  for LET
############################################




##########################

LoadRaw(Filename=str(wb),OutputWorkspace="wb_wksp") # load whitebeam
#LoadEventNexus(Filename='D:\temp\LET00004286.nxs',OutputWorkspace='wb_wksp',SingleBankPixelsOnly='0',LoadMonitors='1',MonitorsAsEvents='1')
#ConjoinWorkspaces(InputWorkspace1='wb_wksp',InputWorkspace2='wb_wksp_monitors')
#Rebin(InputWorkspace='wb_wksp',OutputWorkspace='wb_wksp',Params=[15000,50,40000],PreserveEvents='0')
#Rebin(InputWorkspace='wb_wksp_monitors',OutputWorkspace='wb_wksp_monitors',Params=[15000,50,40000])
#MaskDetectors(Workspace='wb_wksp',SpectraList=sl)
	######################################################################




for run in run_no:     #loop around runs
	LoadEventNexus(Filename='LET0000'+str(run)+'.nxs',OutputWorkspace='w1',SingleBankPixelsOnly='0',LoadMonitors='1',MonitorsAsEvents='1')
	#MaskDetectors(Workspace='w1',SpectraList=sl)	
	Rebin(InputWorkspace='w1_monitors',OutputWorkspace='mon',Params=[1000,100,90000])
	ExtractSingleSpectrum(InputWorkspace='mon',OutputWorkspace='mon',WorkspaceIndex='5')  #extract monitor 6
	ConvertToMatrixWorkspace(InputWorkspace='mon',OutputWorkspace='mon')
	ConvertUnits(InputWorkspace='mon',OutputWorkspace='mon',Target='Energy')
	NormaliseByCurrent(InputWorkspace='mon',OutputWorkspace='mon')     #monitor 6 converted to energy and normalised
	ConjoinWorkspaces(InputWorkspace1='w1',InputWorkspace2='w1_monitors')
	##################################
	# this section finds all the transmitted incident energies
	if len(ei) == 0:
		for x in range(0,15):
			Max(InputWorkspace='mon',OutputWorkspace='maxval')
			mv=mtd['maxval']
			if mv.dataY(0)[0] >= 250:
				min=mv.dataX(0)[0] -0.02
				max=mv.dataX(0)[1] +0.02
				RemoveBins(InputWorkspace='mon',OutputWorkspace='mon',XMin=min,XMax=max)
				ei.append(mv.dataX(0)[0])
	ei.sort()     #sorts energies into order
	ei = [ '%.2f' % elem for elem in ei ]  
	print 'energies transmitted are:'
	print (ei)

	for energy in ei:
		energy=float(energy)
		print (energy)
		emin=0.2*energy   #minimum energy is with 80% energy loss
		lam=(81.81/energy)**0.5
		lam_max=(81.81/emin)**0.5
		tsam=252.82*lam*25   #time at sample
		tmon2=252.82*lam*23.5 #time to monitor 6 on LET
		tmax=tsam+(252.82*lam_max*4.1) #maximum time to measure inelastic signal to
		t_elastic=tsam+(252.82*lam*4.1)   #maximum time of elastic signal
		tbin=[int(tmon2),1.6,int(tmax)]
		Rebin(InputWorkspace='w1',OutputWorkspace='w1reb',Params=tbin,PreserveEvents='1')	

		
		energybin=[ebin[0]*energy,ebin[1]*energy,ebin[2]*energy]
		energybin = [ '%.4f' % elem for elem in energybin ]  
		ebinstring=str(energybin[0])+','+str(energybin[1])+','+str(energybin[2])
		print ebinstring
		out=iliad("wb_wksp","w1reb",energy,ebinstring,mapping,bleed=False,norm_method='current',det_cal_file='det_corrected7.dat',detector_van_range=[0.1,500],bkgd_range=[34000,35000])
		#out=iliad("wb_wksp","w1reb",energy,ebinstring,mapping,bleed=False,norm_method='current',det_cal_file='det_corrected7.dat',detector_van_range=[0.5,200],diag_remove_zero=True)
		SaveNXSPE(out,'LET'+str(run)+'_'+str(energy)+'mev.nxspe')
		#SaveSPE(out,'/home/let/let_data_share/LET'+str(run)+'_'+str(energy)+'mev.spe')
