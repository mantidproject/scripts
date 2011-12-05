import numpy as np
maskfile='/SNS/SEQ/shared/2011_B/mask_top_bottom.nxs'
V_file='/SNS/SEQ/shared/2011_B/V_files/'+'SEQ_11865_event.nxs'
runs=[range(11900,11904)]
Eguess=70																								#initial energy guess
Erange="-50,0.25,65.0"																						#Energy bins:    Emin,Estep,Emax
datadir="/SNS/SEQ/IPTS-5788/data/"																		#Data directory	
outdir="/SNS/SEQ/IPTS-5788/shared/spes/"	#Output directory
fname_red_sum="horace_info.csv"
fout_prefix="Ei_70_fine_100K"
ang_offset=0.0
angle_name='CCR13VRot'  #Name of the angle to read
maskandnormalize=True	#flag to do the masking and normalization to Vanadium
flag_spe=False                  # flag to generate an spe file
flag_nxspe=True               #flag to generate an nxspe file
single_crystal=False
run_info=[]
def GetEiT0(ws_name,EiGuess):
	"""
	Function to get Ei and  -T0
	
	"""
	alg=GetEi(InputWorkspace=ws_name,Monitor1Spec="1",Monitor2Spec="2",EnergyEstimate=EiGuess)				#Run GetEi algorithm
	[Ei,Tzero]=[float(alg.getPropertyValue("IncidentEnergy")),-float(alg.getPropertyValue("Tzero"))]					#Extract incident energy and T0
	return [Ei,Tzero]
def LoadPathMaker(runs,folder,prefix,suffix):
	"""	
	Function to create paths to files from runnumbers
	return a list of lists with the path, and a corrected list of runs. Files in the inner lists are added together
	side effects: none
	"""
	path=[]
	newruns=[]
	try:
		len(runs)
	except:
		runs=[runs]
	for r in runs:
		try:
			len(r)
		except:
			r=[r]
		temppath=[]
		tempnewruns=[]
		for i in range(len(r)):
			temppath.append(folder+prefix+str(r[i])+suffix)      
			tempnewruns.append(r[i])
			if (not(os.path.isfile(temppath[i]))):
				raise IOError(temppath[i]+" not found")
		path.append(temppath)
		newruns.append(tempnewruns)
	return [path,newruns]
def CreateMasksAndVanadiumNormalization(vanfile,maskfile=''):		
	if (not(os.path.isfile(outdir+"van.nx5"))):
		LoadEventNexus(Filename=vanfile,OutputWorkspace="VAN")
		
		Rebin(InputWorkspace="VAN",OutputWorkspace="VAN",Params="2000,6000,8000",PreserveEvents=False)			#integrate all events between 2000 and 8000 microseconds
		NormaliseByCurrent(InputWorkspace="VAN",OutputWorkspace="VAN")											#normalize by proton charge
		MedianDetectorTest(InputWorkspace="VAN",OutputWorkspace="MASK")			#determine which detectors to mask, and store them in the "MASK" workspace
		if len(maskfile)>0:
		    LoadNexus(Filename=maskfile,OutputWorkspace="temp_mask")
		    MaskDetectors(Workspace="MASK",MaskedWorkspace="temp_mask")		    #add detectors masked in "temp_mask" to "MASK"
		    DeleteWorkspace(Workspace="temp_mask")
		MaskDetectors(Workspace="VAN",MaskedWorkspace="MASK")												#Mask "VAN". This prevents dividing by 0		
		DeleteWorkspace(Workspace="MASK")																	#Mask is carried by VAN workspace
		SaveNexus(InputWorkspace="VAN",Filename=outdir+"van.nx5")
	else:
		LoadNexus(Filename=outdir+"van.nx5",OutputWorkspace="VAN")
		
if (maskandnormalize):
	CreateMasksAndVanadiumNormalization(V_file,maskfile=maskfile)														#Creates two worspaces, one for Vanadium normalization, one for masking

[paths,runs]=LoadPathMaker(runs,datadir,'SEQ_','_event.nxs')
for flist,rlist,i in zip(paths,runs,range(len(paths))):	
	print i+1,"/",len(paths)
	for f,j in zip(flist,range(len(flist))):
		if (j==0):
			LoadEventNexus(Filename=f,OutputWorkspace="IWS")						#Load an event Nexus file
			LoadNexusMonitors(Filename=f,OutputWorkspace="monitor_ws")			#Load monitors    
			#ChangeLogTime(InputWorkspace="IWS",OutputWorkspace="IWS",LogName=angle_name,TimeOffset="859.0")
		else:
			LoadEventNexus(Filename=f,OutputWorkspace="IWS_temp")					#Load an event Nexus file
			LoadNexusMonitors(Filename=f,OutputWorkspace="monitor_ws_temp")			#Load monitors   
		        #ChangeLogTime(InputWorkspace="IWS_temp",OutputWorkspace="IWS_temp",LogName=angle_name,TimeOffset="859.0")	
			Plus(LHSWorkspace="IWS",RHSWorkspace="IWS_temp",OutputWorkspace="IWS")	#Add events to the original workspcace 
			Plus(LHSWorkspace="monitor_ws",RHSWorkspace="monitor_ws_temp",OutputWorkspace="monitor_ws")	#Add  monitors to the original monitor workspcace 
			DeleteWorkspace("IWS_temp")                                                                 			#cleanup
			DeleteWorkspace("monitor_ws_temp")                                                           		#cleanup
	FilterBadPulses(InputWorkspace="IWS",OutputWorkspace = "IWS",LowerCutoff = 50)	
	[Efixed,T0]=GetEiT0("monitor_ws",Eguess)																	#Get Ei and -T0 using the function defined before
	ChangeBinOffset(InputWorkspace="IWS",OutputWorkspace="OWS",Offset=T0)									#Change all TOF by -T0
	NormaliseByCurrent(InputWorkspace="OWS",OutputWorkspace="OWS")												#normalize by proton charge
	ConvertUnits(InputWorkspace="OWS",OutputWorkspace="OWS",Target="Wavelength",EMode="Direct",EFixed=Efixed)	#The algorithm for He3 tube efficiency requires wavelength units
	He3TubeEfficiency(InputWorkspace="OWS",OutputWorkspace="OWS")												#Apply correction due to absorption in He3
	ConvertUnits(InputWorkspace="OWS",OutputWorkspace="OWS",Target="DeltaE",EMode="Direct",EFixed=Efixed)		#Switch  to energy transfer
	CorrectKiKf(InputWorkspace="OWS",OutputWorkspace="OWS")
	Rebin(InputWorkspace="OWS",OutputWorkspace="OWST",Params=Erange,PreserveEvents=False)
	ConvertToDistribution(Workspace="OWST")																		#Convert to differential cross section by dividing by the energy bin width
	if (maskandnormalize):
		MaskDetectors(Workspace="OWST",MaskedWorkspace="VAN")													#apply overall mask
		Divide(LHSWorkspace="OWST",RHSWorkspace="VAN",OutputWorkspace="OWST")									#normalize by Vanadium, if desired
	if single_crystal:
		w=mtd["OWST"]
		psi=np.array(w.getRun()[angle_name].value).mean()+ang_offset
	else:
		psi=0.0
	fname_out="%s%s%g"%(outdir,fout_prefix,psi)
	fname_out=fname_out.replace('.','p')
	rn_infostrng='%s,%g,%g,%g\n'%(fname_out+".spe",Efixed,T0,psi)
	run_info.append(rn_infostrng)
	if flag_spe:
		SaveSPE(InputWorkspace="OWST",Filename=fname_out+".spe")	#save the data. 
	if flag_nxspe:
		SaveNXSPE(InputWorkspace="OWST",Filename=fname_out+".nxspe",Efixed=Efixed,psi=psi,ki_over_kf_scaling=True)
	mtd.sendLogMessage("finished run "+str(rlist))
fhand=open(outdir+fname_red_sum,'w')
fhand.writelines(['runnum,Efixed,T0,psi\n','n,meV,mus,deg.\n'])
fhand.writelines(run_info)
fhand.close()	
