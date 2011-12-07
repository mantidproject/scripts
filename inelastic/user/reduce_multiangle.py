import numpy as np
#maskfile='/SNS/SEQ/shared/2011_B/mask_top_bottom
maskfile=''
V_file='/SNS/ARCS/IPTS-3824/0/12730/NeXus/ARCS_12730_event.nxs'


runs=[12924,12925,12926,12927,12928,12929,12930,12932]
Eguess=130										#initial energy guess
Erange="-50,1,125.0"							#Energy bins:    Emin,Estep,Emax
datadir="/SNS/ARCS/IPTS-4607/0/"					#Data directory	
outdir="/SNS/ARCS/IPTS-4607/shared/mantid/"		#Output directory
ang_offset=88.0
angle_name='CCR12Rot'  							#Name of the angle to read
angle_step=1
maskandnormalize=True							#flag to do the masking and normalization to Vanadium
flag_spe=False                  							# flag to generate an spe file
flag_nxspe=True               							#flag to generate an nxspe file

def GetEiT0(ws_name,EiGuess):
	"""
	Function to get Ei and  -T0
	
	"""
	alg=GetEi(InputWorkspace=ws_name,Monitor1Spec="1",Monitor2Spec="2",EnergyEstimate=EiGuess)		#Run GetEi algorithm
	[Ei,Tzero]=[float(alg.getPropertyValue("IncidentEnergy")),-float(alg.getPropertyValue("Tzero"))]		#Extract incident energy and T0
	return [Ei,Tzero]

def CreateMasksAndVanadiumNormalization(vanfile,maskfile=''):		
	if (not(os.path.isfile(outdir+"van.nx5"))):
		LoadEventNexus(Filename=vanfile,OutputWorkspace="VAN")
		Rebin(InputWorkspace="VAN",OutputWorkspace="VAN",Params="3000,2000,5000",PreserveEvents=False)		#integrate all events between 3000 and 5000 microseconds
		NormaliseByCurrent(InputWorkspace="VAN",OutputWorkspace="VAN")										#normalize by proton charge
		MedianDetectorTest(InputWorkspace="VAN",OutputWorkspace="MASK")										#determine which detectors to mask, and store them in the "MASK" workspace
		if len(maskfile)>0:
		    LoadNexus(Filename=maskfile,OutputWorkspace="temp_mask")
		    MaskDetectors(Workspace="MASK",MaskedWorkspace="temp_mask")		    								#add detectors masked in "temp_mask" to "MASK"
		    DeleteWorkspace(Workspace="temp_mask")
		MaskDetectors(Workspace="VAN",MaskedWorkspace="MASK")												#Mask "VAN". This prevents dividing by 0		
		DeleteWorkspace(Workspace="MASK")																	#Mask is carried by VAN workspace
		SaveNexus(InputWorkspace="VAN",Filename=outdir+"van.nx5")
	else:
		LoadNexus(Filename=outdir+"van.nx5",OutputWorkspace="VAN")											#if Vanadium is processed, just load it
		
if (maskandnormalize):
	CreateMasksAndVanadiumNormalization(V_file,maskfile=maskfile)												#Creates two worspaces, one for Vanadium normalization, one for masking


for r,i in zip(runs,range(len(runs))):
	f=datadir+str(r)+'/NeXus/ARCS_'+str(r)+'_event.nxs'
	if (i==0):
		LoadEventNexus(Filename=f,OutputWorkspace="IWS")					#Load an event Nexus file
		LoadNexusMonitors(Filename=f,OutputWorkspace="monitor_ws")			#Load monitors    
	else:
		LoadEventNexus(Filename=f,OutputWorkspace="IWS_temp")					#Load an event Nexus file
		LoadNexusMonitors(Filename=f,OutputWorkspace="monitor_ws_temp")			#Load monitors   
		#ChangeLogTime(InputWorkspace="IWS_temp",OutputWorkspace="IWS_temp",LogName=angle_name,TimeOffset="859.0")	
		Plus(LHSWorkspace="IWS",RHSWorkspace="IWS_temp",OutputWorkspace="IWS")	#Add events to the original workspcace 
		Plus(LHSWorkspace="monitor_ws",RHSWorkspace="monitor_ws_temp",OutputWorkspace="monitor_ws")	#Add  monitors to the original monitor workspcace 
		DeleteWorkspace("IWS_temp")                                                                 			#cleanup
		DeleteWorkspace("monitor_ws_temp")                                                           		#cleanup
FilterBadPulses(InputWorkspace="IWS",OutputWorkspace = "IWS",LowerCutoff = 50)	
[Efixed,T0]=GetEiT0("monitor_ws",Eguess)

ChangeBinOffset(InputWorkspace="IWS",OutputWorkspace="IWS",Offset=T0)				#Change all TOF by -T0
w=mtd["IWS"]
psi=np.array(w.getRun()[angle_name].value)  #Add angle offset when you write the file.
psimin=min(psi)
psimax=max(psi)
psiarray=np.arange(psimin+0.5*angle_step,psimax,angle_step)

for ang in psiarray:
	FilterByLogValue(InputWorkspace="IWS",OutputWorkspace="OWS",LogName=angle_name,MinimumValue=ang-0.5*angle_step,MaximumValue=ang+0.5*angle_step)
	NormaliseByCurrent(InputWorkspace="OWS",OutputWorkspace="OWS")												#normalize by proton charge
	ConvertUnits(InputWorkspace="OWS",OutputWorkspace="OWS",Target="Wavelength",EMode="Direct",EFixed=Efixed)		#The algorithm for He3 tube efficiency requires wavelength units
	He3TubeEfficiency(InputWorkspace="OWS",OutputWorkspace="OWS")												#Apply correction due to absorption in He3
	ConvertUnits(InputWorkspace="OWS",OutputWorkspace="OWS",Target="DeltaE",EMode="Direct",EFixed=Efixed)			#Switch  to energy transfer
	CorrectKiKf(InputWorkspace="OWS",OutputWorkspace="OWS")														#ki/kf
	Rebin(InputWorkspace="OWS",OutputWorkspace="OWSH",Params=Erange,PreserveEvents=False)							# forget events at this point
	ConvertToDistribution(Workspace="OWSH")																		#Convert to differential cross section by dividing by the energy bin width
	if (maskandnormalize):
		MaskDetectors(Workspace="OWSH",MaskedWorkspace="VAN")													#apply overall mask
		Divide(LHSWorkspace="OWSH",RHSWorkspace="VAN",OutputWorkspace="OWSH")									#divide by vanadium
	fname_out=outdir+"Angle_"+str(ang+ang_offset)+"_deg"
	if flag_spe:
		SaveSPE(InputWorkspace="OWSH",Filename=fname_out+".spe")	#save the data. 
	if flag_nxspe:
		SaveNXSPE(InputWorkspace="OWSH",Filename=fname_out+".nxspe",Efixed=Efixed,psi=ang + ang_offset,KiOverKfScaling=True)
	mtd.sendLogMessage("finished angle "+str(ang))



