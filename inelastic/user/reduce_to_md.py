import os

#maskfile='/SNS/SEQ/shared/2011_B/mask_40.nxs'
#V_file='/SNS/SEQ/IPTS-5170/data/SEQ_11089_event.nxs'
V_file='/SNS/SEQ/shared/2011_B/V_files/'+'SEQ_10901_event.nxs'
#runs =  range(11016,11089)+range(11090,11106)
#list if runs that have the same energy range(first_run,last_run+1)
runs =  range(11499,11561)
#runs =  range(11499,11511)
#initial energy guess
Eguess=50
#Energy bins:    Emin,Estep,Emax
Erange="-10.0,0.2,45.0"
#Data directory	
datadir="/SNS/SEQ/IPTS-4783/data/"
#Output directory
outdir="/SNS/users/2zr/data/Mantid/"
#flag to do the masking and normalization to Vanadium
maskandnormalize=True
#GetEi and T0 function: note that in different versions of Mantid
#the Monitor1Spec and monitor2Spec might be "0" and "1"
def GetEiT0(ws_name,EiGuess):
    # Run GetEi algorithm
    [Ei, MonPeak,
     MonIndex, Tzero] = GetEi(InputWorkspace=ws_name, Monitor1Spec="1",
                              Monitor2Spec="2", EnergyEstimate=EiGuess)
    # Extract incident energy and T0
    return [Ei,Tzero]	

def CreateMasksAndVanadiumNormalization(vanfile,maskfile=''):		
    if (not(os.path.isfile(outdir+"van.nx5"))):
        LoadEventNexus(Filename=vanfile, OutputWorkspace="VAN")
        # integrate all events between 2000 and 8000 microseconds
        Rebin(InputWorkspace="VAN", OutputWorkspace="VAN",
              Params="2000,6000,8000", PreserveEvents=False)
        # normalize by proton charge
        NormaliseByCurrent(InputWorkspace="VAN",OutputWorkspace="VAN")
        # determine which detectors to mask, and store them in
        # the "MASK" workspace       
        MedianDetectorTest(InputWorkspace="VAN",
                           OutputWorkspace="MASK")			
        if len(maskfile)>0:
            LoadNexus(Filename=maskfile,
                      OutputWorkspace="temp_mask")
            # add detectors masked in "temp_mask" to "MASK"
            MaskDetectors(Workspace="MASK",
                          MaskedWorkspace="temp_mask")
            DeleteWorkspace(Workspace="temp_mask")
            # Mask "VAN". This prevents dividing by 0		
            MaskDetectors(Workspace="VAN", MaskedWorkspace="MASK")
            # Mask is carried by VAN workspace     
            DeleteWorkspace(Workspace="MASK")
            SaveNexus(InputWorkspace="VAN",
                      Filename=outdir+"van.nx5")
    else:
        LoadNexus(Filename=outdir+"van.nx5", OutputWorkspace="VAN")
        
if (maskandnormalize):
    # Creates two workspaces, one for Vanadium normalization, one for
    # masking
    CreateMasksAndVanadiumNormalization(V_file)
    
run_info=[]	
for i,run in enumerate(runs):
    filename=datadir+"SEQ_"+str(run)+"_event.nxs"
    # Load an event Nexus file
    LoadEventNexus(Filename=filename,OutputWorkspace="IWS")
    # FilterBadPulses(InputWorkspace="IWS",
    # OutputWorkspace = "IWS",LowerCutoff = 50)
    LoadNexusMonitors(Filename=filename, OutputWorkspace="MonWS")
    # Get Ei and -T0 using the function defined before
    [Efixed,T0]=GetEiT0("MonWS", Eguess)
    # Change all TOF by -T0
    ChangeBinOffset(InputWorkspace="IWS", OutputWorkspace="OWS", Offset=T0)
    # normalize by proton charge
    NormaliseByCurrent(InputWorkspace="OWS",OutputWorkspace="OWS")
    # The algorithm for He3 tube efficiency requires wavelength units
    ConvertUnits(InputWorkspace="OWS", OutputWorkspace="OWS", 
                 Target="Wavelength", EMode="Direct", EFixed=Efixed)	
    # Apply correction due to absorption in He3
    He3TubeEfficiency(InputWorkspace="OWS", OutputWorkspace="OWS")
    # Switch  to energy transfer
    ConvertUnits(InputWorkspace="OWS", OutputWorkspace="OWS",
                 Target="DeltaE", EMode="Direct", EFixed=Efixed)		
    # Apply k_i/k_f factor
    CorrectKiKf(InputWorkspace="OWS", OutputWorkspace="OWS")
    # Make sure the bins are correct
    #Rebin(InputWorkspace="OWS", OutputWorkspace="OWS",
    #      Params=Erange, PreserveEvents=False)
    # Convert to differential cross section by dividing by the
    # energy bin width
    #ConvertToDistribution(Workspace="OWS")
    if (maskandnormalize):
        # apply overall mask
        MaskDetectors(Workspace="OWS", MaskedWorkspace="VAN")
        # normalize by Vanadium, if desired
        Divide(LHSWorkspace="OWS", RHSWorkspace="VAN",
               OutputWorkspace="OWS")

    # Rename workspace to something meaningful
    ws_name = "SEQ_"+str(run)
    RenameWorkspace(InputWorkspace="OWS",  OutputWorkspace=ws_name)
    # Need to fix the goniometer angle by 49.73 degrees
    w = mtd[ws_name]
    psi = w.getRun()["CCR13VRot"].getStatistics().mean + 49.73
    AddSampleLog(Workspace=ws_name, LogName="CCR13VRot_Fixed",
                 LogType="Number Series", LogText=str(psi))
    # Set the Goiniometer information
    SetGoniometer(Workspace=ws_name, Axis0="CCR13VRot_Fixed,0,1,0,1")
    # Set the information for the UB matrix
    SetUB(Workspace=ws_name,
          a=3.643, b=3.643, c=5.781, alpha=90, beta=90, gamma=120,
          u='1,1,0', v='0,0,1')
    # Create the MDEventWorkspace
    md_ws_name=ws_name+'_md'
    ConvertToMDEvents(InputWorkspace=ws_name, OutputWorkspace=md_ws_name,
                      QDimensions='QhQkQl', MinValues='-5,-5,-5,-10',
                      MaxValues='5,5,5,45', MaxRecursionDepth='1')
    #ConvertToMDEvents(InputWorkspace=ws_name, OutputWorkspace=md_ws_name,
    #                  QDimensions='QhQkQl', MinValues='0,-0.75,-0.9,17',
    #                  MaxValues='0.5,-0.25,0.9,22', MaxRecursionDepth='1')
    # Remove SPE workspace
    DeleteWorkspace(Workspace=ws_name)
    # Save MDEventWorkspace to a file
    md_path = os.path.join(outdir, md_ws_name+'.nxs')
    SaveMD(InputWorkspace=md_ws_name, Filename=md_path)
    # Save filename for later merge
    run_info.append(md_path)
    # Delete MDEventWorkspace for safety
    DeleteWorkspace(Workspace=md_ws_name)
    logger.notice("Finished run " + str(run))

# Create merged MD dataset
final_md_name='SEQ_Gd'
#final_md_name='SEQ_Gd_small'
MergeMDFiles(Filenames=",".join(run_info),
             OutputFilename=os.path.join(outdir, final_md_name+'.nxs'),
             OutputWorkspace=final_md_name)
