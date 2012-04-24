import os
import sys
import shutil
sys.path.append("/opt/Mantid/bin")
#from MantidFramework import mtd
#mtd.initialise()
from mantid.simpleapi import *
from numpy import *
from string import *
from MaskBTP import *

class V_norm_obj(object):
   def __init__(self,vanfile,wkspcnm,wlstr,outdir,maskfile='',maskw=None,ld_saved_fl=True):
	 self.outdir=outdir
	 self.vanfile=vanfile
	 self.maskfile=maskfile
	 self.maskw=maskw
	 self.wlstr=wlstr
	 self.ld_saved_fl=ld_saved_fl
	 self.wkspcnm=wkspcnm
   def CreateMasksAndVanadiumNormalization(self):
        """
	"""
	if (os.path.isfile(self.outdir+"van.nx5") and (self.ld_saved_fl)):
	   LoadNexus(Filename=outdir+"van.nx5",OutputWorkspace="VAN")
	else:   
	   LoadEventNexus(Filename=self.vanfile,OutputWorkspace="VAN")
	   ConvertUnits(InputWorkspace="VAN",OutputWorkspace="VAN",Target="Wavelength",EMode="Elastic")
	   Rebin(InputWorkspace="VAN",OutputWorkspace="VAN",Params=self.wlstr,PreserveEvents=False)			#integrate all events between 2000 and 8000 microseconds
	   ConvertToDistribution("VAN")
	   NormaliseByCurrent(InputWorkspace="VAN",OutputWorkspace="VAN")											#normalize by proton charge
	   MedianDetectorTest(InputWorkspace="VAN",OutputWorkspace="MASK")			#determine which detectors to mask, and store them in the "MASK" workspace
	   if self.maskw!=None:
		   MaskDetectors(Workspace="MASK",MaskedWorkspace=self.maskw)
	   if len(self.maskfile)>0:
	      LoadNexus(Filename=self.maskfile,OutputWorkspace="temp_mask")
	      MaskDetectors(Workspace="MASK",MaskedWorkspace="temp_mask")		    #add detectors masked in "temp_mask" to "MASK"
	      DeleteWorkspace(Workspace="temp_mask")
	   MaskDetectors(Workspace="VAN",MaskedWorkspace="MASK")												#Mask "VAN". This prevents dividing by 0		
	   DeleteWorkspace(Workspace="MASK")																	#Mask is carried by VAN workspace
	   SaveNexus(InputWorkspace="VAN",Filename=self.outdir+"van.nx5")
	CloneWorkspace(InputWorkspace="VAN",OutputWorkspace=self.wkspcnm)
	DeleteWorkspace(Workspace="VAN")



def GetEiT0(ws_name,EiGuess):
	
	alg=GetEi(InputWorkspace=ws_name,Monitor1Spec="1",Monitor2Spec="2",EnergyEstimate=float(EiGuess))				#Run GetEi algorithm
	Ei=alg[0]
	Tzero=-alg[3]					#Extract incident energy and T0
	return [Ei,Tzero]

def load_n_sum(runs,datadir,wkspcnm):
   """
   given a list of runs in datadir this will sum the runs and reduce them into a
   workspace wkspcnm
   """
   for i,run in enumerate(runs):
	filename=datadir+"SEQ_"+str(run)+"_event.nxs"
	LoadEventNexus(Filename=filename,OutputWorkspace="IWS")												#Load an event Nexus file
	FilterBadPulses(InputWorkspace="IWS",OutputWorkspace = "IWS",LowerCutoff = 50)
	LoadNexusMonitors(Filename=filename,OutputWorkspace="MonWS")
	if i ==0:
		CloneWorkspace(InputWorkspace="IWS",OutputWorkspace="IWS_sum")
		CloneWorkspace(InputWorkspace="MonWS",OutputWorkspace = "MonWS_sum")
	else:
		Plus(LHSWorkspace = "IWS_sum",RHSWorkspace="IWS",OutputWorkspace = "IWS_sum")
		Plus(LHSWorkspace = "MonWS_sum",RHSWorkspace ="MonWS",OutputWorkspace = "MonWS_sum")
   Mons=mtd['MonWS_sum']		
   Eguess=array(Mons.getRun()['EnergyRequest'].value).mean()
   logger.notice("Ei (Guess) = " + str(Eguess))
   [Efixed,T0]=GetEiT0("MonWS_sum",Eguess)
   Estep=Efixed*0.02
   Elow=-Efixed*0.5
   Ehigh=Efixed*0.95
   Erange='%g,%g,%g'%(Elow,Estep,Ehigh)													#Get Ei and -T0 using the function defined before
   ChangeBinOffset(InputWorkspace="IWS_sum",OutputWorkspace="OWS",Offset=T0)    #Change all TOF by -T0
   DeleteWorkspace(Workspace="IWS_sum")
   DeleteWorkspace(Workspace="IWS")
   DeleteWorkspace(Workspace="MonWS")
   NormaliseByCurrent(InputWorkspace="OWS",OutputWorkspace="OWS")												#normalize by proton charge
   ConvertUnits(InputWorkspace="OWS",OutputWorkspace="OWS",Target="Wavelength",EMode="Direct",EFixed=Efixed)	#The algorithm for He3 tube efficiency requires wavelength units
   He3TubeEfficiency(InputWorkspace="OWS",OutputWorkspace="OWS")												#Apply correction due to absorption in He3
   ConvertUnits(InputWorkspace="OWS",OutputWorkspace="OWS",Target="DeltaE",EMode="Direct",EFixed=Efixed)		#Switch  to energy transfer
   CorrectKiKf(InputWorkspace="OWS",OutputWorkspace="OWS")													#Apply k_i/k_f factor
   Rebin(InputWorkspace="OWS",OutputWorkspace="OWS",Params=Erange,PreserveEvents=False)						#Make sure the bins are correct
   #ConvertToDistribution(Workspace="OWS") 		
   CloneWorkspace(InputWorkspace="OWS",OutputWorkspace=wkspcnm)
   CloneWorkspace(InputWorkspace="MonWS_sum",OutputWorkspace=wkspcnm+'_mon') 
   DeleteWorkspace(Workspace="OWS")
   DeleteWorkspace(Workspace="MonWS_sum")
   return [Efixed,T0]
   
   
def autoloader(filename,wkspcnm):
     """
     function for autoloader
     """
     LoadEventNexus(Filename=filename,OutputWorkspace="IWS")												#Load an event Nexus file
     FilterBadPulses(InputWorkspace="IWS",OutputWorkspace = "IWS",LowerCutoff = 50)
     LoadNexusMonitors(Filename=filename,OutputWorkspace="MonWS")
     Mons=mtd['IWS']		
     Eguess=array(Mons.getRun()['EnergyRequest'].value).mean()
     runnum=str(Mons.getRun()['run_number'].value)
     [Efixed,T0]=GetEiT0("MonWS",Eguess)
     Estep=Efixed*0.02
     Elow=-Efixed*0.5
     Ehigh=Efixed*0.95
     Erange='%g,%g,%g'%(Elow,Estep,Ehigh)
     in_proc_1("IWS",Erange,wkspcnm,Efixed,T0)
     return [Efixed,T0,runnum]

def in_proc_1(IWS,Erange,wkspcnm,Efixed,T0):
    """
    normalize by current, correct for Tube efficiency, compensate for ki/kf
    """
    ChangeBinOffset(InputWorkspace=IWS,OutputWorkspace="OWS",Offset=T0)
    NormaliseByCurrent(InputWorkspace="OWS",OutputWorkspace="OWS")												#normalize by proton charge
    ConvertUnits(InputWorkspace="OWS",OutputWorkspace="OWS",Target="Wavelength",EMode="Direct",EFixed=Efixed)	#The algorithm for He3 tube efficiency requires wavelength units
    He3TubeEfficiency(InputWorkspace="OWS",OutputWorkspace="OWS")												#Apply correction due to absorption in He3
    ConvertUnits(InputWorkspace="OWS",OutputWorkspace="OWS",Target="DeltaE",EMode="Direct",EFixed=Efixed)		#Switch  to energy transfer
    CorrectKiKf(InputWorkspace="OWS",OutputWorkspace="OWS")													#Apply k_i/k_f factor
    Rebin(InputWorkspace="OWS",OutputWorkspace="OWS",Params=Erange,PreserveEvents=False)						#Make sure the bins are correct
    #ConvertToDistribution(Workspace="OWS") 		
    CloneWorkspace(InputWorkspace="OWS",OutputWorkspace=wkspcnm)
    DeleteWorkspace(Workspace="OWS")

#---------------------------------- THINGS TO CHANGE---------------------------------------------------------

autoreduceflag= True  # if this flag is true,  1 file will be reduced following the auto reduce options.
#inmaskfile='/SNS/SEQ/shared/2011_B/mask_40.nxs'
#dat_runs=[[17624],[2],[3]]     This command will take runs 17624, 2 and 3 and reduce them individually.  

#dat_runs=expand_dims(arange(18141,18240),2)   If all runs are in their own inner square brackets, they will not be summed together. 
# If a series of runs are in individual square brackets, they will be added together in reduction phase e.g. [[2,3]] means 2 and 3 will be added together and then reduced.
dat_runs=[[19589],[19599],[19665]]
#dat_runs = expand_dims(arange(19637,19670),2);


autoreducefile=''
inmaskfile =''# '/SNS/SEQ/shared/masks/mask.nx5'

filename = sys.argv[1]
outdir = sys.argv[2]

eventFile = os.path.split(filename)[-1]
datadir = filename.replace(eventFile, '')
#runNumber = eventFile.split('_')[1]

#Eguess=50.	determined in Load nsum											#initial energy guess
#Erange="-30,0.2,48"																						#Energy bins:    Emin,Estep,Emax
#datadir="/SNS/SEQ/IPTS-6292/data/"#Data directory	
V_dir = "/SNS/SEQ/shared/2012_A/V_files/"
#outdir="/SNS/SEQ/IPTS-6292/shared/fixed/"
outpre="Auto_reduced"														#Output directory
pow_flag=True
ang_name="CCR13VRot"
NXSPE_flag=True
SPE_flag=False
MaskBTP(Instrument='SEQUOIA',Bank="58,62,98,99,100,101,102,118,128,141")
MaskBTP(Workspace="temporaryWorkspaceForMasking",Pixel="1,2,3,4,5,6,7,8,121,122,123,124,125,126,127,128")
MaskBTP(Workspace="temporaryWorkspaceForMasking",Bank="59",Tube="6")
MaskBTP(Workspace="temporaryWorkspaceForMasking",Bank="61",Tube="4")
MaskBTP(Workspace="temporaryWorkspaceForMasking",Bank="96",Tube="8")
MaskBTP(Workspace="temporaryWorkspaceForMasking",Bank="136",Tube="4")
CropWorkspace(InputWorkspace="temporaryWorkspaceForMasking",OutputWorkspace="temporaryWorkspaceForMasking",StartWorkspaceIndex=2)

V_norm = V_norm_obj(V_dir+'SEQ_15884_event.nxs',"V_norm","0.3,0.9,1.2",outdir,maskfile=inmaskfile,maskw="temporaryWorkspaceForMasking",ld_saved_fl=True)
V_norm.CreateMasksAndVanadiumNormalization()
run_info=[]
if autoreduceflag:
  [Ef,T0,runnum]=autoloader(filename,"data_w")
  Divide(LHSWorkspace="data_w", RHSWorkspace="V_norm",OutputWorkspace="data_w")
  ConvertToDistribution("data_w")
  outfile=outpre+'_'+runnum
  SaveNXSPE(InputWorkspace="data_w", Filename= outdir+outfile+".nxspe",Efixed=Ef,psi=str(0),KiOverKfScaling=True)
else:  
  for idx, inruns in enumerate(dat_runs):
      [Ef,T0]=load_n_sum(inruns,datadir,"data_w")
      Divide(LHSWorkspace="data_w", RHSWorkspace="V_norm",OutputWorkspace="data_w")
      if pow_flag:
	    ang=str(0)
      else:
           ang=str(mean(mtd['data_w'].run().get(ang_name).value))
	   angstr=str(ang).replace('.','p')
      outfile=outpre+'_'+str(Eguess).replace('.','p')+'_'+angstr+'_'+str(inruns[0])
      ConvertToDistribution("data_w")
      if NXSPE_flag:
	   SaveNXSPE(InputWorkspace="data_w", Filename= outdir+outfile+".nxspe",Efixed=Ef,psi=ang,KiOverKfScaling=True) 
      if SPE_flag:
             SaveSPE(InputWorkspace="data_w",Filename=outdir+outfile+".spe")
	     temp_string = outdir+outfile+".spe"
	     rn_infostrng='%s,%g,%g,%s\n'%(temp_string,Ef,T0,ang)
	     run_info.append(rn_infostrng)
  if SPE_flag:
	fhand=open(outdir+'Gd_300K_50meV_1.csv','w')
	fhand.writelines(['runnum,Efixed,T0,psi\n','n,meV,mus,deg.\n'])
	fhand.writelines(run_info)
	fhand.close()	

