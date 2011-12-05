from string import *
from numpy import arange
from os.path import exists

#Script to reduce a series of data files from within Mantid.
#Will make one .nxspe file for each block of data files in the loop.
#Developed by Andrei Savici and Georg Ehlers.
#Requires the September 2011 release of Mantid (or later) to run.

#This reads blocks of data files.

# first and last run numbers
firstrun=44005
lastrun=44013
runstep=9	# how many files are in one block
num2add=9	# how many files should be added for one output file
runlist=arange(firstrun,lastrun+1,runstep)
# data and output directories
inputdir="/SNS/CNCS/IPTS-5439/data/"
outputdir="/SNS/CNCS/IPTS-5439/shared/NXSPE/"
# want vanadium normalization? say True or False, provide ranges
# (expects a monochromatic or white beam vanadium run)
do_vana=True
vanadiumdir="/SNS/CNCS/shared/Masks_Mantid/35029/NeXus/"
vanadium_filename=vanadiumdir+"CNCS_35029_event.nxs"
vana_lowrange=2000
vana_hirange=18000
# want to subtract the time independent background?
# say True or False, provide ranges
do_tib=True
tibmin=19000
tibmax=20000
tibstep=tibmax-tibmin
tibpar=str(tibmin)+","+str(tibstep)+","+str(tibmax)
# energy used for the runs, desired binning
T0="-22.8"
Efixed="12.07"
bins="-100.0,0.02,12."
# hard mask
maskfile="/SNS/CNCS/shared/Masks_Mantid/CNCS_Mantid_8pixel.mask"
#maskfile="/SNS/CNCS/shared/Masks_Mantid/CNCS_Mantid_8pixel_bothSides.mask"
# sample rotation device name
rotationdevice="SERotator2"
# say if you want ASCII spe & phx output
want_ASCII=False
# say if you want a powder averaged output
do_powder=True
anglemin=0.				#minumum angle
anglemax=140.				#maximum angle
anglestep=0.25				#angle step - this can be fine tuned for pixel arc over detectors

# no need to change anything below here
def GetMask(maskfile):
    f = open(maskfile, 'r')
    masklist = []
    for line in f:
        if line.startswith('-'):
            continue
        parts = line.split()
        masklist.extend(parts)
    return ','.join(masklist)

# Function to load white beam vanadium
def LoadVana(vanadium_filename,vana_lowrange,vana_hirange):
	if (exists(outputdir+"van.nx5")):
		LoadNexus(Filename=outputdir+"van.nx5",OutputWorkspace="vanadium.nxs")
	else:
		vana_step=vana_hirange-vana_lowrange
		LoadEventNexus(Filename=vanadium_filename,OutputWorkspace="vanadium.nxs")
		Rebin(InputWorkspace="vanadium.nxs",OutputWorkspace="vanadium_temp.nxs",Params=[vana_lowrange,vana_step,vana_hirange],PreserveEvents=False)
		NormaliseByCurrent(InputWorkspace="vanadium_temp.nxs",OutputWorkspace="vanadium.nxs")
		MedianDetectorTest(InputWorkspace="vanadium.nxs",OutputWorkspace="vanadium_mask")
		SaveNexus(InputWorkspace="vanadium.nxs",Filename=outputdir+"van.nx5")

#Function to create a map of detectors corresponding to angles in a certain range	
def createanglelist(ws,amin,amax,astep):
	bin_angles=arange(amin+astep*0.5,amax+astep*0.5,astep)
	a=[[] for i in range(len(bin_angles))] #list of list with detector IDs
	w=mtd[ws]
	origin = w.getInstrument().getSample().getPos()
	for i in range(w.getNumberHistograms()):
		ang=w.getDetector(i).getTwoTheta(origin,V3D(0,0,1))*180/math.pi
		index=int((ang-amin)/astep)
		if (index>=0) and (index<len(a)) and ((w.getDetector(i).getID())>0):
			a[index].append(w.getDetector(i).getID())
	#create lists with angles and detector ID only for bins where there are detectors 
	ang_list=[]
	detIDlist=[]
	for elem,ang in zip(a,bin_angles):
		if len(elem)>0:
			detIDlist.append(elem)
			ang_list.append(ang)
	# file with grouping information
	f=open(outputdir+"group.map",'w')
	print >>f,len(ang_list)
	for i in range(len(ang_list)):
		print >>f,i
		print >>f,len(detIDlist[i])
		mystring=str(detIDlist[i]).strip(']').strip('[')
		mystring=mystring.replace(',','')
		print >>f,mystring
	f.close()
	#print ang_list
	# par file
	f=open(outputdir+"group.par",'w')
	print >>f,len(ang_list)
	for i in range(len(ang_list)):
		print >>f,3.5,ang_list[i],0.0,1.0,1.0,1
	f.close()
	return [ang_list,detIDlist]
		
#mapping=createanglelist("vanadium.nxs",anglemin,anglemax,anglestep)
#print GetMask(maskfile)

if (do_vana):
	LoadVana(vanadium_filename,vana_lowrange,vana_hirange)

for run in runlist:
	Input_filename=inputdir+"CNCS_"+str(run)+"_event.nxs"
	print "load file",Input_filename,"to start"
	if (do_powder):
		Output_filename=outputdir+"CNCSpowder_"+str(run)+".nxspe"
	else:
		Output_filename=outputdir+"CNCS_"+str(run)+".nxspe"
	LoadEventNexus(Filename=Input_filename,OutputWorkspace="data.nxs")
	if (num2add>1):
		for nf in range(1,num2add,1):
			Input_filename=inputdir+"CNCS_"+str(run+nf)+"_event.nxs"
			print "load file",Input_filename,"to add"
			LoadEventNexus(Filename=Input_filename,OutputWorkspace="data1.nxs")
			Plus("data.nxs","data1.nxs","data.nxs")
	ChangeBinOffset(InputWorkspace="data.nxs",OutputWorkspace="reduced",Offset=T0)
	Rebin(InputWorkspace="reduced",OutputWorkspace="background_origin_ws",Params=tibpar)
	ConvertUnits(InputWorkspace="reduced",OutputWorkspace="reduced",Target="DeltaE",EMode="Direct",EFixed=Efixed)
	Rebin(InputWorkspace="reduced",OutputWorkspace="temp_reduced",Params=bins,PreserveEvents=False)
	RenameWorkspace(InputWorkspace="temp_reduced",OutputWorkspace="reduced")
	ConvertUnits(InputWorkspace="reduced",OutputWorkspace="reduced",Target="TOF",EMode="Direct",EFixed=Efixed)
	#here tib is subtracted
	if (do_tib):
		ConvertToDistribution(Workspace="reduced")
		FlatBackground(InputWorkspace="background_origin_ws",OutputWorkspace="background_ws",StartX=tibmin,EndX=tibmax,Mode="Mean",OutputMode="Return Background")
		ConvertToDistribution(Workspace="background_ws")
		Minus(LHSWorkspace="reduced",RHSWorkspace="background_ws",OutputWorkspace="reduced")
		ConvertFromDistribution(Workspace="reduced")
	#done with tib
	NormaliseByCurrent(InputWorkspace="reduced",OutputWorkspace="reduced")
	ConvertUnits(InputWorkspace="reduced",OutputWorkspace="reduced",Target="DeltaE",EMode="Direct")
	Rebin(InputWorkspace="reduced",OutputWorkspace="reduced",Params=bins,PreserveEvents=False)
	ConvertUnits(InputWorkspace="reduced",OutputWorkspace="reduced",Target="Wavelength",EMode="Direct",EFixed=Efixed)
	He3TubeEfficiency(InputWorkspace="reduced",OutputWorkspace="reduced")
	ConvertUnits(InputWorkspace="reduced",OutputWorkspace="reduced",Target="DeltaE",EMode="Direct",EFixed=Efixed)
	CorrectKiKf(InputWorkspace="reduced",OutputWorkspace="reduced")
	Rebin(InputWorkspace="reduced",OutputWorkspace="reduced",Params=bins)
	ConvertToDistribution(Workspace="reduced")
	MaskDetectors(Workspace="reduced",DetectorList=GetMask(maskfile))
	# normalize by vanadium
	if (do_vana):
		MaskDetectors(Workspace="vanadium.nxs",DetectorList=GetMask(maskfile))
		#MaskDetectors(Workspace="vanadium.nxs",MaskedWorkspace="vanadium_mask")
		Divide(LHSWorkspace="reduced",RHSWorkspace="vanadium.nxs",OutputWorkspace="reduced")
		MaskDetectors(Workspace="reduced",DetectorList=GetMask(maskfile))
		#MaskDetectors(Workspace="reduced",MaskedWorkspace="vanadium_mask")
	# powder average
	if (do_powder):
		if (run==firstrun):
			mapping=createanglelist("reduced",anglemin,anglemax,anglestep)
			print "number of powder spectra",len(mapping[0])
		GroupDetectors(InputWorkspace="reduced",OutputWorkspace="reduced",MapFile=outputdir+"group.map",Behaviour="Sum")
		SolidAngle(InputWorkspace="reduced",OutputWorkspace="sa")
		Divide(LHSWorkspace="reduced",RHSWorkspace="sa",OutputWorkspace="reduced")
	# finish up and save to file
	w=mtd["reduced"]
	psi=w.getRun()[rotationdevice].value
	psi=psi[0]
	print "save file",Output_filename
	if (do_powder):
		SaveNXSPE(InputWorkspace="reduced",Filename=Output_filename,Efixed=Efixed,psi=psi,ki_over_kf_scaling=True,ParFile=outputdir+"group.par")
	else:
		SaveNXSPE(InputWorkspace="reduced",Filename=Output_filename,Efixed=Efixed,psi=psi,ki_over_kf_scaling=True)
	if (want_ASCII):
		Output_spe_filename=outputdir+"CNCS_"+str(run)+".spe"
		SaveSPE(InputWorkspace="reduced",Filename=Output_spe_filename)
		if (run==firstrun):
			Output_phx_filename=outputdir+"CNCS_"+str(run)+".phx"
			SavePHX(InputWorkspace="reduced",Filename=Output_phx_filename)
		
