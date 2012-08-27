import os
import sys
import shutil

#mantid_root = "/opt/Mantid"
mantid_root = "/SNS/users/scu/build/mantid/adara-demo"
mantid_bin = sys.path.append(os.path.join(mantid_root, "bin"))

from mantid.simpleapi import *

nexus_file=sys.argv[1]
output_directory=sys.argv[2]

## For testing
#nexus_file="/SNS/HYSA/IPTS-8072/nexus/HYSA_67.nxs.h5"
#output_directory="/SNS/HYSA/shared/autoreduce/testoutput/"

filename = os.path.split(nexus_file)[-1]
#run_number = filename.split('_')[1]
run_number = os.path.splitext(os.path.splitext(filename.split('_')[1])[0])[0]

autows = "__auto_ws"
backgroundws = "__background_ws"
tibws = "__tib_ws"

processed_filename = os.path.join(output_directory, "HYSA_" + run_number + "_spe.nxs")
nxspe_filename=os.path.join(output_directory, "HYSA_" + run_number + ".nxspe")

# Load the data
LoadEventNexus(Filename=nexus_file, OutputWorkspace=autows)

# Get Ei from file
Ei=mtd[autows].getRun()['EnergyRequest'].value[0]
#Tzero = 25.0 + 85.0/(1+(Ei/27.0)**4)
Tzero = -1.0*(4.0 + 107.0/(1.0+(Ei/31.0)**3))
#print "Ei=",Ei
#print "Tzero=",Tzero

# Work out some energy bins
emin = -(2.0*Ei)
emax = Ei*0.9
estep = (3.0*Ei)/250
energy_bins = "%f,%f,%f" % (emin,estep,emax)

ChangeBinOffset(InputWorkspace=autows,OutputWorkspace=autows,Offset=Tzero)
Rebin(InputWorkspace=autows,OutputWorkspace=backgroundws,Params='6000,2000,8000')
ConvertUnits(InputWorkspace=autows,OutputWorkspace=autows,Target='DeltaE',EMode='Direct',EFixed=Ei)
Rebin(InputWorkspace=autows,OutputWorkspace=autows,Params=energy_bins,PreserveEvents='0')
ConvertToDistribution(Workspace=autows)
FlatBackground(InputWorkspace=backgroundws,OutputWorkspace=tibws,StartX='6000',EndX='8000',Mode='Mean',OutputMode='Return Background')
ConvertToDistribution(Workspace=tibws)
Minus(LHSWorkspace=autows,RHSWorkspace=tibws,OutputWorkspace=autows)
ConvertFromDistribution(Workspace=autows)
NormaliseByCurrent(InputWorkspace=autows,OutputWorkspace=autows)
ConvertUnits(InputWorkspace=autows,OutputWorkspace=autows,Target='DeltaE',EMode='Direct')
Rebin(InputWorkspace=autows,OutputWorkspace=autows,Params=energy_bins,PreserveEvents='0')
ConvertUnits(InputWorkspace=autows,OutputWorkspace=autows,Target='Wavelength',EMode='Direct',EFixed=Ei)
He3TubeEfficiency(InputWorkspace=autows,OutputWorkspace=autows)
ConvertUnits(InputWorkspace=autows,OutputWorkspace=autows,Target='DeltaE',EMode='Direct',EFixed=Ei)
CorrectKiKf(InputWorkspace=autows,OutputWorkspace=autows)
Rebin(InputWorkspace=autows,OutputWorkspace=autows,Params=energy_bins)
GroupDetectors(InputWorkspace=autows,OutputWorkspace=autows,MapFile=r'/SNS/HYSA/shared/autoreduce/hys_group4.map',Behaviour='Average')
ConvertToDistribution(Workspace=autows)
# Save a file
SaveNexus(Filename=processed_filename, InputWorkspace=autows)
SaveNXSPE(Filename=nxspe_filename, InputWorkspace=autows)



