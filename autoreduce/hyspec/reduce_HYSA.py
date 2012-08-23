import os
import sys
import shutil

mantid_root = "/opt/Mantid"
mantid_bin = sys.path.append(os.path.join(mantid_root, "bin"))

from mantid.simpleapi import *

nexus_file=sys.argv[1]
output_directory=sys.argv[2]

filename = os.path.split(nexus_file)[-1]
run_number = filename.split('_')[1]

autows = "__auto_ws"
backgroundws = "__background_ws"
tibws = "__tib_ws"

# TODO: Get this value from the NeXus File.
Ei = 70.0
Tzero = -26.840661475400001
energy_bins = "-100,0.5,60"

processed_filename = os.path.join(output_directory, "HYSA_" + run_number + "_spe.nxs")

LoadEventNexus(Filename=nexus_file,OutputWorkspace=autows)
ChangeBinOffset(InputWorkspace=autows,OutputWorkspace=autows,Offset=Tzero)
Rebin(InputWorkspace=autows,OutputWorkspace=backgroundws,Params='6000,2000,8000')
ConvertUnits(InputWorkspace=autows,OutputWorkspace=autows,Target='DeltaE',EMode='Direct',EFixed=Ei)
Rebin(InputWorkspace=autows,OutputWorkspace=autows,Params=energy_bins,PreserveEvents='0')
ConvertToDistribution(Workspace=autows)
FlatBackground(InputWorkspace=background_ws,OutputWorkspace=tibws,StartX='6000',EndX='8000',Mode='Mean',OutputMode='Return Background')
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
GroupDetectors(InputWorkspace=autows,OutputWorkspace=autows,MapFile=r'/SNS/HYS/shared/hys_group/hys_group4.map',Behaviour='Average')
ConvertToDistribution(Workspace=autows)
# Save a file
SaveNexus(Filename=processed_filename, InputWorkspace=autows) 




