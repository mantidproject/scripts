import os
import sys
import shutil

mantid_root = "/opt/mantidnightly"
mantid_bin = sys.path.append(os.path.join(mantid_root, "bin"))

from mantid.simpleapi import *
#from MantidFramework import mtd
#mtd.initialise()

nexus_file=sys.argv[1]
output_directory=sys.argv[2]

filename = os.path.split(nexus_file)[-1]
run_number = filename.split('_')[1]

autows = "__auto_ws"
autows_monitor = autows + "_monitor"


dave_grp_filename = os.path.join(output_directory, "BASIS_" + run_number + "_1run.dat")
processed_filename = os.path.join(output_directory, "bss" + run_number + "_silicon111_sqw.nxs")

Load(Filename=nexus_file, OutputWorkspace=autows)
LoadMask(Instrument='BASIS', OutputWorkspace='BASIS_MASK', InputFile='/SNS/BSS/shared/autoreduce/BASIS_Mask.xml')
MaskDetectors(Workspace=autows, MaskedWorkspace='BASIS_MASK')
ModeratorTzero(InputWorkspace=autows,OutputWorkspace=autows)
LoadParameterFile(Workspace=autows, Filename=os.path.join(mantid_root, 'instrument', 'BASIS_silicon_111_Parameters.xml'))
LoadNexusMonitors(Filename=nexus_file, OutputWorkspace=autows_monitor)
Rebin(InputWorkspace=autows_monitor,OutputWorkspace=autows_monitor,Params='10')
ConvertUnits(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor, Target='Wavelength')
OneMinusExponentialCor(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor, C='0.20749999999999999', C1='0.001276')
Scale(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor, Factor='9.9999999999999995e-07')
ConvertUnits(InputWorkspace=autows, OutputWorkspace=autows, Target='Wavelength', EMode='Indirect')
RebinToWorkspace(WorkspaceToRebin=autows, WorkspaceToMatch=autows_monitor, OutputWorkspace=autows)
Divide(LHSWorkspace=autows, RHSWorkspace=autows_monitor,  OutputWorkspace=autows)
ConvertUnits(InputWorkspace=autows, OutputWorkspace=autows, Target='DeltaE', EMode='Indirect')
CorrectKiKf(InputWorkspace=autows, OutputWorkspace=autows,EMode='Indirect')

#RenameWorkspace(InputWorkspace=autows,OutputWorkspace='bss20200_silicon111_red')
Rebin(InputWorkspace=autows, OutputWorkspace=autows, Params='-0.12,0.0004,0.12')
#GroupDetectors(InputWorkspace=autows, OutputWorkspace=autows, MapFile='/SNS/BSS/shared/autoreduce/BASIS_Grouping.xml', Behaviour='Sum')
SofQW3(InputWorkspace=autows, OutputWorkspace=autows+'_sqw', QAxisBinning='0.2,0.2,2.0', EMode='Indirect', EFixed='2.082')
SaveDaveGrp(Filename=dave_grp_filename, InputWorkspace=autows+'_sqw', ToMicroEV=True)
SaveNexus(Filename=processed_filename, InputWorkspace=autows+'_sqw') 
