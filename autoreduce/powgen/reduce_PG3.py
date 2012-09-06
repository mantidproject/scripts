import os
import sys
import shutil 
sys.path.append("/opt/Mantid/bin")
from MantidFramework import mtd
mtd.initialize()
from mantidsimple import *

cal_dir = "/SNS/PG3/2012_2_11A_CAL/"
cal_file  = os.path.join(cal_dir, "PG3_FERNS_d10805_2012_08_29.cal")
char_file = os.path.join(cal_dir, "PG3_characterization_2012_08_27-HR.txt")
#MODE = 0664

#from mantidsimple import *

eventFileAbs=sys.argv[1]
outputDir=sys.argv[2]

eventFile = os.path.split(eventFileAbs)[-1]
nexusDir = eventFileAbs.replace(eventFile, '')
runNumber = eventFile.split('_')[1]
configService = mtd.getSettings()
dataSearchPath = configService.getDataSearchDirs()
dataSearchPath.append(nexusDir)
configService.setDataSearchDirs(dataSearchPath)

SNSPowderReduction(Instrument="PG3", RunNumber=runNumber, Extension="_event.nxs",
                   PreserveEvents=True,PushDataPositive="AddMinimum",
                   CalibrationFile=cal_file, CharacterizationRunsFile=char_file,
                   LowResRef=15000, RemovePromptPulseWidth=50,
                   Binning=-0.0008, BinInDspace=True, FilterBadPulses=True,
                   SaveAs="gsas and fullprof", OutputDirectory=outputDir,
                   NormalizeByCurrent=True, FinalDataUnits="dSpacing")


#dirList=os.listdir(outputDir)
#for fname in dirList:
#  os.chmod(os.path.join(outputdir, fname), MODE)

#fileName= "PG3_" + runNumber + "_REDUCED.gsa"
#outputFile=outputDir+"/"+fileName
#f=open(outputFile, 'w')
#f.write('POWGEN auto data reduction results')
#
#fileName2= "PG3_" + runNumber + "_REDUCED2.gsa"
#outputFile2=outputDir+"/"+fileName2
#f2=open(outputFile2, 'w')
#f2.write('More POWGEN auto data reduction results')
#
#print outputFile + " is created"
#print outputFile2 + " is created"
