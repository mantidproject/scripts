from qtiGenie import *
from PySlice2 import *
import time
#from genie_init import *


LoadEmptyInstrument(Filename=r'/Users/jon/Desktop/MARI_Definition.xml',OutputWorkspace='emptyMari',DetectorValue='0',MonitorValue='0')

def spectrumToMantid(specNum):
	ExtractSingleSpectrum(InputWorkspace='emptyMari',OutputWorkspace='singleMari',WorkspaceIndex=str(specnum-1))
	nspec=1
	tmp=get_spectrum(specNum)
	x=tmp[0]
	dat=tmp[1]
	err=sqrt(dat)
	wkspOut=CreateWorkspace(x,dat,err,nspec,'TOF',Distribution=False,ParentWorkspace='singleMari')
	ConvertUnits(InputWorkspace='wkspOut',OutputWorkspace='wkspOut',Target='Wavelength')


