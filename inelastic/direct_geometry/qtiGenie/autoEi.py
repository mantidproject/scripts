# takes wkspin as string of run number 
# returns ei, as value and rebin parameners as string
import numpy as np
from peakdet import *
from mantid.simpleapi import *


def autoEi(WkspIn):
	monspec='2'
	#	
	#WkspIn='18315'

	Load(Filename=WkspIn,OutputWorkspace='tmp_Monitors',Cache=r'Never',LoadLogFiles='0',SpectrumMin=monspec,SpectrumMax=monspec)

	ConvertToDistribution(Workspace='tmp_Monitors')
	ConvertUnits(InputWorkspace='tmp_Monitors',OutputWorkspace='tmp_Monitors',Target='Energy')
	Rebin(InputWorkspace='tmp_Monitors',OutputWorkspace='tmp_Monitors',Params='2,1,2000',PreserveEvents='0')
	#FindPeaks(InputWorkspace='tmp_Monitors',PeakFunction='Voigt',BackgroundType='Flat',HighBackground='0',PeaksList='out')
	#ConvertTableToMatrixWorkspace(InputWorkspace='out',OutputWorkspace='out',ColumnX='centre',ColumnY='height',ColumnE='height')

	dat=mtd['tmp_Monitors']
	DeleteWorkspace('tmp_Monitors')

	x=dat.extractX()		
	y=dat.extractY()

	y=y[0]
	x=x[0]

	xx=(x[1:len(x)]+x[0:len(x)-1])/2

	indat=np.zeros((len(y),2))

	indat[:,0]=xx
	indat[:,1]=y


	[maximal,miniaml]=peakdet(y, y[0]/2 ,xx)
	#print maximal.shape

	#dat=mtd['out']


	max1=maximal[:,1].argmax()

	ei=maximal[max1,0]
	print ei

	rebin=str(-ei*.5)+','+str(ei*(2.5e-3))+','+str(ei*.95)
	#print rebin
	#print type(rebin)
	return ei,rebin
	

#ei,rebin=autoEi('18249')
