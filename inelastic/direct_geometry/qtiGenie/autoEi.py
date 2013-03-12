# takes wkspin as string of run number 
# returns ei, as value and rebin parameners as string
import numpy as np
from peakdet import *
from mantid.simpleapi import *


def autoEi(WkspIn):
	#monitor spectrum in spectrum number
	monspec='2'
	#load single spectrum for monitor
	Load(Filename=WkspIn,OutputWorkspace='tmp_Monitors',Cache=r'Never',LoadLogFiles='0',SpectrumMin=monspec,SpectrumMax=monspec)
	ConvertToDistribution(Workspace='tmp_Monitors')
	#convert to energy & rebin to make life simple later the rebin is required for mari as monitor 2 always has a sit load of gammas at low tofs
	#these need to be ignored rebin can be removed for monitor 3 but the peak find is less reliable as there is always more noise in the m3
	ConvertUnits(InputWorkspace='tmp_Monitors',OutputWorkspace='tmp_Monitors',Target='Energy')
	Rebin(InputWorkspace='tmp_Monitors',OutputWorkspace='tmp_Monitors',Params='1,1,2000',PreserveEvents='0')
	
	#extract x and y data from specrum & delete mantid wksp goto point data in x
	dat=mtd['tmp_Monitors']
	x=dat.extractX()		
	y=dat.extractY()
	DeleteWorkspace('tmp_Monitors')
	y=y[0]
	x=x[0]

	xx=(x[1:len(x)]+x[0:len(x)-1])/2

	[maximal,miniaml]=peakdet(y, y[0]/2 ,xx)
	#find the biggest peak in the returned list of max data
	max1=maximal[:,1].argmax()
	#get the correspoiding energy
	ei=maximal[max1,0]
	print ei
	#generate some rebin parameters
	rebin=str(-ei*.5)+','+str(ei*(2.5e-3))+','+str(ei*.95)
	#print rebin
	#print type(rebin)
	return ei,rebin
	

#ei,rebin=autoEi('18322')
