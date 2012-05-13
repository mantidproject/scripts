from qtiGenie import *
from PySlice2 import *
import time
from genie_init import *




class mariControl:
	
	def __init__(self):
		self.whitebeamfile=''
		self.runfile='00000'
		self.ei=0
		self.rebinInput=[]
		self.delQ=0
		self.delE=0
		self.statlevel=2
		self.inst='mar'
		self.mapfile='mari_res'
		self.ext='.raw'
		self.tmpfile='c:\data\mar00000.raw'
		self.cryoblock='T_Head'
		self.sampleblock='T_sample'
		iliad_setup(self.inst)
	
	def count(self,uamps=0,updateTime=2):
	
		totalcurrent=uamps
		currentUamps=0
		starttime=time.time()
		begin()
		while currentUamps<totalcurrent:
			currentUamps=get_uamps()

			time.sleep(updateTime)
		
			elaspedtime=time.time()-starttime
			print 'Integrated current = ',currentUamps,' Runtime= ',int(elaspedtime),' secs'
	
		end()

	def settemp(Temp=0,offset=2,runControlOffset=2)
		cset(self.cryoblock=Temp-offeset)
		cset(self.sampleblock=Temp,runcontrol=True,lowlimit=Temp-runControlOffset, highlimit=Temp+runControlOffset)
		print 'Setting ',self.cryoblock,'to ',Temp-offset,'K'
		print 'Setting ',self.sampleblock,'to ',Temp,'K'
		print 'Setting runcontol to +/-',runControlOffset,'K'
		print 'Waiting for ',self.cryoblock, ' to reach ',Temp-offset,'+/-',runControlOffset,'K'
		
		cset(self.cryoblock=Temp-offeset, wait=True,lowlimit=Temp-runControlOffset, highlimit=Temp+runControlOffset)
		
		
		
	def countStats(self,statlevel=0,updateTime=2,uamps=0):
	
		totalcurrent=uamps
		currentUamps=0
		currentstatistic=1e10
		starttime=time.time()
		begin()
		while currentstatistic > statlevel:
			currentUamps=get_uamps()
			

			

				
			time.sleep(updateTime)
			snapshot_crpt(self.tmpfile)#update from dae and save file as mar00000.raw
			cal_file=self.whitebeamfile
			LoadRaw(Filename=str(self.whitebeamfile),OutputWorkspace="wb_wksp",LoadLogFiles="0")
	
			qbin=str(self.rebinInput[0])+','+str(self.rebinInput[1]-self.rebinInput[0])+','+str(self.rebinInput[1])
			ebin=str(self.rebinInput[2])+','+str(self.rebinInput[3]-self.rebinInput[2])+','+str(self.rebinInput[3])

			LoadRaw(Filename=str(self.runfile),OutputWorkspace="run_wksp",LoadLogFiles="0")	
			data=iliad("wb_wksp","run_wksp",self.ei,ebin,self.mapfile,det_cal_file=cal_file,norm_method='current',fixei=True,)


			w2=data2D('data')
			w2.rebinproj(qbin)
			currentstatistic=w2.integrate(self.rebinInput[0],self.rebinInput[1],self.rebinInput[2],self.rebinInput[3])
			currentstatistic=w2.percentError(self.rebinInput[0],self.rebinInput[1],self.rebinInput[2],self.rebinInput[3])
			elaspedtime=time.time()-starttime
			print 'ROI statistic = ',str(currentstatistic),' %: Runtime= ',int(elaspedtime),' secs Integrated current = ',currentUamps,' Uamp'
			
			if uamps>0 and currentUamps>uamps:
				print 'Integrated current limit reached before statistic limit= ',currentUamps,' Runtime= ',int(elaspedtime),' secs' 
				break

				
		end()
	

	def defineRoi(self,Qmin,Qmax,Emin,Emax):

		self.rebinInput=[]	
		self.rebinInput.append(Qmin)
		self.rebinInput.append(Qmax)
		self.rebinInput.append(Emin)
		self.rebinInput.append(Emax)
		print 'region of interest for control set to ', self.rebinInput

		

	def setWhiteBeam(self,whiteBeamRunNumber):
		self.whitebeamfile=''
		self.whitebeamfile=str(whiteBeamRunNumber)
		print 'whitebeam file for control set to ',self.whitebeamfile

	

	def setEiForControl(self,ei):
		self.ei=0.0
		self.ei=ei
		print 'incident energy for control set to ',self.ei=ei,'meV'

#init the class once and define some alias commands for the class
run=mariControl()

defineRoi=run.defineRoi
setWhiteBeam=run.setWhiteBeam
seteiforcontrol=run.setEiForControl
count=run.count
countstats=run.countStats

settemp=run.settemp

defineRoi(5.0,5.5,50.0,60.0)
#run.runfile='16654' # this would be the hard coded way of setting the runfile for reduction
setEiForControl(100)
run.count(1,2)
run.countStats(4,5,0)