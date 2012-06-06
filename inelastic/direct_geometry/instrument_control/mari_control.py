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
		self.fix=False
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
	def changetitle(self,title):
		change(Title=title)

	def settemp(self,Temp=0.0,offset=2.0,runControlOffset=2.0):
		
		settemphead=Temp-offset
		
		lowlimittemphead=settemphead-runControlOffset
		highlimittemphead=settemphead+runControlOffset
		
		lowlimittempsample=Temp-runControlOffset
		highlimittempsample=Temp+runControlOffset
		cset(T_head=settemphead)
		cset(T_sample=Temp,runcontrol=True,lowlimit=lowlimittempsample, highlimit=highlimittempsample)
		print 'Setting ',self.cryoblock,'to ',Temp-offset,'K'
		print 'Setting ',self.sampleblock,'to ',Temp,'K'
		print 'Setting runcontol to +/-',runControlOffset,'K'
		print 'Waiting for ',self.cryoblock, ' to reach ',Temp-offset,'+/-',runControlOffset,'K'
		
		cset(T_head=settemphead, wait=True,lowlimit=lowlimittemphead, highlimit=highlimittemphead)
		
		
	def updatestore(self):
		update()
		snapshot_crpt(self.tmpfile)#update from dae and save file as mar00000.raw
		print 'current run saved to ',self.tmpfile
		
	def countStats(self,statlevel=0,updateTime=2,uamps=0):
	
		totalcurrent=uamps
		currentUamps=0
		currentstatistic=1e10
		starttime=time.time()
		begin()
		eltime=[]
		stats=[]
		while currentstatistic > statlevel:
			

			time.sleep(updateTime)
			currentUamps=get_uamps()
			if currentUamps > 0:
				update()
				snapshot_crpt(self.tmpfile)#update from dae and save file as mar00000.raw
				cal_file=self.whitebeamfile
				LoadRaw(Filename=str(self.whitebeamfile),OutputWorkspace="wb_wksp",LoadLogFiles="0")
		
				qbin=str(self.rebinInput[0])+','+str(self.rebinInput[1]-self.rebinInput[0])+','+str(self.rebinInput[1])
				ebin=str(self.rebinInput[2])+','+str(self.rebinInput[3]-self.rebinInput[2])+','+str(self.rebinInput[3])
	
				LoadRaw(Filename=str(self.runfile),OutputWorkspace="run_wksp",LoadLogFiles="0")
				if self.fix==True:
					data=iliad("wb_wksp","run_wksp",self.ei,ebin,self.mapfile,det_cal_file=cal_file,norm_method='current',fixei=True)
				else:
					data=iliad("wb_wksp","run_wksp",self.ei,ebin,self.mapfile,det_cal_file=cal_file,norm_method='current')
	
				w2=data2D('data')
				w2.rebinproj(qbin)
				currentstatistic=w2.integrate(self.rebinInput[0],self.rebinInput[1],self.rebinInput[2],self.rebinInput[3])
				currentstatistic=w2.percentError(self.rebinInput[0],self.rebinInput[1],self.rebinInput[2],self.rebinInput[3])
				elaspedtime=time.time()-starttime
				eltime.append(elaspedtime/60)
				stats.append(currentstatistic)
				statsWksp=CreateWorkspace(eltime,stats,stats,1)
				print 'ROI statistic = ',str(currentstatistic),' %: Runtime= ',int(elaspedtime),' secs Integrated current = ',currentUamps,' Uamp'
				halfstat=int(elaspedtime*4)/3600
				print 'Runtime to halve current statistic on ROI =',halfstat,'hours'
				if uamps>0 and currentUamps>uamps:
					print 'Integrated current limit reached before statistic limit= ',currentUamps,' Runtime= ',int(elaspedtime),' secs' 
					break
			else:
				print 'No beam since start of run'

				
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

	def seteifixed(self):
		self.fix=True
		print 'Ei will be fixed to ', self.ei,'meV'
		
	def seteifree(self):
		self.fix=False
		print 'Ei will be calculated'
		
	def setEiForControl(self,ei):
		self.ei=0.0
		self.ei=ei
		print 'incident energy for control set to ',self.ei,'meV'

#init the class once and define some alias commands for the class
run=mariControl()

defineRoi=run.defineRoi
setWhiteBeam=run.setWhiteBeam
seteiforcontrol=run.setEiForControl
count=run.count
countstats=run.countStats
settemp=run.settemp
changetitle=run.changetitle
updatestore=run.updatestore
seteifixed=run.seteifixed
seteifree=run.seteifree
#defineRoi(5.0,5.5,50.0,60.0)
#run.runfile='16654' # this would be the hard coded way of setting the runfile for reduction
#seteiforcontrol(100)
#run.count(1,2)
#countstats(4,20,5)
#abort()