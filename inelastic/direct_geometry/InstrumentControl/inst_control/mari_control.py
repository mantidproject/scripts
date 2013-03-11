import sys
sys.path.append("C:\\LabVIEW Modules\\dae\\genie_python")

from genie_init import *
import time
from qtiGenie import *
from PySlice2 import *



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
		print get_uamps()
	def reduction_ei(self,*args):
			inp=args[0]
			self.ei=float(inp[0])            
	def list(self):
	#returns a list of displayable command names must be a nicer way todo this!
		return ['set_ei ei freq',
		'count uamps',
		'waitfor uamps',
		'begin',
		'end',
		'settemp temp',
		'changetitle runtitle',
		'cset Command',
		'scanTemp Tstart Tstep Tend uamps Title',
		'setWhiteBeam WBrunNumber',
		'reduction_ei eiForReduction',
		'defineRoi Qmin Qmax Emin Emax',
		'countStats statlevel uamps',
		]
		
	def set_ei(self,*args):
		#ei,freq
		inp=args[0]
		ei=int(inp[0])
		freq=int(inp[1])
		print 'setenergy'
		print 'ei',ei
		print 'freq',freq
		self.ei=0.0
		self.ei=ei
		self.seteifree()
		#set ei here
		set_ei(ei,freq)
		
	def count(self,*args):
	# uamps=0,updateTime=2
		inp=args[0]
		uamps=float(inp[0])
		updateTime=2
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
		
	def changetitle(self,*args):
		title=args[0]
		print title
		settitle=''
		for i in range(len(title)):
			settitle=settitle+' '+title[i]
			
		change(Title=settitle)

	def cset(self,*args):
		command=args[0]
		
		command=str(command[0])
		print command
		commandstr='cset('+command+')'
		eval(commandstr)
	def scanTemp(self,*args):
		inp=args[0]
		print len(inp)
		Tstart=float(inp[0])
		Tstep=float(inp[1])
		Tend=float(inp[2])
		uamps=float(inp[3])
		Title=''
		for i in range(4,len(inp)):
			Title=Title+inp[i]
		for T in range(int(Tstart),int(Tend+Tstep),int(Tstep)):
			settitle=Title+' Temp='+str(T)+'K'
			change(Title=settitle)
			self.settemp([T,0])
			self.count([uamps,0])

			
		
	def settemp(self,*args):
		#Temp=0.0,offset=2.0,runControlOffset=2.0
		inp=args[0]
		Temp=float(inp[0])
		offset=2.0
		runControlOffset=2.0
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
		
		
	def updatestore(self,*args):
		update()
		snapshot_crpt(self.tmpfile)#update from dae and save file as mar00000.raw
		print 'current run saved to ',self.tmpfile
		
	def countStats(self,*args):
		#statlevel=0,uamps=0,updateTime=2
		inp=args[0]
		statlevel=float(inp[0])
		uamps=float(inp[1])
		updateTime=30
		
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
				halfstat=(elaspedtime*4.0)/3600.0
				print 'Runtime to halve current statistic on ROI =',halfstat,'hours'
				if uamps>0 and currentUamps>uamps:
					print 'Integrated current limit reached before statistic limit= ',currentUamps,' Runtime= ',int(elaspedtime),' secs' 
					break
			else:
				print 'No beam since start of run'

				
		end()
	

	def defineRoi(self,*args):
		#Qmin,Qmax,Emin,Emax
		inp=args[0]
		Qmin=float(inp[0])
		Qmax=float(inp[1])
		Emin=float(inp[2])
		Emax=float(inp[3])
		self.rebinInput=[]	
		self.rebinInput.append(Qmin)
		self.rebinInput.append(Qmax)
		self.rebinInput.append(Emin)
		self.rebinInput.append(Emax)
		print 'region of interest for control set to ', self.rebinInput


	def setWhiteBeam(self,*args):
		#whiteBeamRunNumber
		inp=args[0]
		whiteBeamRunNumber=inp[0]
		self.whitebeamfile=''
		self.whitebeamfile=str(whiteBeamRunNumber)
		print 'whitebeam file for control set to ',self.whitebeamfile

	def seteifixed(self,*args):
		self.fix=True
		print 'Ei will be fixed to ', self.ei,'meV'
		
	def seteifree(self,*args):
		self.fix=False
		print 'Ei will be calculated'
    
	def begin(self,*args):
		begin()
            
	def end(self,*args):
		end()
        
	def setEiForControl(self,*args):
		#ei
		inp=args[0]
		ei=float(inp[0])
		self.ei=0.0
		self.ei=ei
		print 'incident energy for control set to ',self.ei,'meV'

#init the class once and define some alias commands for the class
#run=mariControl()

#defineRoi=run.defineRoi
#setWhiteBeam=run.setWhiteBeam
#seteiforcontrol=run.setEiForControl
#count=run.count
#countstats=run.countStats
#settemp=run.settemp
#changetitle=run.changetitle
#updatestore=run.updatestore
#seteifixed=run.seteifixed
#seteifree=run.seteifree
#defineRoi(5.0,5.5,50.0,60.0)
#run.runfile='16654' # this would be the hard coded way of setting the runfile for reduction
#seteiforcontrol(100)
#run.count(1,2)
#countstats(4,20,5)
#abort()