from utils import *
#from mantidsimple import *
import MantidFramework 
MantidFramework.mtd.initialise()
#from DirectEnergyConversion import *
import time as time
from mantidplotpy import *
import dgreduce
import inspect
import numpy
from mantidplot import *
from mantid import *
from mantid.simpleapi import *
import pprint
from numpy import *
#import qti



class data2D:
	global active_layer, graph,figure_dict,fignum
	def __init__(self,wksp):
		self.data=0
		self.data_disp=0
		self.init_params()
		self.wksp_name=mtd[wksp] #new api means this is actually the workspace
	def init_params(self):
		self.linecol=2
		self.figure_dict={}
		self.fignum=0

	def rebinproj(self,qbin,**kwargs):
		#two workspaces one for the cutting and one for the 2D display
		wksp_in=self.wksp_name
		wksp_new_name=wksp_in.getName()+'_SQW'
		wksp_newDisp=wksp_new_name+'display'
		if kwargs.has_key('fast'):
			print 'Using fast rebin'
			tmpdat=self.sqwfast(self.wksp_name,qbin)
			Transpose(InputWorkspace='tmpdat',OutputWorkspace=wksp_newDisp)
			RenameWorkspace(InputWorkspace=tmpdat,OutputWorkspace=wksp_new_name)
			
			self.data=mtd[wksp_new_name]
			self.data_disp=mtd[wksp_newDisp]
		else:
			print 'Using intersecting area rebin'
			tmpdat=self.sqw(self.wksp_name,qbin)
			
			Transpose(InputWorkspace='tmpdat',OutputWorkspace=wksp_newDisp)
			RenameWorkspace(InputWorkspace='tmpdat',OutputWorkspace=wksp_new_name)
			
			self.data=mtd[wksp_new_name]
			self.data_disp=mtd[wksp_newDisp]
	
	def smooth(self,*args):
		"""
		smooth(n, zscale)
		smooths the display workspace with a n times a Gaussian blur if scipy is available or the smoothdata alg 
		which is n times adjacent averaging. 
		"""
		n=args[0]
		
		try:
			from DgFuncs import gaussSmooth
			name=self.wksp_name.getName()+'_SQWdisplay'#self.data_disp.getName()
			ww=gaussSmooth(name,n)
			RenameWorkspace(InputWorkspace=ww,OutputWorkspace=name)
			self.data_disp=mtd[name]
			if len(args)==1:
				self.display()
			if len(args)==2:
				self.display(args[1])
				
		except:
			if n<=2:
				n=3
			print "Failure of Gaussian smooth revert to averaging"	
			SmoothData(InputWorkspace=self.data_disp,OutputWorkspace=self.data_disp,NPoints=n)
			if len(args)==1:
				self.display()
			if len(args)==2:
				self.display(args[1])
	
	def boseFac(self,T):
		from DgFuncs import bose
		
		name=self.wksp_name.getName()+'_SQWdisplay'#self.data_disp.getName()
		name2=self.wksp_name.getName()+'_SQW'#self.data.getName()
		#dat=mtd[name]
		Transpose(InputWorkspace=name2,OutputWorkspace=name2)
		ww=bose(name2,T)
		RenameWorkspace(InputWorkspace=ww,OutputWorkspace=name)
		Transpose(InputWorkspace=ww,OutputWorkspace=name2)
		self.data_disp=mtd[name]
		self.data=mtd[name2]
	
	def units(self,unit):
		from DgFuncs import perCm
		from DgFuncs import meV
		
		if unit =='cm':
			name=self.wksp_name.getName()+'_SQWdisplay'#self.data_disp.getName()
			name2=self.wksp_name.getName()+'_SQW'#self.data.getName()
			#dat=mtd[name]
			Transpose(InputWorkspace=name2,OutputWorkspace=name2)
			ww=perCm(name2)
			RenameWorkspace(InputWorkspace=ww,OutputWorkspace=name)
			Transpose(InputWorkspace=ww,OutputWorkspace=name2)
		
		if unit =='mev':
			name=self.wksp_name.getName()+'_SQWdisplay'#self.data_disp.getName()
			name2=self.wksp_name.getName()+'_SQW'#self.data.getName()
			#dat=mtd[name]
			Transpose(InputWorkspace=name2,OutputWorkspace=name2)
			ww=meV(name2)
			RenameWorkspace(InputWorkspace=ww,OutputWorkspace=name)
			Transpose(InputWorkspace=ww,OutputWorkspace=name2)
			
		self.data_disp=mtd[name]
		self.data=mtd[name2]
	
	def gofe(self,T,dbwFac):
		from DgFuncs import gofe
		
		name=self.wksp_name.getName()+'_SQWdisplay'#self.data_disp.getName()
		name2=self.wksp_name.getName()+'_SQW'#self.data.getName()
		#dat=mtd[name]
		Transpose(InputWorkspace=name2,OutputWorkspace=name2)
		ww=gofe(name2,T,dbwFac)
		RenameWorkspace(InputWorkspace=ww,OutputWorkspace=name)
		Transpose(InputWorkspace=ww,OutputWorkspace=name2)
		self.data_disp=mtd[name]
		self.data=mtd[name2]
	
	def byE(self,T,dbwFac):
		from DgFuncs import gofe
		
		name=self.wksp_name.getName()+'_SQWdisplay'#self.data_disp.getName()
		name2=self.wksp_name.getName()+'_SQW'#self.data.getName()
		#dat=mtd[name]
		Transpose(InputWorkspace=name2,OutputWorkspace=name2)
		ww=byE(name2)
		RenameWorkspace(InputWorkspace=ww,OutputWorkspace=name)
		Transpose(InputWorkspace=ww,OutputWorkspace=name2)
		self.data_disp=mtd[name]
		self.data=mtd[name2]
		
	def colour(self,col):
		active_layer.setCurveLineColor(1,col)
	
	def ECut(self,Qmin,Qmax,Emin,delE,Emax,**kwargs):
		if kwargs.has_key('over'):
			self.cut(self.data,Qmin,Qmax,Emin,delE,Emax,along='e',over=True)
		else:
			self.cut(self.data,Qmin,Qmax,Emin,delE,Emax,along='e')
	
	def QCut(self,Emin,Emax,Qmin,delQ,Qmax,**kwargs):
		if kwargs.has_key('over'):
			self.cut(self.data,Emin,Emax,Qmin,delQ,Qmax,along='q',over=True)
		else:
			self.cut(self.data,Emin,Emax,Qmin,delQ,Qmax,along='q')
			
	
	def xLims(self,*args):
		if len(args) == 2:
			min=args[0]
			max=args[1]
			index=len(self.figure_dict)
		if len(args) ==3:
			min=args[0]
			max=args[1]
			index=args[2]
		lgr=self.figure_dict.get(index)[1]
		active_layer = lgr.activeLayer()
		active_layer.setScale(2,min,max)#set xscale
		
	
	def yLims(self,*args):
		
		
		if len(args) == 2:
			min=args[0]
			max=args[1]
			index=len(self.figure_dict)
		if len(args) ==3:
			min=args[0]
			max=args[1]
			index=args[2]
		
		lgr=self.figure_dict.get(index)[1]
		active_layer = lgr.activeLayer()
		active_layer.setScale(0,min,max)#set yscale
	
	def display(self,*args):
		data=self.data_disp	
		print type(data)
		print self.wksp_name.getName()
		#dat=self.data.importMatrixWorkspace().plotGraph2D()
		#out=importMatrixWorkspace(data.getName())
		out=importMatrixWorkspace(self.wksp_name.getName()+'_SQWdisplay')
		g2D=out.plotGraph2D()
		#set z scale 
		if len(args)==1:
			ll=g2D.activeLayer()
			ll.setAxisScale(Layer.Right,0,args[0])
		if len(args)==2:
			ll=g2D.activeLayer()
			ll.setAxisScale(Layer.Right,args[0],args[1])
	def figures(self):
		
		pprint.pprint(self.figure_dict)
	
	def mergeFigs(self,index1,index2):
		lgr1=self.figure_dict.get(index1)[1]
		lgr2=self.figure_dict.get(index2)[1]
		mergePlots(lgr1,lgr2)
		#need to clean up
		
	def clearFigures(self):
		
		self.figure_dict={}
	
	def sqw(self,wksp_in,qbin):
		"""
		convert to SmodQw assume direct geom requires string of rebin parameters
		sqw(w1,'0,.1,12')
		"""
		
		n,r=lhs('both')
		wksp_out=r[0]
		#ei= (wksp_in.getSampleDetails().getLogData("Ei").value)
		#wksp_in=mtd[wksp_in]
		SofQW2(InputWorkspace=wksp_in,OutputWorkSpace=wksp_out,QAxisBinning=qbin,EMode="Direct")
		
		CloneWorkspace(InputWorkspace=wksp_in,OutputWorkspace='tmp')
		CreateSingleValuedWorkspace(OutputWorkspace='scale',DataValue='0',ErrorValue='0')
		Multiply(LHSWorkspace='tmp',RHSWorkspace='scale',OutputWorkspace='tmp')
		CreateSingleValuedWorkspace(OutputWorkspace='scale2',DataValue='1',ErrorValue='0')
		Plus(LHSWorkspace='tmp',RHSWorkspace='scale2',OutputWorkspace='tmp')
		
		SofQW2(InputWorkspace='tmp',OutputWorkspace='tmp',QAxisBinning=qbin,EMode='Direct')
		SetUncertaintiesToZero(InputWorkSpace='tmp',OutputWorkSpace='tmp')
		Divide(LHSWorkspace=wksp_out,RHSWorkspace='tmp',OutputWorkspace=wksp_out)
		DeleteWorkspace('tmp')
		DeleteWorkspace('scale')
		DeleteWorkspace('scale2')
		return mtd[wksp_out]
		
	def sqwfast(self,wksp_in,qbin):
		"""
		convert to SmodQw assume direct geom requires string of rebin parameters
		sqw(w1,'0,.1,12')
		"""
		try:
			n,r=lhs('both')
			wksp_out=r[0]
			#ei= (wksp_in.getSampleDetails().getLogData("Ei").value)
			SofQW(wksp_in,wksp_out,QAxisBinning=qbin,EMode="Direct",EFixed=str(ei))
			return mtd[wksp_out]
		except:
			print 'no output workpsace defined'


	def cut(self,wksp,intMin,intMax,minX,delX,maxX,**kwargs):
		'''
		Take a cut out of a wksp created by sqw
		assumes |Q| is X and Energy is Y
		keywords
		along=q
		along=e
		'''
		
		if kwargs.has_key('along') and kwargs.get('along')=='q' or kwargs.has_key('along') and kwargs.get('along')=='Q':
		
			#axis2 is |Q|, axis1 is energy transfer
			cut_name='Cut from '+str(minX)+' to '+str(maxX)+' integrating between '+str(intMin)+' and '+str(intMax)+' meV'
			Rebin2D(InputWorkspace=wksp,OutputWorkspace=cut_name,Axis1Binning=str(intMin)+','+str(intMax-intMin)+','+str(intMax),Axis2Binning=str(minX)+','+str(delX)+','+str(maxX))
			ReplaceSpecialValues(InputWorkspace=cut_name,OutputWorkspace=cut_name,NaNValue='0',InfinityValue='0')
			Transpose(InputWorkspace=cut_name,OutputWorkspace=cut_name)

			
			if kwargs.has_key('over'):
				plotover=self.fignum #assumes that the plot over will be the last plotted figure
				print plotover
				lgr=self.figure_dict.get(self.fignum)[1]#gets the handle of the last fig from the dict
				
				self.fignum=self.fignum+1
				graph=plotSpectrum(cut_name,0,1,1)
				plot=[cut_name,graph]
				self.figure_dict.setdefault(self.fignum,plot)
				mergePlots(graph,lgr)
				
			else:
				self.fignum=self.fignum+1
				graph=plotSpectrum(cut_name,0,1,1)
				plot=[cut_name,graph]
				self.figure_dict.setdefault(self.fignum,plot)
				
				active_layer = graph.activeLayer()
			
		if kwargs.has_key('along') and kwargs.get('along')=='e' or kwargs.has_key('along') and kwargs.get('along')=='E':
		
			cut_name='Cut from '+str(minX)+' to '+str(maxX)+' integrating between '+str(intMin)+' and '+str(intMax)+' A^-1'
			#axis2 is |Q|, axis1 is energy transfer
			ReplaceSpecialValues(InputWorkspace=wksp,OutputWorkspace='tmp',NaNValue='0',InfinityValue='0')
			Rebin2D(InputWorkspace='tmp',OutputWorkspace=cut_name,Axis1Binning=str(minX)+','+str(delX)+','+str(maxX),Axis2Binning=str(intMin)+','+str(intMax-intMin)+','+str(intMax))
			DeleteWorkspace('tmp')
			
			if kwargs.has_key('over'):
				plotover=self.fignum #assumes that the plot over will be the last plotted figure
				print plotover
				lgr=self.figure_dict.get(self.fignum)[1]#gets the handle of the last fig from the dict
				
				self.fignum=self.fignum+1
				graph=plotSpectrum(cut_name,0,1,1)
				plot=[cut_name,graph]
				self.figure_dict.setdefault(self.fignum,plot)
				mergePlots(graph,lgr)
				
			else:
				self.fignum=self.fignum+1
				graph=plotSpectrum(cut_name,0,1,1)
				plot=[cut_name,graph]
				self.figure_dict.setdefault(self.fignum,plot)
				
				active_layer = graph.activeLayer()
				

#def transpose(wksp_in):
#	"""
#	transpose workspace
#	"""
#	n,r=lhs('both')
#	wksp_out=r[0]
#	Transpose(InputWorkspace=wksp_in,OutputWorkspace=wksp_out)
#	return mtd[wksp_out]

def createqtiTable(*args):
#create a qti table of length arg1 with name arg0
	if len(args)==0:
		out=qti.app.newTable()
	if len(args)==2:
		out=qti.app.newTable(args[0],args[1],3)
		out.setColumnRole(3, 5)
	return out 