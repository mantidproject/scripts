from utils import *
from mantidsimple import *
import MantidFramework 
MantidFramework.mtd.initialise()
from mantidplot import *
from DirectEnergyConversion import *
import time as time
import dgreduce
import inspect
import numpy
import qti



class data2D:
	global active_layer, graph
	def __init__(self,wksp):
		self.data=0
		self.init_params()
		self.wksp_name=wksp
	def init_params(self):
		self.linecol=2
		

	def rebinproj(self,qbin,**kwargs):
		wksp_in=self.wksp_name
		wksp_new_name=wksp_in.getName()+'_SQW'
		if kwargs.has_key('fast'):
			print 'Using fast rebin'
			tmp=sqwfast(wksp_in,qbin)
			RenameWorkspace(tmp,wksp_new_name)
			self.data=mtd[wksp_new_name]
		else:
			print 'Using intersecting area rebin'
			tmp=sqw(wksp_in,qbin)
			RenameWorkspace(tmp,wksp_new_name)
			self.data=mtd[wksp_new_name]
		
	def colour(self,col):
		active_layer.setCurveLineColor(1,col)
	
	def SimpleCutE(self,Qmin,Qmax):
		cutsimple(self.data,Qmin,Qmax,along='e')
	
	def SimpleCutQ(self,Emin,Emax):
		cutsimple(self.data,Emin,Emax,along='q')
	
	def CutAlongE(self,Qmin,Qmax,Emin,delE,Emax,**kwargs):
		if kwargs.has_key('over'):
			cut(self.data,Qmin,Qmax,Emin,delE,Emax,along='e',over=True)
		else:
			cut(self.data,Qmin,Qmax,Emin,delE,Emax,along='e')
	
	def CutAlongQ(self,Emin,Emax,Qmin,delQ,Qmax,**kwargs):
		if kwargs.has_key('over'):
			cut(self.data,Emin,Emax,Qmin,delQ,Qmax,along='q',over=True)
		else:
			cut(self.data,Emin,Emax,Qmin,delQ,Qmax,along='q')
			
	
	def XaxisLims(self,xscale_min,xscale_max):
		setscale(xscale_min,xscale_max,'x')
	
	def YaxisLims(self,yscale_min,yscale_max):
		setscale(yscale_min,yscale_max,'y')
	
	def display(self,*args):
		data=self.data	
		#dat=self.data.importMatrixWorkspace().plotGraph2D()
		out=importMatrixWorkspace(data.getName())
		g2D=out.plotGraph2D()
		#set z scale 
		if len(args)==1:
			ll=g2D.activeLayer()
			ll.setAxisScale(Layer.Right,0,args[0])
		if len(args)==2:
			ll=g2D.activeLayer()
			ll.setAxisScale(Layer.Right,args[0],args[1])
			
def sqw(wksp_in,qbin):
	"""
	convert to SmodQw assume direct geom requires string of rebin parameters
	sqw(w1,'0,.1,12')
	"""
	try:
		n,r=lhs('both')
		wksp_out=r[0]
		ei= (wksp_in.getSampleDetails().getLogData("Ei").value)

		SofQW2(wksp_in,wksp_out,QAxisBinning=qbin,EMode="Direct",EFixed=str(ei))
		Transpose(InputWorkspace=wksp_out,OutputWorkspace=wksp_out)
		return mtd[wksp_out]
	except:
		print 'no output workpsace defined'

def sqwfast(wksp_in,qbin):
	"""
	convert to SmodQw assume direct geom requires string of rebin parameters
	sqw(w1,'0,.1,12')
	"""
	try:
		n,r=lhs('both')
		wksp_out=r[0]
		ei= (wksp_in.getSampleDetails().getLogData("Ei").value)

		SofQW(wksp_in,wksp_out,QAxisBinning=qbin,EMode="Direct",EFixed=str(ei))
		Transpose(InputWorkspace=wksp_out,OutputWorkspace=wksp_out)
		return mtd[wksp_out]
	except:
		print 'no output workpsace defined'
	
def cutsimple(wksp,min,max,**kwargs):
	'''
	Take a cut out of a wksp created by sqw
	assumes |Q| is X and Energy is Y
	keywords
	along=q
	along=e
	'''
	global active_layer, graph
	
	if kwargs.has_key('along') and kwargs.get('along')=='e' or kwargs.has_key('along') and kwargs.get('along')=='E':
	
		#for e cut x axis of 2d dataset must be in q i.e. label must be 1/Angstrom as a string
		if wksp.getAxis(0).getUnit().label()=='1/Angstrom':
			print 'dimension correct'
		if wksp.getAxis(0).getUnit().label()=='meV':
			print 'dimesnion incorrect transposing workspace'
			Transpose(InputWorkspace=wksp,OutputWorkspace=wksp)
		cut_name='Cut along E integrating between '+str(min)+' and '+str(max)+' A^-1'
		Rebin(InputWorkspace=wksp,OutputWorkspace=cut_name,Params=str(min)+','+str(max-min)+','+str(max),PreserveEvents='0')
		Transpose(InputWorkspace=cut_name,OutputWorkspace=cut_name)	
		plotSpectrum(cut_name,0,1)
		
		
	
	if kwargs.has_key('along') and kwargs.get('along')=='q' or kwargs.has_key('along') and kwargs.get('along')=='Q':
		
		#for q cut x axis of 2d dataset must be in meV i.e. label must be meV as a string
		if wksp.getAxis(0).getUnit().label()=='meV':
			print 'dimension correct'
		if wksp.getAxis(0).getUnit().label()=='1/Angstrom':
			print 'dimesnion incorrect transposing workspace'
			Transpose(InputWorkspace=wksp,OutputWorkspace=wksp)
		
		cut_name='Cut along |Q| integrating between '+str(min)+' and '+str(max)+' meV'
		Rebin(InputWorkspace=wksp,OutputWorkspace=cut_name,Params=str(min)+','+str(max-min)+','+str(max),PreserveEvents='0')
		Transpose(InputWorkspace=cut_name,OutputWorkspace=cut_name)	
		plotSpectrum(cut_name,0,1)

def cut(wksp,intMin,intMax,minX,delX,maxX,**kwargs):
	'''
	Take a cut out of a wksp created by sqw
	assumes |Q| is X and Energy is Y
	keywords
	along=q
	along=e
	'''
	global active_layer, graph
	
	if kwargs.has_key('along') and kwargs.get('along')=='e' or kwargs.has_key('along') and kwargs.get('along')=='E':
	
		#for e cut x axis of 2d dataset must be in q i.e. label must be 1/Angstrom as a string
		if wksp.getAxis(0).getUnit().label()=='1/Angstrom':
			print 'dimension correct'
		if wksp.getAxis(0).getUnit().label()=='meV':
			print 'dimesnion incorrect transposing workspace'
			Transpose(InputWorkspace=wksp,OutputWorkspace=wksp)
		cut_name='Cut along E integrating between '+str(intMin)+' and '+str(intMax)+' A^-1'
		Rebin(InputWorkspace=wksp,OutputWorkspace=cut_name,Params=str(intMin)+','+str(intMax-intMin)+','+str(intMax),PreserveEvents='0')
		Transpose(InputWorkspace=cut_name,OutputWorkspace=cut_name)	
		Rebin(cut_name,cut_name,str(minX)+','+str(delX)+','+str(maxX),PreserveEvents='0')
		
		if kwargs.has_key('over'):
			cut_name=mtd[cut_name]
			table=creatQtiTableFromWorkSpace(cut_name,0)
			active_layer.addCurves(table,(1,2,3),2)
		else:
			graph=plotSpectrum(cut_name,0,1)
			active_layer = graph.activeLayer()
		
	
	if kwargs.has_key('along') and kwargs.get('along')=='q' or kwargs.has_key('along') and kwargs.get('along')=='Q':
		
		#for q cut x axis of 2d dataset must be in meV i.e. label must be meV as a string
		if wksp.getAxis(0).getUnit().label()=='meV':
			print 'dimension correct'
		if wksp.getAxis(0).getUnit().label()=='1/Angstrom':
			print 'dimesnion incorrect transposing workspace'
			Transpose(InputWorkspace=wksp,OutputWorkspace=wksp)
		
		cut_name='Cut along |Q| integrating between '+str(intMin)+' and '+str(intMax)+' meV'
		Rebin(InputWorkspace=wksp,OutputWorkspace=cut_name,Params=str(intMin)+','+str(intMax-intMin)+','+str(intMax),PreserveEvents='0')
		Transpose(InputWorkspace=cut_name,OutputWorkspace=cut_name)	
		Rebin(cut_name,cut_name,str(minX)+','+str(delX)+','+str(maxX),PreserveEvents='0')
		
		if kwargs.has_key('over'):
			cut_name=mtd[cut_name]
			table=creatQtiTableFromWorkSpace(cut_name,0)
			active_layer.addCurves(table,(1,2,3),2)
		else:
			graph=plotSpectrum(cut_name,0,1)
			active_layer = graph.activeLayer()
		
		
def plotFromQtiTable(tb_in,labelx,labely,title):
	global active_layer, graph
	graph=qti.app.plot(tb_in,(1,2,3),2)
	active_layer = graph.activeLayer()
	#l.setLineColor(1, 2)
	active_layer.setTitle(title)
	active_layer.setAxisTitle(0, labelx)
	active_layer.setAxisTitle(2, labely)
	#l.setScale(0,0,50)# set y scale
	#l.setScale(2,0,50)#set xscale

def plotOverFromQtiTable(tb_in,labelx,labely,title):
	global active_layer, graph
	active_layer.addCurves(tb_in,(1,2,3),2)
	#l.setLineColor(1, 2)
	#l.setScale(0,0,50)# set y scale
	#l.setScale(2,0,50)#set xscale

#def makeCurrent():
#	global active_layer, graph
#	gr=qti.app.currentGraph()

def creatQtiTableFromWorkSpace(wksp_in,index):
	out=qti.app.newTable(wksp_in.getName(),wksp_in.getNumberBins()+1,3)
	out.setColumnRole(3, 5)
	
	xx=wksp_in.readX(index)
	yy=wksp_in.readY(index)
	ee=wksp_in.readE(index)
	
	jj=1
	for i in range(wksp_in.getNumberBins()-1):
		out.setCell(1,jj,xx[jj])
		out.setCell(2,jj,yy[jj])
		out.setCell(3,jj,ee[jj])
		jj=jj+1
	return out

def setscale(min,max,dir):
	if dir == 'x':
		active_layer.setScale(2,min,max)#set xscale
	if dir == 'y':
		active_layer.setScale(0,min,max)#set xscale

def transpose(wksp_in):
	"""
	transpose workspace
	"""
	n,r=lhs('both')
	wksp_out=r[0]
	Transpose(InputWorkspace=wksp_in,OutputWorkspace=wksp_out)
	return mtd[wksp_out]

def createqtiTable(*args):
#create a qti table of length arg1 with name arg0
	if len(args)==0:
		out=qti.app.newTable()
	if len(args)==2:
		out=qti.app.newTable(args[0],args[1],3)
		out.setColumnRole(3, 5)
	return out 