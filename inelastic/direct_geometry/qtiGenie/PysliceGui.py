
import sys
from PySliceUI2 import Ui_MainWindow
from PyQt4 import QtCore, uic,QtGui
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
from PySlice2 import *
from qtiGenie import *
from autoEi import *
#class MainWindow(QtGui.QMainWindow, Ui_MainWindow):

class MainWindow(QtGui.QMainWindow):

	def __init__(self, parent=None):
		QtGui.QMainWindow.__init__(self,parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		QtCore.QObject.connect(self.ui.plot, QtCore.SIGNAL("clicked()"), self.plot )
		QtCore.QObject.connect(self.ui.Reduce, QtCore.SIGNAL("clicked()"), self.Reduce )
		QtCore.QObject.connect(self.ui.plotOver, QtCore.SIGNAL("clicked()"), self.plotOver )
		QtCore.QObject.connect(self.ui.calcProj, QtCore.SIGNAL("clicked()"), self.calcProj )
		QtCore.QObject.connect(self.ui.display, QtCore.SIGNAL("clicked()"), self.display )
		QtCore.QObject.connect(self.ui.refresh, QtCore.SIGNAL("clicked()"), self.refreshlist )
		QtCore.QObject.connect(self.ui.cutMin, QtCore.SIGNAL("returnPressed()"), self.cutMin)
		QtCore.QObject.connect(self.ui.cutMax, QtCore.SIGNAL("returnPressed()"), self.cutMax)
		QtCore.QObject.connect(self.ui.delta, QtCore.SIGNAL("returnPressed()"), self.delta)
		QtCore.QObject.connect(self.ui.intMin, QtCore.SIGNAL("returnPressed()"), self.intMin)
		QtCore.QObject.connect(self.ui.intMax, QtCore.SIGNAL("returnPressed()"), self.intMax)
		QtCore.QObject.connect(self.ui.WkspIn, QtCore.SIGNAL("activated()"), self.setwksp)
		QtCore.QObject.connect(self.ui.actionDensityOfstates, QtCore.SIGNAL("triggered()"), self.pdos )
		QtCore.QObject.connect(self.ui.actionBoseFactor, QtCore.SIGNAL("triggered()"), self.bose )
		QtCore.QObject.connect(self.ui.removeWksp, QtCore.SIGNAL("clicked()"), self.removeWksp )
		QtCore.QObject.connect(self.ui.removePlot, QtCore.SIGNAL("clicked()"), self.removePlot )
		
		self.smooth=0
		self.cutMin=0
		self.cutMax=0
		self.delta=0
		self.intMin=0
		self.intMax=0
		self.scalemin=''
		self.scalemax=''
		self.selectedWksp=''
		self.selectedFig=''
		self.data_dict={}
		self.masterFigureDict={}
		self.fignum=0
		self.currentWkspNames=[]
		
		self.ui.axis.insertItem(1,'|Q|')
		self.ui.axis.insertItem(2,'E')
		
		self.currentWkspNames=mtd.getObjectNames()
		iter=0
		for item in self.currentWkspNames: 
			self.ui.WkspIn.insertItem(iter,item)
			iter=iter+1
		
	 
 	def Reduce(self):
 		WB=self.ui.DetVanNum.text()
 		Run=self.ui.RunNum.text()
 		inst='mar'
		iliad_setup(inst)
		ext='.raw'
		mapfile='mari_res2012'
		cal_file='MAR'+WB+'.raw'
		ei,rebin_params=autoEi(str(Run))
		w1=iliad(WB,Run,ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current')
		RenameWorkspace(InputWorkspace='w1',OutputWorkspace='mar'+str(Run))
		
 	def setwksp(self):
 		print 'high'
 	
 	def removeWksp(self):
 		tmp=self.ui.wkspList.selectedItems()[0].text()
 		row= self.ui.wkspList.currentRow()
 		self.ui.wkspList.takeItem(row)
 		Wksp = tmp.split(':')[0]
 		
 		tmpdata=self.data_dict.get(str(Wksp))
 		DeleteWorkspace(tmpdata.data)
 		DeleteWorkspace(tmpdata.data_disp)
 		del self.data_dict[str(Wksp)]
 	
 	def removePlot(self):
 		tmp=self.ui.FigureList.selectedItems()[0].text()
 		row= self.ui.FigureList.currentRow()
 		self.ui.FigureList.takeItem(row)
 		fig = tmp.split(':')[0]
 		del self.masterFigureDict[int(fig)]
 		
 	def bose(self):
 		temp=float(self.ui.Temperature.text())
 		
 		if temp >=0.0:
 			print 'Correct for Bose factor'
 			tmp=self.ui.wkspList.selectedItems()[0].text()
 			row= self.ui.wkspList.currentRow()
 			self.ui.wkspList.takeItem(row)
 			
			self.selectedWksp = tmp.split(':')[0]
			string=tmp+'Bose Factor Applied with T='+" %.2f" % temp+'K'
			self.ui.wkspList.insertItem(row,string)
			data=self.data_dict.get(str(self.selectedWksp))
			data.boseFac(temp)
 		else:
 			print 'Enter temperature in kelvin for run'
 			
 		
 	def pdos(self):
 		temp=float(self.ui.Temperature.text())
 		dbwfac=float(self.ui.DBWFac.text())
 		
 		if temp >=0.0:
 			print 'Calculate Density of States'
 			tmp=self.ui.wkspList.selectedItems()[0].text()
 			row= self.ui.wkspList.currentRow()
 			self.ui.wkspList.takeItem(row)
			self.selectedWksp = tmp.split(':')[0]
			string=tmp+'Density of states calculated with T='+" %.2f" % temp+'K'
			self.ui.wkspList.insertItem(row,string)
			data=self.data_dict.get(str(self.selectedWksp))
			data.gofe(temp,dbwfac)	
		else:
 			print 'Enter temperature in kelvin for run'
 			
 	def refreshlist(self): 
		self.ui.WkspIn.clear() 
 		self.currentWkspNames=mtd.getObjectNames()
		iter=0
		for item in self.currentWkspNames: 
			self.ui.WkspIn.insertItem(iter,item)
			iter=iter+1
			
	def cutMin(self):
		self.cutMin=self.ui.cutMin.text()
		#print self.cutMin
	def cutMax(self):
		self.cutMax=self.ui.cutMax.text()
		#print self.cutMax
	def delta(self):
		self.delta=self.ui.delta.text()
		#print self.delta
	def intMin(self):
		self.intMin=self.ui.intMin.text()
		#print self.intMin
	def intMax(self):
		self.intMax=self.ui.intMax.text()
		#print self.intMax
		
	def plot(self):
		
		#try:
		self.intMax=float(self.ui.intMax.text())
		self.intMin=float(self.ui.intMin.text())
		self.delta=float(self.ui.delta.text())
		self.cutMin=float(self.ui.cutMin.text())
		self.cutMax=float(self.ui.cutMax.text())
		self.scalemin=(self.ui.scalemin.text())
		self.scalemax=(self.ui.scalemax.text())
		
		#print self.cutMin
		#print self.cutMax
		#print self.delta
		#print self.intMin
		#print self.intMax
		#print self.selectedWksp
		try:
			dir=self.ui.axis.currentText()
			tmp=self.ui.wkspList.selectedItems()[0].text()
			self.selectedWksp = tmp.split(':')[0]
			data=self.data_dict.get(str(self.selectedWksp))
		except:
			print 'select workspace to plot from'
		if dir=='E':
			print 'Cut along E'
			self.fignum=self.fignum+1
			data.ECut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,shoelace=True)
			entry=data.figure_dict.get(data.fignum)
			
			fighandle=entry[1]
			if len(self.scalemax) and len(self.scalemin)>0:
				self.scaleY(fighandle)
			
			self.masterFigureDict.setdefault(self.fignum,entry)
			#entry[0]is title entry[1]is handle
			string=str(self.fignum)+':'+entry[0]
			self.ui.FigureList.addItem(string)
			
		if dir =='|Q|':
			print 'Cut along |Q|'
			self.fignum=self.fignum+1
			data.QCut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,shoelace=True)
			entry=data.figure_dict.get(data.fignum)
			fighandle=entry[1]
			if len(self.scalemax) and len(self.scalemin)>0:
				self.scaleY(fighandle)
			self.masterFigureDict.setdefault(self.fignum,entry)
			#entry[0]is title entry[1]is handle
			string=str(self.fignum)+':'+entry[0]
			self.ui.FigureList.addItem(string)
			
		#except:
		 #return
		 
	def plotOver(self):
		
		
		self.intMax=float(self.ui.intMax.text())
		self.intMin=float(self.ui.intMin.text())
		self.delta=float(self.ui.delta.text())
		self.cutMin=float(self.ui.cutMin.text())
		self.cutMax=float(self.ui.cutMax.text())
		self.scalemin=(self.ui.scalemin.text())
		self.scalemax=(self.ui.scalemax.text())
		#print self.cutMin
		#print self.cutMax
		#print self.delta
		#print self.intMin
		#print self.intMax
		#print self.selectedWksp
		dir=self.ui.axis.currentText()
		
		try:
			tmp=self.ui.wkspList.selectedItems()[0].text()
			self.selectedWksp = tmp.split(':')[0]
			data=self.data_dict.get(str(self.selectedWksp))
			tmp=self.ui.FigureList.selectedItems()[0].text()
			self.selectedFig = tmp.split(':')[0]
			#print self.selectedFig
			figHandle=self.masterFigureDict.get(int(self.selectedFig))[1]
			#print self.selectedFig
			#print self.masterFigureDict
			#print figHandle
		except:
			print 'Select both workspace to plot from and figure to plot onto'
		if dir=='E':
			print 'Plot over along E'
			
			data.ECut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,over=True,Handle=figHandle,shoelace=True)
			
			entry=data.figure_dict.get(data.fignum)
			fighandle=entry[1]
			if len(self.scalemax) and len(self.scalemin)>0:
				self.scaleY(fighandle)
			
			del self.masterFigureDict[int(self.selectedFig)]
			self.masterFigureDict.setdefault(int(self.selectedFig),entry)
			
		if dir =='|Q|':
			print 'Plot Over along |Q|'
			data.QCut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,over=True,Handle=figHandle,shoelace=True)
			
			entry=data.figure_dict.get(data.fignum)
			fighandle=entry[1]
			if len(self.scalemax) and len(self.scalemin)>0:
				self.scaleY(fighandle)
			del self.masterFigureDict[int(self.selectedFig)]
			self.masterFigureDict.setdefault(int(self.selectedFig),entry)
			
		
	def scaleY(self,handle):
		min=(self.ui.scalemin.text())
		max=(self.ui.scalemax.text())
		active_layer = handle.activeLayer()
		active_layer.setScale(0,float(min),float(max))
			 
	def display(self):
		smooth= self.ui.smooth.text()
		IMin= self.ui.DispMin.text()
		IMax= self.ui.DispMax.text()
		
		if len(IMax)>0 and len(IMin)==0:
			IMin='0'
		
		try:
			tmp=self.ui.wkspList.selectedItems()[0].text()
			self.selectedWksp = tmp.split(':')[0]
			#print self.selectedWksp
			data=self.data_dict.get(str(self.selectedWksp))
			
			
			if len(smooth)>0 and len(IMin)>0 and len(IMax)>0:
				data.smooth(int(smooth),int(IMin),int(IMax))
				
			if len(smooth)>0 and len(IMin)==0 and len(IMax)==0:
				data.smooth(int(smooth))
				
			if len(smooth)==0 and len(IMin)>0 and len(IMax)>0:	
				data.display(int(IMin),int(IMax))
				
			if len(smooth)==0 and len(IMin)==0 and len(IMax)==0:	
				data.display()
				
		except:
			print 'Select data to display'
				
	def calcProj(self):
	
		wkspName=str(self.ui.WkspIn.currentText())
		
		if self.data_dict.get(str(wkspName)) != None:
			print 'S(qw) data from ',wkspName,'exists, deleting existing instance'
			tmpdata=self.data_dict.get(str(wkspName))
 			DeleteWorkspace(tmpdata.data)
 			DeleteWorkspace(tmpdata.data_disp)
 			del self.data_dict[str(wkspName)]
 			#now deal with the displayed text
 			numLines=self.ui.wkspList.count()
 			for i in range(numLines):
 				tmp= self.ui.wkspList.item(i).text()
				wkspExists=tmp.split(':')[0]
				
				if wkspExists==wkspName:
					self.ui.wkspList.takeItem(i)
					
				
		data=mtd[str(self.ui.WkspIn.currentText())]
		numHist=data.getNumberHistograms()
		minTheta=data.detectorTwoTheta(data.getDetector(0))
		maxTheta=data.detectorTwoTheta(data.getDetector(numHist-1))
		#calculate extents in mod q for the data and generate rebin parameters
		
		ei= data.getRun().getLogData('Ei').value
		ki=sqrt(ei/2.07)
		kf=sqrt((ei-0.0)/2.07)
		Qx=ki*1-cos(minTheta)*kf
		
		Qy=-(sin(minTheta)*cos(0))*ki
		Qz=-(sin(minTheta)*sin(0))*ki
		Qmin=sqrt(Qx**2+Qy**2+Qz**2)
		#print Qx,Qy,Qz
		
		Qx=ki*1-cos(maxTheta)*kf
		Qy=-(sin(maxTheta)*cos(0))*ki
		Qz=-(sin(maxTheta)*sin(0))*ki
		Qmax=sqrt(Qx**2+Qy**2+Qz**2)
		
		deltaQ=(Qmax-Qmin)/numHist
		text=str(self.ui.WkspIn.currentText())
		string1=text+': ei= '+" %.2f" % ei+'meV qmin= '+" %.2f" % Qmin+'qmax = '+" %.2f" % Qmax
		rebin=str(Qmin)+','+str(deltaQ)+','+str(Qmax)
		self.ui.wkspList.addItem(string1)
		#print self.ui.WkspIn.currentText()
		instring=str(self.ui.WkspIn.currentText())
		data=data2D(instring)
		#newVarName=instring+data
		#new_var_name='some_new_variable_name'
		#exec '%s = tmp'%new_var_name
		data.rebinproj(rebin)
		self.data_dict.setdefault(instring,data)
		
def qapp():
	if QtGui.QApplication.instance():
		app = QtGui.QApplication.instance()
	else:
		app = QtGui.QApplication(sys.argv)
	return app
 
app = qapp()
reducer = MainWindow()
reducer.show()
app.exec_()