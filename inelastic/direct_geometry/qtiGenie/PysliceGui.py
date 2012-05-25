
import sys
from pysliceUI import Ui_MainWindow
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
#class MainWindow(QtGui.QMainWindow, Ui_MainWindow):
class MainWindow(QtGui.QMainWindow):

	def __init__(self, parent=None):
		QtGui.QMainWindow.__init__(self,parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		QtCore.QObject.connect(self.ui.plot, QtCore.SIGNAL("clicked()"), self.plot )
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
		
		self.smooth=0
		self.cutMin=0
		self.cutMax=0
		self.delta=0
		self.intMin=0
		self.intMax=0
		self.selectedWksp=''
		self.selectedFig=''
		self.data_dict={}
		self.masterFigureDict={}
		self.fignum=0
		self.ui.axis.insertItem(1,'|Q|')
		self.ui.axis.insertItem(2,'E')
		
		tmpwkps=mtd.getObjectNames()
		iter=1
		for item in tmpwkps: 
			self.ui.WkspIn.insertItem(iter,item)
			iter=iter+1
		
	 
 	def setwksp(self):
 		print 'high'
 		
 	def refreshlist(self):
 		tmpwkps=mtd.getObjectNames()
		iter=1
		for item in tmpwkps: 
			self.ui.WkspIn.insertItem(iter,item)
			iter=iter+1
 		
	def cutMin(self):
		self.cutMin=self.ui.cutMin.text()
		print self.cutMin
	def cutMax(self):
		self.cutMax=self.ui.cutMax.text()
		print self.cutMax
	def delta(self):
		self.delta=self.ui.delta.text()
		print self.delta
	def intMin(self):
		self.intMin=self.ui.intMin.text()
		print self.intMin
	def intMax(self):
		self.intMax=self.ui.intMax.text()
		print self.intMax
		
	def plot(self):
		
		#try:
		self.intMax=float(self.ui.intMax.text())
		self.intMin=float(self.ui.intMin.text())
		self.delta=float(self.ui.delta.text())
		self.cutMin=float(self.ui.cutMin.text())
		self.cutMax=float(self.ui.cutMax.text())
		print self.cutMin
		print self.cutMax
		print self.delta
		print self.intMin
		print self.intMax
		print self.selectedWksp
		dir=self.ui.axis.currentText()
		tmp=self.ui.wkspList.selectedItems()[0].text()
		self.selectedWksp = tmp.split(':')[0]
		data=self.data_dict.get(str(self.selectedWksp))
		if dir=='E':
			print 'Cut along E'
			self.fignum=self.fignum+1
			data.ECut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,shoelace=True)
			entry=data.figure_dict.get(data.fignum)
			self.masterFigureDict.setdefault(self.fignum,entry)
			#entry[0]is title entry[1]is handle
			string=str(self.fignum)+':'+entry[0]
			self.ui.FigureList.addItem(string)
			
		if dir =='|Q|':
			print 'Cut along |Q|'
			self.fignum=self.fignum+1
			data.QCut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,shoelace=True)
			entry=data.figure_dict.get(data.fignum)
			self.masterFigureDict.setdefault(self.fignum,entry)
			#entry[0]is title entry[1]is handle
			string=str(self.fignum)+':'+entry[0]
			self.ui.FigureList.addItem(string)
			
		#except:
		 #return
		 
	def plotOver(self):
		
		#try:
		self.intMax=float(self.ui.intMax.text())
		self.intMin=float(self.ui.intMin.text())
		self.delta=float(self.ui.delta.text())
		self.cutMin=float(self.ui.cutMin.text())
		self.cutMax=float(self.ui.cutMax.text())
		print self.cutMin
		print self.cutMax
		print self.delta
		print self.intMin
		print self.intMax
		print self.selectedWksp
		dir=self.ui.axis.currentText()
		
		tmp=self.ui.wkspList.selectedItems()[0].text()
		self.selectedWksp = tmp.split(':')[0]
		data=self.data_dict.get(str(self.selectedWksp))
		
		if dir=='E':
			print 'Plot over along E'
			tmp=self.ui.FigureList.selectedItems()[0].text()
			self.selectedFig = tmp.split(':')[0]
			print self.selectedFig
			figHandle=self.masterFigureDict.get(int(self.selectedFig))[1]
			print self.selectedFig
			print self.masterFigureDict
			print figHandle
			data.ECut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,over=True,Handle=figHandle,shoelace=True)
			entry=data.figure_dict.get(data.fignum)
			print entry
			del self.masterFigureDict[int(self.selectedFig)]
			self.masterFigureDict.setdefault(int(self.selectedFig),entry)
			#entry[0]is title entry[1]is handle
			#string=str(self.fignum)+':'+entry[0]
			#self.ui.FigureList.addItem(string)
		if dir =='|Q|':
			print 'Plot Over along |Q|'
			tmp=self.ui.FigureList.selectedItems()[0].text()
			self.selectedFig = tmp.split(':')[0]
			print self.selectedFig
			figHandle=self.masterFigureDict.get(int(self.selectedFig))[1]
			print self.selectedFig
			print self.masterFigureDict
			print figHandle
			data.QCut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,over=True,Handle=figHandle,shoelace=True)
			entry=data.figure_dict.get(data.fignum)
			del self.masterFigureDict[int(self.selectedFig)]
			self.masterFigureDict.setdefault(int(self.selectedFig),entry)
			#entry[0]is title entry[1]is handle
			#string=str(self.fignum)+':'+entry[0]
			#self.ui.FigureList.addItem(string)
		#except:
		 #return
		 
	def display(self):
		smooth= self.ui.smooth.text()
		try:
			tmp=self.ui.wkspList.selectedItems()[0].text()
			self.selectedWksp = tmp.split(':')[0]
			print self.selectedWksp
			data=self.data_dict.get(str(self.selectedWksp))

			if len(smooth)>0:
				data.smooth(int(smooth))
			else:	
				data.display()
		except:
			print 'Select data to display'
				
	def calcProj(self):
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
		print Qx,Qy,Qz
		
		Qx=ki*1-cos(maxTheta)*kf
		Qy=-(sin(maxTheta)*cos(0))*ki
		Qz=-(sin(maxTheta)*sin(0))*ki
		Qmax=sqrt(Qx**2+Qy**2+Qz**2)
		
		deltaQ=(Qmax-Qmin)/numHist
		text=str(self.ui.WkspIn.currentText())
		string1=text+': ei= '+str(ei)+'meV qmin= '+str(Qmin)+'qmax = '+str(Qmax)
		rebin=str(Qmin)+','+str(deltaQ)+','+str(Qmax)
		self.ui.wkspList.addItem(string1)
		print self.ui.WkspIn.currentText()
		instring=str(self.ui.WkspIn.currentText())
		data=data2D(instring)
		#newVarName=instring+data
		#new_var_name='some_new_variable_name'
		#exec '%s = tmp'%new_var_name
		data.rebinproj(rebin)
		self.data_dict.setdefault(instring,data)
		
app = QtGui.QApplication(sys.argv)
reducer = MainWindow()
reducer.show()
app.exec_()