
import sys
from PysliceUI import Ui_MainWindow #import line for the UI python class
from PyQt4 import QtCore, uic,QtGui#PyQt imports required (may need more)


import time as time
##import and init the mantid framework
import MantidFramework # import and init the mantid framework 
MantidFramework.mtd.initialise()
from mantidplotpy import *
import inspect
import numpy
from mantidplot import *
from mantid import *
from mantid.simpleapi import *



#class MainWindow(QtGui.QMainWindow, Ui_MainWindow):

class MainWindow(QtGui.QMainWindow):

	def __init__(self, parent=None):
		QtGui.QMainWindow.__init__(self,parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		#connects to a button press
		QtCore.QObject.connect(self.ui.plot, QtCore.SIGNAL("clicked()"), self.plot ) 
		
		#connects to a text box when return is pressed
		QtCore.QObject.connect(self.ui.cutMin, QtCore.SIGNAL("returnPressed()"), self.cutMin)
		
		#connects to a menu item
		QtCore.QObject.connect(self.ui.actionDensityOfstates, QtCore.SIGNAL("triggered()"), self.pdos )
		
		
		
		
		#insert itmes
		self.ui.axis.insertItem(1,'|Q|')
		self.ui.axis.insertItem(2,'E')
		#in a loop
		self.currentWkspNames=mtd.getObjectNames() #get some items 
		iter=0
		for item in self.currentWkspNames: 
			self.ui.WkspIn.insertItem(iter,item)
			iter=iter+1
		
 	
 	def removeWksp(self):
 		tmp=self.ui.wkspList.selectedItems()[0].text() #get selected row text from ui
 		row= self.ui.wkspList.currentRow() #get current row integer from ui 
 		self.ui.wkspList.takeItem(row)#delete item from ui
 		Wksp = tmp.split(':')[0]
 		
 		
		string=tmp+'Bose Factor Applied with T='+" %.2f" % temp+'K'
		self.ui.wkspList.insertItem(row,string)#insert a string to a selected row
		
 		else:
 			print 'Enter temperature in kelvin for run'
 			
 		
 	
 			
 
		self.cutMin=self.ui.cutMin.text() #get text from a uitext box
		self.intMax=float(self.ui.intMax.text())#get text from a text box and cast to a float
		#print self.cutMin
	
	
		#calc modQ for a fixed energy on the elastic line could change for inelastic scattering
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
		
		
#this is required to make the UI run within the mantid environment and remain stable
def qapp():
	if QtGui.QApplication.instance():
		app = QtGui.QApplication.instance()
	else:
		app = QtGui.QApplication(sys.argv)
	return app
 
app = qapp()
reducer = MainWindow()#the main ui class in this file is called MainWindow
reducer.show()
app.exec_()