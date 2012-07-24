import sys
from Ui_MainWindow import Ui_MainWindow #import line for the UI python class
from PyQt4 import QtCore, uic,QtGui#PyQt imports required (may need more)
import sys
import math 
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

class MainWindow(QtGui.QMainWindow):
	needsThetaInputList = ['Momentum transfer (Q Angstroms^-1)', 'd-Spacing (Angstroms)']
	
	needsThetaOutputList = ['Momentum transfer (Q Angstroms^-1)', 'd-Spacing (Angstroms)']
	
	needsFlightPathInputList = ['Time of flight (microseconds)']
	
	needsFlightPathOutputList = ['Time of flight (microseconds)']
	
	
	def thetaEnable (self, enabled):
		self.ui.lineEdit_4.setEnabled(enabled)
		if  enabled == False:
			self.ui.lineEdit_4.clear()
	def flightPathEnable (self, enabled):
		self.ui.lineEdit_3.setEnabled(enabled)
		if  enabled == False:
			self.ui.lineEdit_3.clear()
	def setInstrumentInputs (self):
		print"setInstrumentInputs"
		#disable both
		self.thetaEnable(False)
		self.flightPathEnable(False)
		#get the values of the two unit strings
		inOption=self.ui.inputUnits.currentText()
		outOption=self.ui.outputUnits.currentText()
		#for theta: enable if input or output unit requires it 
		
		if inOption in self.needsThetaInputList:
			self.thetaEnable(True)
		
		if outOption in self.needsThetaOutputList:
			self.thetaEnable(True)
		
		#for flightpath: enable if input or output unit requires it
		
		if inOption in self.needsFlightPathInputList:
			self.flightPathEnable(True)
		
		if outOption in self.needsFlightPathOutputList:
			self.flightPathEnable(True)
			
	def __init__(self, parent=None):
		QtGui.QMainWindow.__init__(self,parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		self.ui.InputVal.setValidator(QtGui.QDoubleValidator())
		self.ui.lineEdit_3.setValidator(QtGui.QDoubleValidator())
		self.ui.lineEdit_4.setValidator(QtGui.QDoubleValidator())
		QtCore.QObject.connect(self.ui.convert, QtCore.SIGNAL("clicked()"), self.convert )
		QtCore.QObject.connect(self.ui.inputUnits, QtCore.SIGNAL("currentIndexChanged(QString)"), self.setInstrumentInputs )
		QtCore.QObject.connect(self.ui.outputUnits, QtCore.SIGNAL("currentIndexChanged(QString)"), self.setInstrumentInputs )
		self.setInstrumentInputs() 
		##defaults
		self.flightpath=-1.0
		self.Theta=-1.0
		self.stage1output=0.0
		self.stage2output=0.0
	def convert(self):
		print 'default flightpath=',self.flightpath
		try:
			inOption=self.ui.inputUnits.currentText()
			outOption=self.ui.outputUnits.currentText()
			if self.ui.lineEdit_3.text() !='':
				self.flightpath = float(self.ui.lineEdit_3.text())
			if self.ui.lineEdit_4.text() !='':
				self.Theta = ((float(self.ui.lineEdit_4.text()) * 3.14159265358979323846264) / 360.0)
		
			print 'input value=',self.ui.InputVal.text()
		
			self.stage1output= self.input2energy(float(self.ui.InputVal.text()),inOption)
			self.stage2output= self.energy2output(self.stage1output,outOption)
		
			self.ui.lineEdit_2.clear()
			self.ui.lineEdit_2.insert(str(self.stage2output))
		except:
			print 'Input numbers please'
			return
	
	def input2energy(self,inputval, inOption):
		e2lam = 81.787
		e2nu = 4.139
		e2v = 0.0000052276
		e2k = 2.717
		e2t = 0.086165
		ed2cm = 0.123975

		pi = 3.14159265358979323846264
		iv2 = inputval ** 2 


		if inOption == 'Wavelength (Angstroms)':
			Energy = e2lam / iv2

		elif inOption == 'Energy  (meV)':
			Energy = inputval
	
	
		elif inOption == 'Nu (THz)':
			Energy = e2nu * inputval

		elif inOption == 'Velocity (m/s)':
			Energy = e2v *iv2

		elif inOption == 'Momentum ( k Angstroms^-1)':
			Energy = e2k*iv2

		elif inOption == 'Temperature (K)':
			Energy = e2t *inputval

		elif inOption == 'Energy (cm^-1)':
			Energy = e2cm * inputval

		elif inOption == 'Momentum transfer (Q Angstroms^-1)':
			
			if self.Theta >=0:
				k = inputval * 0.5 / math.sin(self.Theta)
				Energy = e2k * k * k
			else:
				print 'Theta needs to be defined'
				sys.exit
			print self.Theta
		elif inOption == 'd-Spacing (Angstroms)':
			lam = 2 * inputval * math.sin(self.Theta)
			Energy = e2lam / (lam * lam)

		elif  inOption == 'Time of flight (microseconds)':
			Energy = 1000000 * self.flightpath
			Energy = e2v * Energy *Energy / iv2

		return Energy

	def energy2output(self, Energy, inOption):
		e2lam = 81.787
		e2nu = 4.139
		e2v = 0.0000052276
		e2k = 2.0717
		e2t = 0.086165
		e2cm = 0.123975

		pi = 3.14159265358979323846264
		iv2 = Energy ** 2 
		print(inOption)
		if inOption == 'Wavelength (Angstroms)':
			OutputVal =  (e2lam/ Energy)**0.5

		elif inOption == 'Nu (THz)':

			OutputVal = Energy / e2nu

		elif inOption == 'Velocity (m/s)':
			OutputVal = (Energy / e2v)**0.5

		elif inOption == 'Momentum ( k Angstroms^-1)':
			OutputVal = (Energy / e2k)**0.5

		elif inOption == 'Temperature (K)':
			OutputVal = Energy / e2t

		elif inOption == 'Energy (cm^-1)':
			OutputVal = Energy / e2cm

		elif inOption == 'Momentum transfer (Q Angstroms^-1)':
			k = (Energy / e2k) ** 0.5
			OutputVal = 2 * k * math.sin(self.Theta)

		elif inOption == 'd-Spacing (Angstroms)':
			lam = (e2lam / Energy)**0.5
			OutputVal = lam * 0.5 / math.sin(self.Theta)

		elif inOption == 'Time of flight (microseconds)':
			OutputVal = self.flightpath * 1000 * ((e2v * 1000000 / Energy) ** 0.5)
		
		elif inOption == 'Energy  (meV)':
			OutputVal = Energy
		  
		return OutputVal



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