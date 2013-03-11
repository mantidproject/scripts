import sys
sys.path.append("C:\\LabVIEW Modules\\dae\\genie_python")
from liveScriptUI import Ui_MainWindow
from PyQt4 import QtCore, uic,QtGui

import time as time
import numpy

sys.path.append("C:\\LabVIEW Modules\\dae\\genie_python")

from genie_init import *
import commandset

class MainWindow(QtGui.QMainWindow):

	def __init__(self, parent=None):
		QtGui.QMainWindow.__init__(self,parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		QtCore.QObject.connect(self.ui.run, QtCore.SIGNAL("clicked()"), self.start )
		QtCore.QObject.connect(self.ui.update, QtCore.SIGNAL("clicked()"), self.update )
		QtCore.QObject.connect(self.ui.deleteLine, QtCore.SIGNAL("clicked()"), self.deleteLine )
		QtCore.QObject.connect(self.ui.editLine, QtCore.SIGNAL("clicked()"), self.editLine )
		QtCore.QObject.connect(self.ui.load, QtCore.SIGNAL("triggered()"), self.load )
		QtCore.QObject.connect(self.ui.save, QtCore.SIGNAL("triggered()"), self.save )
		QtCore.QObject.connect(self.ui.saveAs, QtCore.SIGNAL("triggered()"), self.saveAs )
		QtCore.QObject.connect(self.ui.stop, QtCore.SIGNAL("clicked()"), self.stop )
		QtCore.QObject.connect(self.ui.updateInfo, QtCore.SIGNAL("clicked()"), self.updateInfo )
		
		self.ui.commands.addItem('stop')
		self.numLines=0
		self.currentline=0
		self.uampsum=0
		self.integrateduamps=0
		listOcommands=commandset.list()
		iter=0
		print get_uamps()
		for command in listOcommands:
			tmp=self.ui.commands.insertItem(iter,command)
 			line= self.ui.commands.item(iter)
 			line.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEditable|QtCore.Qt.ItemIsDragEnabled|QtCore.Qt.ItemIsDropEnabled|QtCore.Qt.ItemIsUserCheckable|QtCore.Qt.ItemIsEnabled)
			iter=iter+1
			
			
	def start(self):
		self.update()
		self.mainThread = GenericThread(self.startScript)
		self.mainThread.start()
		#self.mainThread = Thread(target=self.run_)
		#self.mainThread.start()

	def stop(self):
		print 'stopping script, killing threads'
		line= self.ui.script.item(self.currentline)
		
		#set to green grey =204,204,204
		brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
		brush.setStyle(QtCore.Qt.SolidPattern)
		line.setBackground(brush)
		self.genericThread.terminate()
		self.mainThread.terminate()
		
	def updateInfo(self):
		lines= self.ui.script.count()
		self.uampsum=0
		
		for i in range(lines):
			line= self.ui.script.item(i)
			tmp= str(line.text())
			command, args=self.parseString(tmp)
			if command =='count':
				self.uampsum=self.uampsum+float(args[0])
		
		string='Total Uamps in script = ' +str((self.uampsum))
		self.ui.runInfo.clear()
		self.ui.runInfo.addItem(string)
		
	def load(self):
		
		self.ui.script.clear() 
		fileName = QtGui.QFileDialog.getOpenFileName();
		print 'Load File', fileName
		file = open(fileName, 'r')
		iter=0
		for line in file:
			
			self.ui.script.insertItem(iter,line.rstrip('\n'))
			iter=iter+1
		file.close()
			
	def save(self):
		fileName = QtGui.QFileDialog.getSaveFileName();
		print  'save file', fileName
		file = open(fileName, 'w')
		self.update()
		for row in range(self.numLines):
			line=str( self.ui.script.item(row).text())+'\n'
			file.write(line)
		
		file.close()
		
	def saveAs(self):
		
		fileName = QtGui.QFileDialog.getOpenFileName();
		print  'save  as file ',fileName
		file = open(fileName, 'w')
		self.update()
		for row in range(self.numLines):
			line=str( self.ui.script.item(row).text())+'\n'
			file.write(line)
		
		file.close()
		
	def parseString(self,instring):	
		
		parse=instring.split()
		command=parse[0]
		inputArgs=[]
		for i in range(1,len(parse)):
			inputArgs.append(parse[i])
		return command, inputArgs
		
	def startScript(self):
		
		startnumlines=self.numLines
		self.currentline=0
		self.integrateduamps=0
		status=1
		#for currentline in range(self.numLines):
		while status==1:
			print 'current line',self.currentline
			line= self.ui.script.item(self.currentline)
			tmp= str(line.text())
			#set to green grey =204,204,204
			brush = QtGui.QBrush(QtGui.QColor(0, 255, 0))
			brush.setStyle(QtCore.Qt.SolidPattern)
			line.setBackground(brush)
			
			command, args=self.parseString(tmp)
			
			#set to take any command from commandset.py and pass it to a thread along with the rest of the arguments
			#self.genericThread = GenericThread(self.count,uamps)
			commandToCall = getattr(commandset, command)
			self.genericThread = GenericThread(commandToCall,args)
			self.genericThread.start()
			self.genericThread.wait()
			#while self.genericThread.isRunning() ==True:
			#	print 'waiting'
			#if self.numLines > startnumlines:
			#	currentline=currentline+1
			#	print 'current line condition',currentline
			if tmp=='stop':
				print 'stopping'
				status=0
			if self.currentline == self.numLines-1:
				print 'stopping'
				status=0
			if command == 'count':
				self.integrateduamps=self.integrateduamps+(float(args[0]))
				percentcomplete = (self.integrateduamps/self.uampsum)*100
				string='Completed = ' +str(percentcomplete)+'% of script'
				self.ui.runInfo.addItem(string)
			brush = QtGui.QBrush(QtGui.QColor(204, 204, 204))
			brush.setStyle(QtCore.Qt.SolidPattern)
			line.setBackground(brush)
			self.currentline=self.currentline+1
		print 'finished'
					
			
	def count(self,uamps=0):
		#runs the script
		#count(10)
		print 'counting for uamps=',uamps
		time.sleep(uamps)
		
		
	def update(self):
		#stops the script
		self.numLines= self.ui.script.count()
		print 'updating now ',self.numLines, 'lines in sequence'
		self.updateInfo()
		
	def deleteLine(self):
		#stops the script
		row= self.ui.script.currentRow()
 		self.ui.script.takeItem(row)
	
	def editLine(self):
		#stops the script
		row= self.ui.script.currentRow()
 		line= self.ui.script.item(row)
 		line.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEditable|QtCore.Qt.ItemIsDragEnabled|QtCore.Qt.ItemIsDropEnabled|QtCore.Qt.ItemIsUserCheckable|QtCore.Qt.ItemIsEnabled)

class GenericThread(QtCore.QThread):
	def __init__(self, function, *args, **kwargs):
		QtCore.QThread.__init__(self)
		self.function = function
		self.args = args
		self.kwargs = kwargs

	def __del__(self):
		self.wait()
	
	def run(self):
		try:
			pythoncom.CoInitializeEx(pythoncom.COINIT_MULTITHREADED)
		except pythoncom.com_error:
			pass
		self.function(*self.args,**self.kwargs)
		return

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