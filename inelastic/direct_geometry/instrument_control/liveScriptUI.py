# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/jon/QtSDK/liveScript.ui'
#
# Created: Wed Jun  6 23:54:07 2012
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(707, 571)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.script = QtGui.QListWidget(self.centralwidget)
        self.script.setGeometry(QtCore.QRect(20, 40, 311, 401))
        self.script.setAcceptDrops(True)
        self.script.setEditTriggers(QtGui.QAbstractItemView.DoubleClicked)
        self.script.setDragEnabled(True)
        self.script.setDragDropMode(QtGui.QAbstractItemView.DragDrop)
        self.script.setDefaultDropAction(QtCore.Qt.CopyAction)
        self.script.setSelectionBehavior(QtGui.QAbstractItemView.SelectItems)
        self.script.setObjectName(_fromUtf8("script"))
        self.commands = QtGui.QListWidget(self.centralwidget)
        self.commands.setGeometry(QtCore.QRect(360, 40, 331, 192))
        self.commands.setAcceptDrops(True)
        self.commands.setEditTriggers(QtGui.QAbstractItemView.DoubleClicked)
        self.commands.setDragEnabled(True)
        self.commands.setDragDropMode(QtGui.QAbstractItemView.DragDrop)
        self.commands.setDefaultDropAction(QtCore.Qt.CopyAction)
        self.commands.setObjectName(_fromUtf8("commands"))
        self.run = QtGui.QPushButton(self.centralwidget)
        self.run.setGeometry(QtCore.QRect(360, 240, 114, 32))
        self.run.setObjectName(_fromUtf8("run"))
        self.update = QtGui.QPushButton(self.centralwidget)
        self.update.setGeometry(QtCore.QRect(220, 450, 114, 32))
        self.update.setObjectName(_fromUtf8("update"))
        self.deleteLine = QtGui.QPushButton(self.centralwidget)
        self.deleteLine.setGeometry(QtCore.QRect(10, 480, 114, 32))
        self.deleteLine.setObjectName(_fromUtf8("deleteLine"))
        self.editLine = QtGui.QPushButton(self.centralwidget)
        self.editLine.setGeometry(QtCore.QRect(10, 450, 114, 32))
        self.editLine.setObjectName(_fromUtf8("editLine"))
        self.stop = QtGui.QPushButton(self.centralwidget)
        self.stop.setGeometry(QtCore.QRect(580, 240, 114, 32))
        self.stop.setObjectName(_fromUtf8("stop"))
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(20, 10, 101, 20))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(362, 10, 111, 20))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.runInfo = QtGui.QListWidget(self.centralwidget)
        self.runInfo.setGeometry(QtCore.QRect(360, 310, 331, 131))
        self.runInfo.setObjectName(_fromUtf8("runInfo"))
        self.updateInfo = QtGui.QPushButton(self.centralwidget)
        self.updateInfo.setGeometry(QtCore.QRect(580, 450, 114, 32))
        self.updateInfo.setObjectName(_fromUtf8("updateInfo"))
        self.label_3 = QtGui.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(360, 280, 151, 20))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 707, 22))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionLoad = QtGui.QAction(MainWindow)
        self.actionLoad.setObjectName(_fromUtf8("actionLoad"))
        self.actionSave = QtGui.QAction(MainWindow)
        self.actionSave.setObjectName(_fromUtf8("actionSave"))
        self.load = QtGui.QAction(MainWindow)
        self.load.setObjectName(_fromUtf8("load"))
        self.save = QtGui.QAction(MainWindow)
        self.save.setObjectName(_fromUtf8("save"))
        self.saveAs = QtGui.QAction(MainWindow)
        self.saveAs.setObjectName(_fromUtf8("saveAs"))
        self.menuFile.addAction(self.load)
        self.menuFile.addAction(self.save)
        self.menuFile.addAction(self.saveAs)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.run.setText(QtGui.QApplication.translate("MainWindow", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.update.setText(QtGui.QApplication.translate("MainWindow", "update", None, QtGui.QApplication.UnicodeUTF8))
        self.deleteLine.setText(QtGui.QApplication.translate("MainWindow", "Delete Line", None, QtGui.QApplication.UnicodeUTF8))
        self.editLine.setText(QtGui.QApplication.translate("MainWindow", "Edit Line", None, QtGui.QApplication.UnicodeUTF8))
        self.stop.setText(QtGui.QApplication.translate("MainWindow", "Stop", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "Sequence", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "Commands", None, QtGui.QApplication.UnicodeUTF8))
        self.updateInfo.setText(QtGui.QApplication.translate("MainWindow", "update", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "Run Infomation", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setTitle(QtGui.QApplication.translate("MainWindow", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.actionLoad.setText(QtGui.QApplication.translate("MainWindow", "load", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave.setText(QtGui.QApplication.translate("MainWindow", "save", None, QtGui.QApplication.UnicodeUTF8))
        self.load.setText(QtGui.QApplication.translate("MainWindow", "Load", None, QtGui.QApplication.UnicodeUTF8))
        self.save.setText(QtGui.QApplication.translate("MainWindow", "Save", None, QtGui.QApplication.UnicodeUTF8))
        self.saveAs.setText(QtGui.QApplication.translate("MainWindow", "Save As", None, QtGui.QApplication.UnicodeUTF8))

