import time
from qtiGenie import *
from PySlice2 import *
import time
from genie_init import *


def list():
	#returns a list of displayable command names must be a nicer way todo this!
	return ['count uamps',
	'settemp temp',
	'changetitle runtitle',
	'set_ei ei freq',
	'begin',
	'end',
	'waitfor uamps']
	
def count(*args):
	#runs the script
	#count(10)
	updateTime=2.0
	inp=args[0]
	uamps=float(inp[0])
	print 'counting for uamps=',uamps
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

def settemp(*args,**kwargs):
	inp=args[0]
	print 'set temp',float(inp[0])

def changetitle(*args):
	inp=args[0]
	print 'title',' '.join(inp)


def set_ei(*args):
	inp=args[0]
	
	print 'setenergy'
	print 'ei',inp[0]
	print 'freq',inp[1]

def begin():
	print 'begin run'
	begin()

def end():
	print 'end run'
	end()
	
def waitfor(*args):
	inp=args[0]
	uamps=float(inp[0])
	print 'waitfor ',uamps,'uamps'
	totalcurrent=uamps
	currentUamps=0
	starttime=time.time()
	updateTime=2.0
	while currentUamps<totalcurrent:
		currentUamps=get_uamps()

		time.sleep(updateTime)
		
		elaspedtime=time.time()-starttime
		print 'Waiting for',uamps, 'uamps, now at = ',currentUamps,' Runtime= ',int(elaspedtime),' secs'
