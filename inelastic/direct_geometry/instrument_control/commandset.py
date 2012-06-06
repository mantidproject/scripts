import time

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
		inp=args[0]
		uamps=float(inp[0])
		print 'counting for uamps=',uamps
		time.sleep(uamps)

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

def end():
	print 'end run'

def waitfor(*args):
	print 'waitfor'
	
