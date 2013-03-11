import sys
sys.path.append("C:\\LabVIEW Modules\\dae\\genie_python")

#from genie_api import *
from genie_init import *
import time
from qtiGenie import *



from mari_control import mariControl


run=mariControl()

defineRoi=run.defineRoi

set_ei=run.set_ei

setWhiteBeam=run.setWhiteBeam

reduction_ei=run.reduction_ei
seteiforcontrol=run.setEiForControl

count=run.count

countStats=run.countStats

settemp=run.settemp

changetitle=run.changetitle

updatestore=run.updatestore

seteifixed=run.seteifixed

seteifree=run.seteifree

list=run.list

scanTemp=run.scanTemp
begin=run.begin
cset=run.cset

end =run.end
print get_uamps()
#def list():
#	#returns a list of displayable command names must be a nicer way todo this!
#	return ['count uamps',
#	'settemp temp',
#	'changetitle runtitle',
#	'set_ei ei freq',
#	'begin',
#	'end',
#	'waitfor uamps']
	
