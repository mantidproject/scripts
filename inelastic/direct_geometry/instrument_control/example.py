from mari_control import *
#begin()
#end()
#abort()
#defineRoi(5.0,5.5,50.0,60.0)
#run.runfile='16654' # this would be the hard coded way of setting the runfile for reduction
#seteiforcontrol(100)
#count(1,2)
#countStats(4,5,0)
#settemp(10)
#set_ei(50,150)
#cset(T_head=100)
#changetitle('Change the run title to this)
#updatestore()
#count(.5)

setWhiteBeam('16993')

ei=50
seteiforcontrol(ei)
defineRoi(0.0,3.0,-5.0,5.0)
title='Test python control MonoVan'+str(ei)+'meV 3Gd'
changetitle(title)
countstats(1,40,50)

defineRoi(0.0,5.0,20.0,40.0)
title='Test python control MonoVan'+str(ei)+'meV 3Gd'
changetitle(title)
countstats(1,60,500)

set_ei(100,300)
ei=100
seteiforcontrol(ei)
defineRoi(3.0,7.0,25,60)
title='Test python control MonoVan'+str(ei)+'meV 6Gd'
changetitle(title)
countstats(1.5,60,500)