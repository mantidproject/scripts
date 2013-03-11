from genie_init import *
#from mari_init import *
import time
begin()
change(title='test ')
waitfor(seconds=60)
updatestore()
abort()

set_ei(100,150)
aa=get_uamps()
print aa
begin()
waitfor(frames=2000)
abort()

cur=0
countto=1000
print 'waiting for ',countto,' frames'
begin()
while cur <=countto:
	cur=get_frames()
	print 'Counted ', cur, 'of ',countto,' frames'
	#time.sleep(.5)
abort()
waitfor(seconds=10)
time.sleep(10)

