import os
import sys
import threading
import time
sys.path.append("/home/wzz/Mantid/Code/debug/bin")

from MantidFramework import mtd
mtd.initialise()
from mantidsimple import *

#************************************************
bankid = 1
maxiteration = 10
seeds = range(0, 20)

max_processes = 20
slurm_queue_name = None
#************************************************

# slurm_queue_name = VALID_QUEUE_NAME  This option is to run on a cluster with a slurm set up 

import sys
myscriptname = sys.argv[0]

class MyThread(threading.Thread):

    def setCommand(self, command=""):
        self.command = command

    def run(self):
        # print "STARTING PROCESS: %s" % (self.command)
        os.system(self.command)


joblist = []
active_list = []
script = "LeBailRandomWalk_Parallel.py"
# script = "testscript.py"
for seed in seeds:
    joblist.append( MyThread() )
    outputfile = "output_seed%d.txt" % (seed)
    cmd = 'python %s %d %d %d > %s' % (script, bankid, seed, maxiteration, outputfile)

    if slurm_queue_name is not None:
        # This section is for running on a cluster 
        console_file = "./" + str(seed) + "_output.txt"
        cmd = "srun -p " + slurm_queue_name + ' --cpus-per-task=3 -J ' + myscriptname + ' -o ' + console_file + ' ' + cmd
        print cmd

    joblist[-1].setCommand(cmd)

all_done = False
while not all_done:
    if len(joblist) > 0 and len(active_list) < max_processes:
        nextthread = joblist.pop(0)
        active_list.append(nextthread)
        nextthread.start()
    time.sleep(0.1)

    for thread in active_list:
        if not thread.isAlive():
            active_list.remove(thread)
            print "Remove thread.  Active List Size = %d " % (len(active_list))

    if len(joblist) == 0 and len(active_list) == 0:
        all_done = True

