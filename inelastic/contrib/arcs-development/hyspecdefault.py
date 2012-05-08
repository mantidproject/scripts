#this is a dictionary of all possible things which have to be parsed and dealt with
#that belong to the instrument

#ALL KEYWORDS FOR THE DEFAULT DICTIONARY SHOULD BE lowercase.

import math

def instrumentparameters():

	globaldict = dict()


	#instrument items
	globaldict['instrument']='HYSPEC'
	globaldict['filterbadpulses'] = False

	#calibration and mask items
	globaldict['vanruns'] = None
	globaldict['units']   = 'TOF'
	globaldict['vanmin']  = None
	globaldict['vanmax']  = None
	globaldict['processedfilename'] = 'vanadium.nx5'
	globaldict['maskfilename'] = None
	globaldict['mask'] = [] #the mask parameters are stored as a list of dictionaries. (check this)
	globaldict['normalizedcalibration'] = False


	#data items
	globaldict['ipts'] = None  #optional, need to make work for ISIS too.
	globaldict['runs'] = []  #returns a list of ints
	globaldict['efixed']= None #number double
	globaldict['t0'] = None #number double
	globaldict['calce'] = False
	globaldict['ei-mon1-spec'] = 1
	globaldict['ei-mon2-spec'] = 2
	globaldict['emin']=None #number double
	globaldict['emax']=None #number double
	globaldict['ebin']=None #number double
	globaldict['qstep']=None #number double
	globaldict['kiokf']=True #can only use boolean here, not 1 or 0.
	globaldict['tibg']=False
	globaldict['tibgstart']=None
	globaldict['tibgstop']= None
	globaldict['grouping']=None
	globaldict['powderanglestep'] = 0.5
	globaldict['goniometermotor']=None  #script takes care of default situations.
	globaldict['goniometermotoroffset']=None
	globaldict['goniometermotoraxis']=["0,1,0"]
	globaldict['goniometermotordirection']=[1]
	globaldict['save']=[]
	globaldict['friendlyname']='HYSPEC'
	globaldict['friendlynamelogs'] = 'run_number' #will get the run_number from the logs, first run number only.

	return globaldict


def t0fromei(ei):
    return 1.0*(25.0 + 85.0 / (1+math.pow((ei/27.0),4.0)))

