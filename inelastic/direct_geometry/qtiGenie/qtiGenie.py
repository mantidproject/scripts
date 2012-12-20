from utils import *
#from mantidsimple import *
import MantidFramework 
from mantid import config
MantidFramework.mtd.initialise()
from DirectEnergyConversion import *
import time as time
import dgreduce
import inspect
import numpy
import PySlice2
import PyChop
# avoid _qti if running outside Mantid.
try:
	import _qti as qti
except ImportError:
  pass 
import os
#########################
#########################
global mpl, wdir, instname,ext,instdae,mon1_spec,mon2_spec,mon3_spec
mpl = 0
if mpl==1:
	print 'Using matplotlib graphics'
if mpl==0:
	print 'Using qtiplot graphics'
wdir=""

instname=""

ext=""

instdae=""

plotmkr=""

plotcolor="red"


save_dir = config.getString('defaultsave.directory')
if len(save_dir) ==0 :
   #set save directory to the directory where qtigenie resides
   save_dir = str(os.path.dirname(inspect.getmodule(dgreduce).__file__))
  
print 'Working directory set to: ',save_dir;
cd(save_dir)

# set default instrument from Mantid configuration
instname = config['default.instrument']
setinst(instname);
print 'Default instrument is set to : ',instname
print 'You can change it by issuing setinst(InstrumentName) command'
print 'where InstrumentName can be MER, MAR, MAP, LET, TSK or XSD'

#set up some 'isis friendly' alias names to dgreduce
iliad_setup=dgreduce.setup
iliad=dgreduce.arb_units
iliad_abs=dgreduce.abs_units
iliad_help=dgreduce.help
iliad_set_calfile = dgreduce.set_cal_file

#######################
#######################

def createqtiTable(*args):
#create a qti table of length arg1 with name arg0
	if len(args)==0:
		out=qti.app.newTable()
	if len(args)==2:
		out=qti.app.newTable(args[0],args[1],3)
		out.setColumnRole(3, 5)
	return out 
					
def fillqtiTable(data):
	return
	
#def arb_units(wb_run,sample_run,ei_guess,rebin,mapfile):
def test(wb_run,sample_run,ei_guess,rebin,mapfile,**kwargs):
	print wb_run,sample_run,ei_guess,rebin,mapfile,kwargs

	for key in kwargs:
		print "another keyword arg: %s: %s" % (key, kwargs[key])
	
	
	if kwargs.has_key('fixei'):
		fix_ei = kwargs.get('fixei')
		print fix_ei
	
	
	
	normalise_method = 'monitor-1'
	background = False
	fix_ei = False
	save_formats = ['.spe']
	#Set parameters for the run
	
	energy_bins = rebin
	background_range=[15000,19000]
	wb_integr_range=[20,100]
	diag_median_rate_limit_hi=3.0
	diag_median_rate_limit_lo=0.1
	bkgd_median_rate_limit=5.0


def listfiles():

	"""
	list function for the working data directory
	"""
	aa= os.listdir(wdir)
	for i in range(0,len(aa)):
		print aa[i]
def print_locals():
	a=locals()
	for i in range(0,len(a)):
		print a[i]

def print_globals():
	a=locals()
	for i in range(0,len(a)):
		print a[i]

##instrument definitions
def setinst(iname):
    """
    setinst('mar')
    setup instrument defaults by reading the instname.txt file
    """
    config['default.instrument'] = iname
    if iname.capitalize()[0:2] =='MER':
        print "Reading params file for merlin from qtiGenie directory"        
        readsetuptxtfile('merlin.txt')

    if iname.capitalize()[0:2]=='MAP':
        print "Reading params file for maps from qtiGenie directory"
        readsetuptxtfile('maps.txt')

    if iname.capitalize()[0:2]=='LET':
        print "Reading params file for LET from qtiGenie directory"
        readsetuptxtfile('let.txt')

    if iname.capitalize()[0:2]=='MAR':
        print "Reading params file for MARI from qtiGenie directory"
        readsetuptxtfile('mari.txt')
    if iname.capitalize()[0:2]=='TSC':
        print "Reading params file for TOSCA from qtiGenie directory"
        readsetuptxtfile('tosca.txt')
    if iname.capitalize()[0:2]=='SXD':
        print "Reading params file for SXD from qtiGenie directory"
        readsetuptxtfile('sxd.txt')

def readsetuptxtfile(fname):
    global wdir, instname,ext,instdae,mon1_spec,mon2_spec,mon3_spec
    
    file_path = str(os.path.dirname(inspect.getmodule(dgreduce).__file__))
    longName = os.path.join(file_path,fname);
    f = open(longName, 'r+')
        
    instring=f.readline()
    instname=instring.split()[1]
    
    instring=f.readline()
    ext=instring.split()[1]
    
    instring=f.readline()
    instdae=instring.split()[1]

    instring=f.readline()
    mon1_spec=int(instring.split()[1])

    instring=f.readline()
    mon2_spec=int(instring.split()[1])

    instring=f.readline()
    mon3_spec=int(instring.split()[1])
    f.close()
    return wdir,instname,ext,instdae,mon1_spec,mon2_spec,mon3_spec
    
def setmon_1_spec(spec):
	"""
	sets the default mon 1 spec to another spectrum
	"""
	global mon1_spec
	mon1_spec=spec
	print 'Monitor one spectrum changed to ',mon1_spec
	return mon1_spec

def showgpath():
    """
    shows the gloabal qtigenie variables
    """
    print '----------------------------------------------------------'
    print 'Global::: ', instname, ': specific variables are:'
    
    print 'Data directory::: ', wdir, ':'
    
    print 'file extension::: ', ext , ':'
    print 'dae path::: ', instdae, ':'
    try:
        mon1_spec
        print 'Monitor 1 spectrum' ,mon1_spec
        print 'Monitor 2 spectrum' ,mon2_spec
        print 'Monitor 3 spectrum' ,mon3_spec
    except NameError:
        print 'No monitors are defined for this function'
    
    getgpath();
    
def getgpath(silent=False):
    """shows global data searh path 
    
    """    
    current_path = config.getDataSearchDirs();
    if not(silent):
        print 'inst_data:::'
        for i in range(0,len(current_path)):
            print '         :::',current_path[i];
            
    return current_path;
def addgpath(path):
    """adds the specified path to data search path
    """
    config.appendDataSearchDir(path)
    
def getspepath(silent=False):
    """Shows the path where qtiGenie and mantid writes its results and temporary working files.
    
    The folder is used for output operations if other path is not specified explicitly by an output operation.    
    if silent opiton is used, the function just returns the path without printing it to std output. 
    """    
    spepath=config.getetString('defaultsave.directory')
    if not(silent):
       print ' Default output data path:',spepath
    
    return spepath;
     
def head(runnumber=0000000,keepWSwithResults=False):
    """Classic head command.
    
    Prints head information defined for the run number specified
    """
    global instname
    
    runnumber=getnumor(runnumber)
    
    runinfo=RawFileInfo(instname+str(runnumber),GetRunParameters=True)

    title =runinfo.getPropertyValue('RunTitle')
    header=runinfo.getPropertyValue('RunHeader')

    paramWSName = 'Raw_RPB'
    temp = mtd[paramWSName] 
    #enddate=temp.getString('r_enddate') 
    #endtime=temp.getString('r_endtime')
  
    print 'RunID\t\t: '+instname+header #+' to '+enddate+endtime
    print 'Title\t\t: '+title     
 
    print 'Protons\t\t:', temp.getDouble('r_gd_prtn_chrg', 0),' uAmps'

    run_length = temp.getInt('r_dur', 0)
    if run_length > 3600 :
        hrs = run_length/3600
        run_length = run_length-hrs*3600
        mins= run_length/60
        sec = run_length - mins*60
        print 'Run duration\t\t:', hrs,' hrs ',mins,' mins ',sec,' sec'
    elif run_length > 60:
        mins= run_length/60
        sec = run_length - mins*60
        print 'Run duration\t\t:',mins,' mins ',sec,' sec' 
    else:
        print 'Run duration\t\t:',run_length,' sec'     
    
    print 'More details availible from Mantid RawFileInfo algorithm'
    if not(keepWSwithResults) :
        mantid.deleteWorkspace(paramWSName)
    
    
#R_dur # r_durunits# r_dur_freq# r_dmp# r_dmp_units# r_dmp_freq#r_freq
#r_gd_prtn_chrg# r_tot_prtn_chrg# r_goodfrm#r_rawfrm# r_dur_wanted#r_dur_secs
#r_mon_sum1# r_mon_sum2#r_mon_sum3# r_enddate#r_endtime#r_prop


def iv(wksp_in):
	"""
	iv(wksp_in)
	initiate the instrument view with workspace wksp_in
	"""
	#bring the mantid instrument view window up
	#need to call with a string??
	wksp=str(wksp_in)
	iv=getInstrumentView(wksp)
	iv.showWindow()

def load(*args):
	"""
	load raw or nxs files
	eg:
	w1=load(1234): will load the file instname1234.raw /nxs
	"""
  # get the lhs of the calling command

	n,r=lhs('both')
	wksp=r[0]
	if (args[0])=='dae': #runnumber == 'dae':
		print 'Access DAE', instdae
		awksp=LoadDAE(instdae,OutputWorkspace=wksp)
		return mtd[wksp]
	else:
	
		runnumber=getnumor(args[0])
		
		
		
		fullname=instname+str(runnumber)
		print fullname
		awksp=Load(fullname,OutputWorkspace=wksp,Cache="Never")
		ConvertToDistribution(wksp)
		runinfo=RawFileInfo(fullname,GetRunParameters=True)
		title=runinfo.getPropertyValue('runtitle')
		return mtd[wksp]

def load_spectra(runnumber,specmin,specmax):
	"""
	load raw or nxs files limited number of spectra
	eg:
	w1=load_spectra(1234,specmin,specmax): will load the file instname1234.raw /nxs
	"""
  # get the lhs of the calling command

	n,r=lhs('both')
	wksp=r[0]
		
	fullname=instname+str(runnumber)
	print fullname
	awksp=Load(fullname,OutputWorkspace=wksp,Cache="Never",SpectrumMin=specmin,SpectrumMax=specmax)
	ConvertToDistribution(wksp)
	runinfo=RawFileInfo(fullname,GetRunParameters=True)
	title=runinfo.getPropertyValue('runtitle')
	return mtd[wksp]		
		

def load_monitors(*args):
	"""
	load monitors from a raw or nxs files
	eg:
	w1=load_monitors(1234): will load the monitors  file instname1234.raw /nxs
	"""
  # load monitors only

	n,r=lhs('both')
	wksp=r[0]

	if len(args)==1: #runnumber == 'dae':

		#print 'Access DAE', instdae

		#mantid.LoadDAE(instdae,OutputWorkspace="data")
		runnumber=getnumor(args[0])
		
		fullname=instname+str(runnumber)
		print fullname

		awksp=Load(fullname,OutputWorkspace=wksp,Cache="Never",LoadMonitors="Separate")
		try: 
			ConvertToDistribution(wksp+'_Monitors')
		except:
			CloneWorkspace(wksp,wksp+'_Monitors')		
				
			
		runinfo=RawFileInfo(fullname,GetRunParameters=True)

		title=runinfo.getPropertyValue('runtitle')
		#clear(wksp)
		DeleteWorkspace(wksp)
		return mtd[wksp+'_Monitors']

	else:
		print 'error'

def getnumor(runnumber):
	"""creates a string runnumber from interger input and pads with zerso to cope with isis 
	   file naminging convention need an additional switch for ts1 and ts2 instruments
	   to cope with the different number of preceeding zeros.
	"""
	if instname=='LET':
		run=str(runnumber)	
		if len(run) == 3:
			runnumber_out='00000'+run
			return runnumber_out	
		if len(run) == 4:
			runnumber_out='0000'+run
			return runnumber_out
		if len(run) == 5:
			runnumber_out='000'+run
			return runnumber_out
		if len(run) == 6:
			runnumber_out='00'+run
			return runnumber_out
		if len(run) == 7:
			runnumber_out='0'+run
			return runnumber_out
		if len(run) == 8:
			runnumber_out=run
			return runnumber_out
	else:
		run=str(runnumber)	
		if len(run) == 3:
			runnumber_out='00'+run
			return runnumber_out	
		if len(run) == 4:
			runnumber_out='0'+run
			return runnumber_out
		if len(run) == 5:
			runnumber_out=run
			return runnumber_out
	
def loadascii(name):
	"""
	loads the x,y,e file as ascii
	w1=loadascii('myfile.txt')
	"""
	#load three col tab delim ascii from current path
	n,r=lhs('both')
	wksp=r[0]
	fullname=os.getcwd()+'/'+name
	print 'Loading ',fullname
	LoadAscii(fullname,OutputWorkspace=wksp,Separator="Tab",Unit="Empty")
	return mtd[wksp]

def ass(wksp):
	"""
	assign a workspace or a raw file as current data for plotting
	ass(w1)
	ass(1234)
	"""
	#define a workspace as current_working_data within the shell
	#simplifies plot etc
	global current_working_data
	if mtd.workspaceExists(str(wksp)):
		current_working_data=wksp
	else:
		tmp=load(wksp)
		current_working_data=tmp
	return current_working_data

def clear(wksp):
	"""
	clear workspace 
	clear(w1)
	"""
	name=wksp.getName()
	mtd.deleteWorkspace(name)
def whos():
	"""
	list all current loaded workspaces
	"""
	names=mtd.getWorkspaceNames()
	print 'Instrument name: ', instname, '\t', 'Data Directory: ',wdir
	print '---------------------------------------------------------------------------'
	print 'WkSp Name', '\t', '\t', '\t', 'Allocated Mem', '\t', 'Title', '\t'
	print '---------------------------------------------------------------------------'
	for i in range(0,len(names)):
		name = names[i]
		tmp=mtd[name]
		mem=tmp.getMemorySize()/1024
		title=tmp.getTitle()
		print name,'\t','\t','\t',mem,'Mb','\t',title
def sync():
	names=mtd.getWorkspaceNames()
	for i in range(0,len(names)):
		name = names[i]
		print name
		#name = mtd[name]
		evalstr='global '+str(name)
		exec(evalstr)
		evalstr1=name+'=mtd['+"'"+str(name)+"'"+']'
		exec(evalstr1)
	return name
	
def default_plotting(inp):
	global mpl	
	print 'set default graphics to matplotlib (1) or qtiplot(0)'
	mpl=inp

def dspacing(wksp_in):
	"""
	convert workspace to dpsacing
	w2=dspacing(w1)
	"""
	#convert input to d spacing and return
	n,r=lhs('both')
	wksp_out=r[0]
	ConvertUnits(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,Target="dSpacing")
	return mtd[wksp_out]

def get_ei(wksp_in,guess):
	"""
	runs getei on a workpsace ei result will be stored in that workspace
	"""
	global mon2_peak
	out=GetEi(wksp_in,mon2_spec,mon3_spec,guess)
	mon2_peak = float(out.getPropertyValue("FirstMonitorPeak"))
	mon2_index = int(out.getPropertyValue("FirstMonitorIndex"))
	ei = (wksp_in.getSampleDetails().getLogData("Ei").value)
	print 'Incident Energy = ', ei, 'meV'
	print 'Peak in monitor',mon2_index, 'at ', mon2_peak ,'usec'
	return ei,mon2_peak
def diag(wb_wksp,run_wksp):
	"""
	runs diag on a white beam and a sample run to list bad spectra
	"""
	#performs the diag test on a WB van and sample run
	#and returns a list of masked detectors
	
	#white beam integrates from 1000-2000 musec
	wb_out=MedianDetectorTest(InputWorkspace=wb_wksp,OutputWorkspace='wb_tmp',StartWorkspaceIndex="1000",EndWorkspaceIndex="2000")
	WB_SpectraList=wb_out.getPropertyValue('BadSpectraNums')
	
	#monorun has threshold set higher and integrate all tof
	run_out=MedianDetectorTest(InputWorkspace=run_wksp,OutputWorkspace='run_tmp',HighThreshold="10")
	Run_SpectraList=run_out.getPropertyValue('BadSpectraNums')
	
	
	bad_dets=maskUnion(Run_SpectraList,WB_SpectraList)
	print 'number of masked spectra = ', len(bad_dets)
	return bad_dets

def maskUnion(a,b):
	a=eval(a)
	b=eval(b)
	aa=list(a)
	bb+list(b)
	return list(set.union(set(aa),set(bb)))

def SaveData(out,run_num):
	"""
	save workspace as spe format
	SaveData(w1,'fname')
	"""
	runnumber=getnumor(run_num)
	fullname=wdir+instname+str(runnumber)+'.spe'
	SaveSPE(out,fullname)




def normalise(*args):
	"""
	normalise 
	w2=normalise(w1,1) <-- normalise to monitor default integration limits of 1000 2000 usec
	w2=normalise(w1,2) <-- normalise to uamp
	w2=normalise(w1,1,1000,2000) <-- normalise to monitor with sepcifies integration limits
	"""
	
	n,r=lhs('both')
	wksp_out=r[0]

	if len(args)==2: #must have two imputs at least input wksp and method
		wksp_in =  args[0]
		method  =  args[1]
		if method == 1:
			print 'Normalise to monitor 1'
			#use defaults for the rest
			mon_spec=mon1_spec
			time_min=1000
			time_max=2000
				
			if wksp_in.isDistribution():
				print 'input wksp is distribution'				
				ConvertFromDistribution(wksp_in)
				NormaliseToMonitor(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,MonitorSpectrum=mon_spec,IntegrationRangeMin=time_min,IntegrationRangeMax=time_max,IncludePartialBins="1")
				#put all data back to distriubution				
				ConvertToDistribution(wksp_out)
				ConvertToDistribution(wksp_in)
			else:
				NormaliseToMonitor(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,MonitorSpectrum=mon_spec,IntegrationRangeMin=time_min,IntegrationRangeMax=time_max,IncludePartialBins="1")

		if method ==2:
			print 'Normalise to current'
			NormaliseByCurrent(wksp_in,wksp_out)
		

	if len(args)==4: #must have two imputs at least input wksp and method
		wksp_in = args[0]
		method  =  args[1]
		#assume normalise to monitor
		if method == 1:
			print 'Normalise to monitor ', args[1],'is spec ',mon1_spec,'between ',args[2],' and ',args[3],' usec'
			#use inputs for the rest
			mon_spec=mon1_spec
			time_min=str(args[2])
			time_max=str(args[3])
			if wksp_in.isDistribution():
				print 'input wksp is distribution'				
				ConvertFromDistribution(wksp_in)
				NormaliseToMonitor(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,MonitorSpectrum=mon_spec,IntegrationRangeMin=time_min,IntegrationRangeMax=time_max,IncludePartialBins="1")
				#put all data back to distriubution				
				ConvertToDistribution(wksp_out)
				ConvertToDistribution(wksp_in)
			else:
				NormaliseToMonitor(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,MonitorSpectrum=mon_spec,IntegrationRangeMin=time_min,IntegrationRangeMax=time_max,IncludePartialBins="1")
		
		
	return mtd[wksp_out]	


def integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax,*args,**kwargs):
	"""
	reads in multiple runs and calculates integral of a region of interest
	as there is no simple method to create a mantid workspace the output is 
	cast as a data_1d object
	integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax)
	or
	integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax,xaxis_start_pos,xstep)
	integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax,xaxis_start_pos,xstep,normalise=mon)
	integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax,xaxis_start_pos,xstep,normalise=mon,mon_range=[1000,2000])
	integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax,xaxis_start_pos,xstep,normalise=uamp)
	"""
	if kwargs.has_key('normalise'):
		norm_method = kwargs.get('normalise')
		if norm_method=='uamp':
			print 'Normalise to uamps'
			norm=2
		elif norm_method=='mon':
			print 'Normalise to monitor 1'
			norm=1
	
	else:
		print 'default uamphr normalisation use keyword normalise to change'
		norm =2

	if kwargs.has_key('normalise') and kwargs.get('normalise')=='mon' and kwargs.has_key('mon_range'):
		mon_int_range=kwargs.get('mon_range')

	jj=1
	outdat=createqtiTable('integral',runstop+1-runstart)
	
	for i in range(runstart,runstop+1):	
		#Extra line added by RAE - check if we are using uamps norm and run has zero uamphr count time
		runnumber=getnumor(i)
		runinfo=RawFileInfo(instname+str(runnumber),GetRunParameters=True)
		temp = mtd.getTableWorkspace('Raw_RPB')
		uamps=temp.getDouble('r_gd_prtn_chrg', 0)
		if uamps==0:
			print 'Integral from run ', i, 'is NaN because zero beam current'
			outdat.setCell(1,jj,jj)
			outdat.setCell(2,jj,0)
			outdat.setCell(3,jj,0)
		else:	
			#tmp=load(i)
			#RAE added:
			tmp=load_spectra(i,specmin,specmax)
			if kwargs.has_key('normalise') and kwargs.get('normalise')=='mon' and kwargs.has_key('mon_range'):		
				out=normalise(tmp,norm,mon_int_range[0],mon_int_range[1])
			else:
				out=normalise(tmp,norm)
			#tmp=integrate(out,tmin,tmax,specmin,specmax)
			tmp=integrate(out,tmin,tmax,1,specmax-specmin+1)
			tmp=sumspec(tmp)
			if len(args)==2:
				outdat.setCell(1,jj,args[0]+args[1]*(jj-1))
			else:
				outdat.setCell(1,jj,jj)
			outdat.setCell(2,jj,tmp.readY(0)[0])
			outdat.setCell(3,jj,tmp.readE(0)[0])
			print 'Integral from run ', i, '=' ,tmp.readY(0)[0],'+/-',tmp.readE(0)[0]
		jj=jj+1
	qti.app.plot(outdat,(1,2,3),2)
	return outdat

def integrate_maps_monitors_over_runs(runstart,runstop,tmin,tmax,mypath):
	"""
	reads in multiple runs and calculates integral of a region of interest
	as there is no simple method to create a mantid workspace the output is 
	cast as a data_1d object
	integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax)
	or
	integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax,xaxis_start_pos,xstep)
	integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax,xaxis_start_pos,xstep,normalise=mon)
	integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax,xaxis_start_pos,xstep,normalise=mon,mon_range=[1000,2000])
	integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax,xaxis_start_pos,xstep,normalise=uamp)
	"""

	norm =2
	jj=1
	outdat=createqtiTable('integral',runstop+1-runstart)
	
	
	for i in range(runstart,runstop+1):	
		#Extra line added by RAE - check if we are using uamps norm and run has zero uamphr count time
		runnumber=getnumor(i)
		runinfo=RawFileInfo(instname+str(runnumber),GetRunParameters=True)
		temp = mtd.getTableWorkspace('Raw_RPB')
		uamps=temp.getDouble('r_gd_prtn_chrg', 0)
		if uamps==0:
			print 'Integral from run ', i, 'is NaN because zero beam current'
			outdat.setCell(1,jj,jj)
			outdat.setCell(2,jj,0)
			outdat.setCell(3,jj,0)
		else:	
			t1=os.path.isfile(mypath+'map'+runnumber+'.raw')
			t2=os.path.isfile(mypath+'map'+runnumber+'.RAW')
			t3=os.path.isfile(mypath+'MAP'+runnumber+'.raw')
			t4=os.path.isfile(mypath+'MAP'+runnumber+'.RAW')
			#check file size
			if t1:
				sz=os.path.getsize(mypath+'map'+runnumber+'.raw')
			elif t2:
				sz=os.path.getsize(mypath+'map'+runnumber+'.RAW')
			elif t3:
				sz=os.path.getsize(mypath+'MAP'+runnumber+'.raw')
			elif t4:
				sz=os.path.getsize(mypath+'MAP'+runnumber+'.RAW')
			
			#Three possibilities, either we are using one to one, or mid-tubes mapping, or a problem. If
			#the file size is not too big, load the whole thing, as it is much easier to get the monitors then
			if sz<50000000:
				tmp=load_monitors(i)
				mflag=True
				badflag=False
			else:
				mflag=False
				try:
					tmp=load_spectra(i,41473,41473)
					specmin=41473
					specmax=41473
					igood=i
					badflag=False
				except:
					try:
						tmp=load_spectra(i,577,577)
						specmin=573
						specmax=573
						badflag=False
					except:
						print('Unexpected item in the bagging area')
						badflag=True
			if badflag:
				print 'Integral from run ', i, 'is NaN because problem getting mon spectra'
				outdat.setCell(1,jj,jj)
				outdat.setCell(2,jj,0)
				outdat.setCell(3,jj,0)
			else:	
				if mflag:
					CloneWorkspace('tmp_Monitors','tmp')			
					#out=normalise(tmp,norm)
					out=tmp/uamps
					tmp=integrate(out,tmin,tmax,1,1)
				else:
					out=normalise(tmp,norm)
					tmp=integrate(out,tmin,tmax,1,specmax-specmin+1)
				tmp=sumspec(tmp)
			
				outdat.setCell(1,jj,jj)
				outdat.setCell(2,jj,tmp.readY(0)[0])
				outdat.setCell(3,jj,tmp.readE(0)[0])
				print 'Integral from run ', i, '=' ,tmp.readY(0)[0],'+/-',tmp.readE(0)[0]
		jj=jj+1
	qti.app.plot(outdat,(1,2,3),2)
	return outdat

def Log(wksp_in):
	"""
	log10 of input workspace
	"""
	n,r=lhs('both')
	wksp_out=r[0]
	Logarithm(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,Natural="0")
	return mtd[wksp_out]
def Ln(wksp_in):
	"""
	loge of input
	"""
	n,r=lhs('both')
	wksp_out=r[0]
	Logarithm(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,Natural="1")
	return mtd[wksp_out]		
def etrans(*args):
	"""
	convert units to energy transfer
	no renormalisation of TOF scale
	convert input to energy transfer and return
	uses the ei stored in the wksp from get_ei as default if ei is specified then the vaule is FIXED to that
	"""	
	n,r=lhs('both')
	wksp_out=r[0]
	
	if len(args)==1:
		wksp_in=args[0]
		ei = float(wksp_in.getSampleDetails().getLogData("Ei").value())
		print 'Converting to energy transfer ei = ', ei,'meV'
		
		ConvertUnits(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,Target="DeltaE",EMode="Direct",EFixed=ei)
	
	if len(args)==2:
		wksp_in=args[0]
		ei=args[1]
		
		print 'Converting to energy transfer ei = ', ei, 'meV'
		
		ConvertUnits(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,Target="DeltaE",EMode="Direct",EFixed=ei)

	
	return mtd[wksp_out]

def avrg_spectra(ws_name,index_min=0,index_max=sys.float_info.max,calc_sigma_avrg=False,include_monitors=False):
    """Get averaged workspace spectrum out of the workpsace spectra in the specified range without using map file.
        
    Usage:
    >>(avrg,stats)=avrg_spectra(ws)
    >>(avrg,stats)=avrg_spectra(ws,index_min,index_max,False,True)
    >>(avrg,stats)=avrg_spectra(ws,index_min,index_max,True,True)    
    >>(avrg,stats)=avrg_spectra(ws,index_min,index_max,calc_sigma_avrg=False,include_monitors=True)    
    Where:
    Input Arguments:
        ws_name          -- the name or the handler of the input workspace to calculae the average spectra. 
        index_min        -- minimal workspace index to sum from default 0 
        index_max        -- maximal workspace index -- maximal number of spectra in the workspace to sum to
                            Detault is max wsorkspace index.
        calc_sigma_avrg  -- calculate the weighted averaged sum. By default this value is False
        include_monitors -- if monitors are in the range of averaging, include them in the sum. 
                            By default, the monitors are excluded.
    Outputs:
       avrg      -- the averaged spectra, the 
       stats     -- 3-element array, which contains number of spectra, contributed into the average,
                    numberOfMaskedSpectra in the input workspace (these spectra were dropped from the sum) and
                    numberOfZeros -- number of zero elements in the input workspace
    """
    # get pointer to the workspace    
    if (type(ws_name) == str):
        ws = mtd[ws_name]
    else:
        ws = ws_name
        
    # check the ws indexes are within limits or defaults
    max_ind_possible = ws.getNumberHistograms()-1
    if (index_max>max_ind_possible):
        index_max = max_ind_possible
    if (index_min < 0):
        index_min = 0
    if (index_max<index_min):
        raise ValueError('Min ws index > Max WS index')        
    
    SumSpectra(InputWorkspace=ws,OutputWorkspace='sumWS',StartWorkspaceIndex=str(index_min),EndWorkspaceIndex=(index_max),
               WeightedSum=calc_sigma_avrg,IncludeMonitors=include_monitors)
    
    
    pOutWS = mtd['sumWS'];
    
    nSpectra       = pOutWS.getRun().getLogData("NumAllSpectra").value
    nMaskedSpectra = pOutWS.getRun().getLogData("NumMaskSpectra").value
    nUsedSpectra   = nSpectra
    nZeroSpectra   = pOutWS.getRun().getLogData("NumZeroSpectra").value
    if (calc_sigma_avrg):
        if nZeroSpectra>0:
           print "->avrg_spectra:: ",nZeroSpectra," spectra out of: ",nUsedSpectra," have have been droped out due to no counts in it"          

        nUsedSpectra-= nZeroSpectra
        if(nUsedSpectra <=0) :
            mtd.deleteWorkspace('sumWS')        
            raise Exception(" no valid spectra found in the workspace")

    spectra = pOutWS.readY(0);  
            
    if len(spectra)==1:
        rez= spectra[0]/nUsedSpectra
    else:
        rez = [None]*len(spectra)           
        for i in range(0,len(spectra)):
            rez[i]=spectra[i]/nUsedSpectra
             
    mtd.deleteWorkspace('sumWS')
    stats =[nUsedSpectra,nMaskedSpectra,nZeroSpectra]
    return (rez,stats);
    
def sumspec(*args):
	"""
	#sums spectra onto a single 1d matrix	
	w2=sumpsec(w1) sum all spec in w1 
	w2=sumpsec(w1,10,100) sum spec 10-->100 in w1
	"""
	n,r=lhs('both')
	wksp_out=r[0]

	if len(args)==1:
		wksp_in=args[0]	
		SumSpectra(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,IncludeMonitors="1")
	if len(args)==3:
		wksp_in=args[0]
		index_1=args[1]-1 #convert to workspace index
		index_2=args[2]-1
		SumSpectra(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,StartWorkspaceIndex=index_1,EndWorkspaceIndex=index_2,IncludeMonitors="1")


	return mtd[wksp_out]

def integrate(*args):
	"""
	#integrates spectra range within limits of x scale 2d	
	w1=integrate(w2) integrate all w2
	w1=integrate(w2,10,1000) integrate w2 within limits of 10 and 1000
	"""
	n,r=lhs('both')
	if len(r)==0:
		wksp_out='tmp_integral'
	if len(r)==1:
		wksp_out=r[0]
	
	#get around bin width issue
	wksp_in=args[0]	
	
	if len(args)==1:
		wksp_in=args[0]	
		Integration(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,IncludePartialBins="1")
	if len(args)==3:
		wksp_in=args[0]
		xrange_lo=args[1]
		xrange_hi=args[2]
		Integration(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,RangeLower=xrange_lo,RangeUpper=xrange_hi,IncludePartialBins="1")

	if len(args)==5:
		wksp_in=args[0]
		xrange_lo=args[1]
		xrange_hi=args[2]
		specrange_lo=args[3]-1#convert to workspace index
		specrange_hi=args[4]-1
		Integration(InputWorkspace=wksp_in,OutputWorkspace=wksp_out,RangeLower=xrange_lo,RangeUpper=xrange_hi,StartWorkspaceIndex=specrange_lo,EndWorkspaceIndex=specrange_hi,IncludePartialBins="1")
	if wksp_out =='tmp_integral':
		print 'Integral = ' ,mtd[wksp_out].readY(0)[0],'+/-',mtd[wksp_out].readE(0)[0]
		
	return mtd[wksp_out]
def transpose(wksp_in):
	"""
	transpose workspace
	"""
	n,r=lhs('both')
	wksp_out=r[0]
	Transpose(InputWorkspace=wksp_in,OutputWorkspace=wksp_out)
	return mtd[wksp_out]

def pwksp(wksp,spec):
	"""
	Plot spectrum from workpsace
	pwksp(w1,10), plot spec 10 from w1
	pwksp(w1,[10,1,100]) plot 10 to 100 from w1
	"""  
	if is_int(spec):
		spec=spec-1	
			
	if is_list(spec):
		spec2matrixsub(spec)
	plotSpectrum(wksp,spec,True)



def is_int(inp):
	#logical type check for integer
	try:
		return int(inp)==inp
	except:
		return False
def is_list(inp):
	#logical type check for list
	try:
		return list(inp)==inp
	except:
		return False		
def spec2matrixsub(inp):
	#convert spectrum number to matric subscipt i.e. -1
	if is_list(inp):
		tmp=inp
		for i in range(len(tmp)):
			tmp[i]=tmp[i]-1
		return tmp
	
                

def getspec(wksp_in,spec):
	"""
	extract specrtum as new workspace
	w2=getspec(w1,10) get spec 10 as w2
	"""
	n,r=lhs('both')
	wksp_out=r[0]
	
	spec=spec-1
	X=list(wksp_in.readX(spec))
	Y=list(wksp_in.readY(spec))
	E=list(wksp_in.readE(spec))
	CreateWorkspace(wksp_out, X, Y, E)
	return mtd[wksp_out]

def smooth(wksp_in,fac):
	"""
	adjacent average smooth of 2d data
	"""
	wksp_out=wksp_in
	SmoothData(wksp_in,wksp_out,NPoints=str(fac))


def pcolor(*args):
	"""
	2d plot using the mantid plot plotting rather that matplotlib
	"""       

	wksp_in=args[0]	
	dat=importMatrixWorkspace(str(wksp_in))
	g2d=dat.plotGraph2D()
		#set z scale 
	if len(args)==2:
		ll=g2d.activeLayer()
		ll.setAxisScale(Layer.Right,0,args[1])
	if len(args)==3:
		ll=g2d.activeLayer()
		ll.setAxisScale(Layer.Right,args[1],args[2])

## MAPS specific scripts ##
def wccr_ang_maps(runnumber):
	"""
	ang=wccr_ang_maps(runnumber)

	If a run title ends with the string wccr=<some angle>, extract that angle
	"""

	n,r=lhs('both')
	ang_out=r[0]	

	runnumber=getnumor(runnumber)	
	runinfo=RawFileInfo(instname+str(runnumber),GetRunParameters=True)
	title=runinfo.getPropertyValue('runtitle')
	#Check that the data do actually come from MAPS here...

	ss=len(title)
	pos=title.find("wccr")
	#Need to catch the possibility here that wccr is not printed in the title

	tt=title[pos+5:ss]
	ang_str=tt.rstrip()
	print "wccr angle is "+ang_str+" degrees"
	ang_out=float(ang_str)
	return ang_out

#########

def dspace_maps(runnumber):
	"""
	wout=dspace_maps(runnumber)

	Equivalent to the old dspace_maps in mgenie. Generates an object and plot for a white
	beam alignment run using mid-tubes mapping, showing intensity on a colour scale with 
	d-spacing and 2theta as the axes
	"""
	n,r=lhs('both')
	wksp_out=r[0]
		
	setinst('maps')
	w1a_=load_spectra(runnumber,1,127)#central strip
	w1b_=load_spectra(runnumber,440,574)
	ConjoinWorkspaces("w1a_","w1b_")
	w2_=dspacing(w1a_)#conjoined wksp is named after the first one
	ConvertSpectrumAxis("w2_","w3_","theta")
	Rebin("w3_","w3a_",[0,0.05,10])
	#w3b_=transpose("w3a_")
	#Rebin("w3b_","w3c_",[0,0.26,60])
	#wout=transpose("w3c_")
	#pcolor("w3c_")
	#Make a copy of wout using CropWorkspace with no additional inputs.
	#CropWorkspace("wout",wksp_out)
	CropWorkspace("w3a_",wksp_out)
	#Get rid of intermediate gubbins:
	DeleteWorkspace("w1a_")
	DeleteWorkspace("w2_")
	DeleteWorkspace("w3_")
	DeleteWorkspace("w3a_")
	#DeleteWorkspace("w3b_")
	#DeleteWorkspace("w3c_")
	#DeleteWorkspace("wout")	
	return mtd[wksp_out]

	
########
def dintegrate(win,d_lo,d_hi):
	"""
	wout=dintegrate(win,d_lo,d_hi)

	Take as input a Workspace2D object created by dspace_maps, and take a cut from it 
	over specified range of d-spacing, d_lo to d_hi. Also plots the output
	"""
	n,r=lhs('both')
	wksp_out=r[0]

	dd=d_hi-d_lo
	#w1=rebin(win,[d_lo,dd,d_hi])
	Rebin(win,"w1",[d_lo,dd,d_hi])
	Transpose("w1","w2")
	Rebin("w2","wout",[0,0.26,60])
	#Make a copy assigned to name of specified output
	CropWorkspace("wout",wksp_out)
	#plotSpectrum("wout",0)
	#Delete intermediate stuff:
	DeleteWorkspace("wout")
	DeleteWorkspace("w1")
	DeleteWorkspace("w2")
	return mtd[wksp_out]

############
def get_maps_spec(win,ang):
	"""
	spectrum_no=get_maps_spec(win,ang)


	Specify a 2theta angle for the workspace2d object win, and the index of the spectrum closest to that
	angle will be returned
	"""
	w2_=transpose(win)
	mylen=len(w2_.readX(0))
	
	for i in range(1,mylen):
		state=(w2_.readX(0)[i])>ang
		if state:
			spec_hi=i
			spec_lo=i-1
			#NB - these are actually workspace indices
			break

	#now find out which is closest:
	d_hi=w2_.readX(0)[spec_hi] - ang
	d_lo=ang - w2_.readX(0)[spec_lo] 
	if d_lo>d_hi:
		specout=spec_hi+1
	elif d_hi>d_lo:
		specout=spec_lo+1
	else:
		specout=spec_lo+1
	#note final statement is to catch the case where angle specified is exactly between two detectors
	angout=w2_.readX(0)[specout-1]
	print "Closest spectrum is no."+str(specout)+" with a scattering angle of "+str(angout)+" degrees"
	DeleteWorkspace("w2_")	
	return specout

##############

def dscan_maps_analysis(run_start,run_end,d_lo,d_hi,specno):
	"""
	wout=dscan_maps_analysis(run_start,run_end,d_lo,d_hi,specno)
	
	Similar to integrate_over_runs, but you specify a single spectrum (determined using get_maps_spec)
	and a lower and upper limit for d-spacing, rather than ToF
	"""
	
	#Initialise qtiTable output:
	jj=1
	outdat=createqtiTable('integral',run_end+1-run_start)	
	
	#Now determine the range of wccr angles from the raw file header info (if a scan command on MAPS was used).
	for i in range(run_start,run_end+1):
		wccr=wccr_ang_maps(i)
		wtmp1_=dspace_maps(i)	
		wtmp2_=sumspec(wtmp1_,specno-1,specno-1)
		ConvertToHistogram("wtmp2_","wtmp3_")
		tmp=integrate("wtmp3_",d_lo,d_hi)
		outdat.setCell(1,jj,wccr)
		outdat.setCell(2,jj,tmp.readY(0)[0])
		outdat.setCell(3,jj,tmp.readE(0)[0])
		print 'Integral from run ', i, '=' ,tmp.readY(0)[0],'+/-',tmp.readE(0)[0]
		jj=jj+1
	qti.app.plot(outdat,(1,2,3),2)
	DeleteWorkspace("wtmp1_")
	DeleteWorkspace("wtmp2_")
	DeleteWorkspace("wtmp3_")
	DeleteWorkspace("tmp")
	
	return outdat

def export_masks(ws,fileName='',returnMasks=False):
    """Exports masks applied to Mantid workspace into the old fashioned ascii msk file with masked spectra numbers
    
      The file is Libisis/Mantid old ISIS format compartible and can be read by libisis or Manid LoadMasks algorithm
	 
	  If optional parameter fileName is present, the masks are saved in the file with this name
	  Otherwise, the file with the name equal to the workspace name and the extension .msk is used.
	  
	  If returnMasks is set to True, the function does not write to file but returns masks array instead
    """
   # get pointer to the workspace    
    if (type(ws) == str):
        pws = mtd[ws]
    else:
        pws = ws
 	
 
    ws_name=pws.getName()       
    nhist = pws.getNumberHistograms()
 
    masks = []
    for i in range(nhist):
        # set provisional spectra ID
        ms = i+1
        try: 
            sp = pws.getSpectrum(i)
            # got real spectra ID, which would correspond real spectra num to spectra ID map
            ms = sp.getSpectrumNo();
        except Exception: 
            print " Can not get spectra No: ",i
            masks.append(ms) 
            continue        
        
        try:
            det = pws.getDetector(i)
        except Exception:
            masks.append(ms)        
            continue
        if det.isMasked():
            masks.append(ms)
 
      

    nMasks = len(masks);
    if nMasks == 0:
        print 'workspace ',ws_name,' have no masked spectra'
        return masks
    print 'workspace ',ws_name,' have ',nMasks,' masked spectra'
    
    filename=''
    if len(fileName)==0 :
        filename=ws_name+'.msk'
    else:
        filename = fileName
        
    if returnMasks :
        return masks
    else:
        writeISISmasks(filename,masks,8)
        
        
def flushOutString(f,OutString,BlockSize,BlockLimit):
    """Internal function for writeISISmasks procedure
    
    """
    BlockSize+=1;
    if BlockSize >= BlockLimit: 
       if len(OutString)>0:
           f.write(OutString+'\n');
       OutString = ''
       BlockSize = 0
    return (f,BlockSize,OutString)

    
def  writeISISmasks(filename,masks,nSpectraInRow=8):
    """ Function writes input array in the form of ISSI mask file array
    
        namely, if one have array 1,2,3,4, 20 30,31,32
        file will have the following ascii stgings:
        1-4 20 30-32
    
    Usage: 
    >>writeISISmasks(fileName,masks)
    where:
    fileName  -- the name of the output file
    masks      -- the array with data
    """
    ext = os.path.splitext(filename)[1]
    if len(ext) == 0 :
        filename=filename+'.msk'

    
    f = open(filename,'w')   
    
    # prepare and write mask data in conventional msk format
    # where adjusted spectra are separated by - sign
    OutString   = ''
    LastSpectraN= ''
    BlockSize = 0;
    iDash = 0;
    im1=masks[0]
    for i in masks:       
        if len(OutString)== 0:
            OutString = str(i)        
            (f,BlockSize,OutString) = flushOutString(f,OutString,BlockSize,nSpectraInRow)
            im1 = i  
            continue
        # if the current spectra is different from the previous one by 1 only, we may want to skip it
        if im1+1 == i :
            LastSpectraN = str(i)
            iDash += 1;
        else :  # it is different and should be dealt separately
            if iDash > 0 :
                OutString = OutString+'-'+LastSpectraN
                iDash = 0
                LastSpectraN=''
                # write the string if it is finished
                (f,BlockSize,OutString) = flushOutString(f,OutString,BlockSize,nSpectraInRow)

      
            if len(OutString) == 0:
                OutString = str(i)
            else:
                OutString = OutString + ' ' + str(i)
            # write the string if it is finished
            (f,BlockSize,OutString) = flushOutString(f,OutString,BlockSize,nSpectraInRow)             
        #endif
      
        # current spectra is the previous now
        im1 = i  
    # end masks loop
    if iDash > 0 :
        OutString = OutString+'-'+LastSpectraN   
    (f,OutString,BlockSize)=flushOutString(f,OutString,BlockSize,0)                   
    f.close();
            
def help(*args):
	if len(args)==0:
		print '!-------------------------------------------------------------------!'    
		print '!                  Mantid Built in Fucntions                        !'
		print '!-------------------------------------------------------------------!'            
		mantidHelp()
        #from inspect import *        
		print '!-------------------------------------------------------------------!'    
		print '!-------------------------------------------------------------------!'    
		print '!                  qtiGenie functions                               !'
		print '!-------------------------------------------------------------------!'    
		print '\t''trim(dat,t1,t2) '
		print '\t''listfiles() '
		print '\t''setinst() '
		print '\t''head(runnumber) '
		print '\t''iv(wksp_in) '
		print '\t''load(*args) '
		print '\t''load_monitors(*args) '
		print '\t''getnumor(runnumber) '
		print '\t''loadascii(name) '
		print '\t''ass(wksp) '	
		print '\t''clear(wksp) '	
		print '\t''whos() '
		print '\t''default_plotting(inp) '
		print '\t''dspacing(wksp_in) '
		print '\t''get_ei(wksp_in,guess) '
		print '\t''normalise(*args) '
		print '\t''rebin(wksp_in,params) '
		print '\t''integrate_over_runs(runstart,runstop,tmin,tmax,specmin,specmax) '	
		print '\t''Log(wksp_in) '	
		print '\t''Ln(wksp_in) '
		print '\t''etrans(*args) '
		print '\t''sumspec(*args) '
		print '\t''integrate(*args) '
		print '\t''transpose(wksp_in) '
		print '\t''pwksp(wksp,spec) '
		print '\t''changecolour(*args) '
		print '\t''changemarker(*args) '
		print '\t''p(spec) '
		print '\t''pe(spec) '
		print '\t''getspec(spec) '
		print '\t''get2d() '
		print '\t''p2d(*args) '
		print '\t''psurf(*args) '
		print '\t''plus(a,b) '
		print '\t''minus(a,b) '
		print '\t''mult(a,b) '
		print '\t''div(a,b) '         
		print '\t''help() '
		print '--------------------------------------------------------------------'	
		print 'qtiGenie classes'
		print '--------------------------------------------------------------------'
		print 'class Data_1D() '
		print '--------------------------------------------------------------------'
		print '\t''__init__(Data_1D) '  
		print '\t''plot(wksp) '
		print '\t''plotwe(wksp) '
		print '\t''Add_1d(wksp,factor) '
		print '\t''Minus_1D(wksp,factor) '
		print '\t''Multiply_1d(wksp,factor) '
		print '\t''Divide_1d(wksp,factor) '
		print '\t''integrate(wksp,t1,t2) '
		print '\t''xscale(dat,t1,t2) '
		print '--------------------------------------------------------------------'	
		print 'class Data_2D()'
		print '--------------------------------------------------------------------'
		print '\t''__init__(Data_2D) '
		print '\t''sumspec(wksp,*args) '
		print '\t''sum2d(wsk,s1,s2) '
	else:
		execstr='print '+str(args[0])+'.__doc__'
		exec(execstr)
