from utils import *
#from mantidsimple import *
import MantidFramework 
MantidFramework.mtd.initialise()
from DirectEnergyConversion import *
import CommonFunctions as common
import time as time
import numpy


def setup(instname):
	"""
	setup('mar')
	setup redcution parameters from instname_parameter.xml file
	"""
	global reducer, inst_name,van_mass,bleed_switch,rate,pixels
	if instname=='MAR' or instname=='mar':
		print 'setup mari'
		inst_name='MAR'
		reducer = setup_reducer('MARI')
		bleed_switch=False
		rate=0.0
		pixels=0
		van_mass=32.58
	elif instname=='MER' or instname=='mer':
		print 'setup merlin'
		inst_name='MER'
		reducer = setup_reducer('MERLIN')
		bleed_switch=True
		rate=0.01
		pixels=80
		van_mass=32.62
	elif instname=='MAP' or instname=='map':
		print 'setup maps'
		inst_name='MAP'
		reducer = setup_reducer('MAPS')
		bleed_switch=False
		rate=0.0
		pixels=0.0
		van_mass=32.58
	elif instname=='LET' or instname=='let':
		print 'setup let'
		inst_name='LET'
		reducer = setup_reducer('LET')
		bleed_switch=True
		rate=0.01
		pixels=80
		van_mass=32.58
	elif instname=='ARCS' or instname=='arcs':
		print 'setup Arcs'
		inst_name='ARC'
		reducer = setup_reducer('ARCS')
		bleed_switch=False
		rate=0.01
		pixels=80
		van_mass=32.58
	elif instname=='SEQ' or instname=='seq':
		print 'setup Sequoia'
		inst_name='SEQ'
		reducer = setup_reducer('SEQUOIA')
		bleed_switch=False
		rate=0.01
		pixels=80
		van_mass=32.58
	elif instname=='CNCS' or instname=='cncs':
		print 'setup cncs'
		inst_name='SEQ'
		reducer = setup_reducer('CNCS')
		bleed_switch=False
		rate=0.01
		pixels=80
		van_mass=32.58
	elif instname=='HYSPEC' or instname=='hyspec':
		print 'setup hyspec'
		inst_name='SEQ'
		reducer = setup_reducer('HYSPEC')
		bleed_switch=False
		rate=0.01
		pixels=80
		van_mass=32.58
	else:
		print 'Instrument name not defined'
		return

def help(*args):
	print 'available keywords for dgreduce with default values'
	print 'norm_method- normalistion monitor-1,monitor-2,uamph'
	print  "normalise_method = 'monitor-1'"
	print 'background subtract background True/False'
	print	"background = False"
	print 'fixei- fix ei to input value True/False'
	print	"fix_ei = False"
	print 'save_format- output file type'
	print	"save_formats = ['.spe']"
	print 'bkgd_range- integration range of background list or True for auto'
	print	"background_range=[15000,19000]"
	print 'detector_van_range- integratin in E(mev) for detector vanadium data'
	print	"wb_integr_range=[20,100]"
	print 'diag_sigma- diag sigma override'
	print	"diag_sigma=3.0"
	print 'abs_units_van_range- integration for absolute units correction vanadium data'
	print	"abs_units_van_range=[-40,40]"

def arb_units(wb_run,sample_run,ei_guess,rebin,map_file,**kwargs):
	"""
	arb_units(wb_run,sample_run,ei_guess,rebin,mapfile,**kwargs)
	
	dgreduce.arb_units(1000,10001,80,'-10,.1,70','mari_res', additional keywords as required)
	
	dgreduce.arb_units(1000,10001,80,'-10,.1,70','mari_res',fixei=True)
	
	Available keywords
	norm_method =[monitor-1],[monitor-2][Current]
	background  =False , True
	fixei 		=False , True
	save_format	=['.spe'],['.nxspe'],'none'
	detector_van_range			=[20,40] in mev
	
	bkgd_range	=[15000,19000]	:integration range for background tests
	
	second_white	- If provided an additional set of tests is performed on this. (default = None)
	hard_mask  		- A file specifying those spectra that should be masked without testing (default=None)
	tiny        	- Minimum threshold for acceptance (default = 1e-10)
	large        	- Maximum threshold for acceptance (default = 1e10)
	bkgd_range 		- A list of two numbers indicating the background range (default=instrument defaults)
	diag_van_median_rate_limit_lo  	- Lower bound defining outliers as fraction of median value (default = 0.01)
	diag_van_median_rate_limit_hi 	- Upper bound defining outliers as fraction of median value (default = 100.)
	diag_van_median_sigma_lo      	- Fraction of median to consider counting low for the white beam diag (default = 0.1)
	diag_van_median_sigma_hi      	- Fraction of median to consider counting high for the white beam diag (default = 1.5)
	diag_van_sig  - Error criterion as a multiple of error bar i.e. to fail the test, the magnitude of the
			  		difference with respect to the median value must also exceed this number of error bars (default=0.0)
	diag_remove_zero 				- If true then zeroes in the vanadium data will count as failed (default = True)
	diag_samp_samp_median_sigma_lo  - Fraction of median to consider counting low for the white beam diag (default = 0)
	diag_samp_samp_median_sigma_hi  - Fraction of median to consider counting high for the white beam diag (default = 2.0)
	diag_samp_sig  					- Error criterion as a multiple of error bar i.e. to fail the test, the magnitude of the"
			  						  difference with respect to the median value must also exceed this number of error bars (default=3.3)
	variation  		-The number of medians the ratio of the first/second white beam can deviate from
			   		the average by (default=1.1)
	bleed_test 		- If true then the CreatePSDBleedMask algorithm is run
	bleed_maxrate 	- If the bleed test is on then this is the maximum framerate allowed in a tube
	bleed_pixels 	- If the bleed test is on then this is the number of pixels ignored within the
				 	bleed test diagnostic
	print_results - If True then the results are printed to the screen
	
	diag_remove_zero =True, False (default):Diag zero counts in background range
	bleed=True , turn bleed correction on and off on by default for Merlin and LET
	
	sum	=True,False(default) , sum multiple files
	
	det_cal_file= a valid detector block file and path or a raw file. Setting this
				  will use the detector calibraion from the specified file NOT the
				  input raw file
	mask_run = RunNumber to use for diag instead of the input run number
	
	one2one =True, False :Reduction will not use a mapping file
	
	hardmaskPlus=Filename :load a hardmarkfile and apply together with diag mask
	
	hardmaskOnly=Filename :load a hardmask and use as only mask
	
	"""
	global reducer, rm_zero,inst_name,van_mass,bleed_switch,rate,pixels
	print 'DGreduce run for ',inst_name,'run number ',sample_run
	try:
		n,r=lhs('both')
		wksp_out=r[0]
	except:
		if sample_run == 0:
			#deal with the current run being parsed as 0 rather than 00000
			sample_run='00000'
		wksp_out=inst_name+str(sample_run)+'.spe'
		if kwargs.has_key('sum') and kwargs.get('sum')==True:
			wksp_out=inst_name+str(sample_run[0])+'sum'+'.spe'
		
	start_time=time.time()
	
	if sample_run=='00000' and mtd.workspaceExists(inst_name+'00000.raw')==True:
		print 'Deleteing previous instance of temp data'
		DeleteWorkspace(inst_name+'00000.raw')
	
	#repopulate defualts
	if kwargs.has_key('norm_method'):
		reducer.normalise_method = kwargs.get('norm_method')
		print 'Setting normalisation method to ', kwargs.get('norm_method')
	else:
		reducer.normalise_method = 'monitor-1'
	if kwargs.has_key('mask_run'):
		mask_run = kwargs.get('mask_run')
		print 'Using run ', kwargs.get('mask_run'),' for diag'
	else:
		mask_run=sample_run
	
	if kwargs.has_key('background'):
		reducer.background = kwargs.get('background')
		print 'Setting background option to ', kwargs.get('background')
	else:
		reducer.background = False
	
	if kwargs.has_key('fixei'):
		reducer.fix_ei = kwargs.get('fixei')
		print 'Setting fixei to ', kwargs.get('fixei')
	else:
		reducer.fix_ei = False
	
	if kwargs.has_key('save_format'):
		reducer.save_formats = kwargs.get('save_format')
		print 'Setting save format to ', kwargs.get('save_format')
	else:
		reducer.save_formats = ['.spe']
	#Set parameters for the run
	
	if kwargs.has_key('detector_van_range'):
		reducer.wb_integr_range = kwargs.get('detector_van_range')
		print 'Setting detector van int range to ', kwargs.get('detector_van_range')
	else:
		reducer.wb_integr_range=[20,100]
	#-------------DIAG------------------------
	if kwargs.has_key('bkgd_range'):
		background_range = kwargs.get('bkgd_range')
		print 'Setting background intergration to ', kwargs.get('bkgd_range')
	else:
		background_range=[15000,19000]
	
	if kwargs.has_key('tiny'):
		tinyval = kwargs.get('tiny')
		print 'Setting tiny ratelimit to ', kwargs.get('tiny')
	else:
		tinyval=1e-10
		
	if kwargs.has_key('large'):
		largeval = kwargs.get('large')
		print 'Setting large limit to ', kwargs.get('large')
	else:
		largeval=1e10
	
	if kwargs.has_key('diag_remove_zero'):
		sampzero = kwargs.get('diag_remove_zero')
		print 'Setting diag to reject zero backgrounds '
	else:
		sampzero =False
		
	if kwargs.has_key('diag_van_median_rate_limit_hi'):
		vanouthi = kwargs.get('diag_van_median_rate_limit_hi')
		print 'Setting diag_van_median_rate_limit_hi to ', kwargs.get('diag_van_median_rate_limit_hi')
	else:
		vanouthi=100
	
	if kwargs.has_key('diag_van_median_rate_limit_lo'):
		vanoutlo = kwargs.get('diag_van_median_rate_limit_lo')
		print 'Setting diag_van_median_rate_limit_lo to ', kwargs.get('diag_van_median_rate_limit_lo')
	else:
		vanoutlo=0.01
		
	if kwargs.has_key('diag_van_median_sigma_lo'):
		vanlo = kwargs.get('diag_van_median_sigma_lo')
		print 'Setting diag_van_median_sigma_lo to ', kwargs.get('diag_van_median_sigma_lo')
	else:
		vanlo=0.1
		
	if kwargs.has_key('diag_van_median_sigma_hi'):
		vanhi = kwargs.get('diag_van_median_sigma_hi')
		print 'Setting diag_van_median_sigma_hi to ', kwargs.get('diag_van_median_sigma_hi')
	else:
		vanhi=1.5
		
	if kwargs.has_key('diag_van_median_sigma'):
		vansig = kwargs.get('diag_van_median_sigma')
		print 'Setting diag_van_median_sigma to ', kwargs.get('diag_van_median_sigma')
	else:
		vansig=0.0
		
	if kwargs.has_key('diag_samp_median_sigma_lo'):
		samplo = kwargs.get('diag_samp_median_sigma_lo')
		print 'Setting diag_samp_median_sigma_lo to ', kwargs.get('diag_samp_median_sigma_lo')
	else:
		samplo=0.0
		
	if kwargs.has_key('diag_samp_median_sigma_hi'):
		samphi = kwargs.get('diag_samp_median_sigma_hi')
		print 'Setting diag_samp_median_sigma_hi to ', kwargs.get('diag_samp_median_sigma_hi')
	else:
		samphi=2.0
		
	if kwargs.has_key('diag_samp_median_sigma'):
		sampsig = kwargs.get('diag_samp_median_sigma')
		print 'Setting diag_samp_median_sigma to ', kwargs.get('diag_samp_median_sigma')
	else:
		sampsig=3.0
	
	if kwargs.has_key('bleed'):
		bleed_switch = kwargs.get('bleed')
		print 'Setting bleed ', kwargs.get('bleed')
	else:
		print 'bleed set to default'
	#---------------END of DIAG--------------------
	if kwargs.has_key('det_cal_file'):
		reducer.det_cal_file = kwargs.get('det_cal_file')
		reducer.relocate_dets = True
		print 'Setting detector calibration file to ', kwargs.get('det_cal_file')
	else:
		print 'Setting detector calibration to detector block info from ', sample_run
		reducer.det_cal_file =None
		reducer.relocate_dets = False
	
	if mtd.workspaceExists(str(sample_run))==True and kwargs.has_key('det_cal_file')==False:
		print 'For data input type: workspace detector calibration must be specified'
		print 'use Keyword det_cal_file with a valid detctor file or run number'
		return
	
	
	
	if kwargs.has_key('one2one'):
		reducer.map_file =None
		print 'one2one selected'
		
	else:
		reducer.map_file = map_file+'.map'

	reducer.energy_bins = rebin
	
	if float(str.split(rebin,',')[2])>=float(ei_guess):
		print 'error rebin range exceeds ei'
		return
	
	print 'output will be normalised to', reducer.normalise_method
	if (numpy.size(sample_run)) > 1 and kwargs.has_key('sum') and kwargs.get('sum')==True:
		#this sums the runs together before passing the summed file to the rest of the reduction
		#this circumvents the inbuilt method of summing which fails to sum the files for diag
		
		sumfilename=str(sample_run[0])+'sum'
		accum=sum_files(sumfilename, sample_run)
		#the D.E.C. tries to be too clever so we have to fool it into thinking the raw file is already exists as a workpsace
		RenameWorkspace(accum,inst_name+str(sample_run[0])+'.raw')
		sample_run=sample_run[0]
	
	if kwargs.has_key('hardmaskPlus'):
		HardMaskFile = kwargs.get('hardmaskPlus')
		print 'Use hardmask from ', HardMaskFile
		#hardMaskSpec=common.load_mask(HardMaskFile)
		#MaskDetectors(Workspace='masking',SpectraList=hardMaskSpec)
	else:
		HardMaskFile=None
	
	if kwargs.has_key('hardmaskOnly'):
		totalmask = kwargs.get('hardmaskOnly')
		print 'Using hardmask from ', totalmask
		#next stable version can replace this with loadmask algoritum
		specs=diag_load_mask(totalmask)
		CloneWorkspace(InputWorkspace=sample_run,OutputWorkspace='mask_wksp')
		MaskDetectors(Workspace='mask_wksp',SpectraList=specs)
		masking=mtd['mask_wksp']
	else:
	
		masking = reducer.diagnose(wb_run, 
			sample=mask_run,
			second_white = None,
			tiny=tinyval, 
			huge=largeval, 
			van_out_lo=vanoutlo,
			van_out_hi=vanouthi,
			van_lo=vanlo,
			van_hi=vanhi,
			van_sig=vansig,
			samp_zero=sampzero,
			samp_lo=samplo,
			samp_hi=samphi,
			samp_sig=sampsig,
			bkgd_range=background_range, 
			variation=1.1,
			print_results=True,
			bleed_test=bleed_switch,
			bleed_maxrate=rate,
			bleed_pixels=pixels,
			hard_mask=HardMaskFile)
			
	reducer.spectra_masks=masking
	fail_list=get_failed_spectra_list(masking)
	
	print 'Diag found ', len(fail_list),'bad spectra'
	
	#Run the conversion
	deltaE_wkspace = reducer.convert_to_energy(sample_run, ei_guess, wb_run)
	end_time=time.time()
	results_name=str(sample_run)+'.spe'
	
	ei= (deltaE_wkspace.getSampleDetails().getLogData("Ei").value)
	
	if mtd.workspaceExists('_wksp.spe-white')==True:
		DeleteWorkspace('_wksp.spe-white')
	
	if mtd.workspaceExists(results_name)==False:
		RenameWorkspace(deltaE_wkspace,results_name)
	
	print 'Incident energy found ',ei,' meV'
	print 'Elapsed time =',end_time-start_time, 's'
	#get the name that convert to energy will use
	
	
	RenameWorkspace(results_name,wksp_out)
	
	return mtd[wksp_out]
    
	
def abs_units_old(wb_run,sample_run,mono_van,wb_mono,samp_rmm,samp_mass,ei_guess,rebin,map_file,monovan_mapfile,**kwargs):
	"""	
	dgreduce.abs_units(wb_run,sample_run,mono_van,wb_mono,samp_rmm,samp_mass,ei_guess,rebin,map_file,monovan_mapfile,keyword arguments)
	
	Available keywords
	norm_method =[monitor-1],[monitor-2][Current]
	background  =False , True
	fixei 		=False , True
	save_format	=['.spe'],['.nxspe']
	detector_van_range			=[20,40] in mev
	
	bkgd_range	=[15000,19000]	:integration range for background tests
	
	second_white	- If provided an additional set of tests is performed on this. (default = None)
	hard_mask  		- A file specifying those spectra that should be masked without testing (default=None)
	tiny        	- Minimum threshold for acceptance (default = 1e-10)
	large        	- Maximum threshold for acceptance (default = 1e10)
	bkgd_range 		- A list of two numbers indicating the background range (default=instrument defaults)
	diag_van_median_rate_limit_lo  	- Lower bound defining outliers as fraction of median value (default = 0.01)
	diag_van_median_rate_limit_hi 	- Upper bound defining outliers as fraction of median value (default = 100.)
	diag_van_median_sigma_lo      	- Fraction of median to consider counting low for the white beam diag (default = 0.1)
	diag_van_median_sigma_hi      	- Fraction of median to consider counting high for the white beam diag (default = 1.5)
	diag_van_sig  - Error criterion as a multiple of error bar i.e. to fail the test, the magnitude of the
			  		difference with respect to the median value must also exceed this number of error bars (default=0.0)
	diag_remove_zero 				- If true then zeroes in the vanadium data will count as failed (default = True)
	diag_samp_samp_median_sigma_lo  - Fraction of median to consider counting low for the white beam diag (default = 0)
	diag_samp_samp_median_sigma_hi  - Fraction of median to consider counting high for the white beam diag (default = 2.0)
	diag_samp_sig  					- Error criterion as a multiple of error bar i.e. to fail the test, the magnitude of the"
			  						  difference with respect to the median value must also exceed this number of error bars (default=3.3)
	variation  		-The number of medians the ratio of the first/second white beam can deviate from
			   		the average by (default=1.1)
	bleed_test 		- If true then the CreatePSDBleedMask algorithm is run
	bleed_maxrate 	- If the bleed test is on then this is the maximum framerate allowed in a tube
	bleed_pixels 	- If the bleed test is on then this is the number of pixels ignored within the
				 	bleed test diagnostic
	print_results - If True then the results are printed to the screen
	
	diag_remove_zero =True, False (default):Diag zero counts in background range
	
	bleed=True , turn bleed correction on and off on by default for Merlin and LET
	
	sum	=True,False(default) , sum multiple files

	det_cal_file= a valid detector block file and path or a raw file. Setting this
				  will use the detector calibraion from the specified file NOT the
				  input raw file
	mask_run = RunNumber to use for diag instead of the input run number
	
	one2one =True, False :Reduction will not use a mapping file
	
	hardmaskPlus=Filename :load a hardmarkfile and apply together with diag mask
	
	hardmaskOnly=Filename :load a hardmask and use as only mask
	
	use_sam_msk_on_monovan=False This will set the total mask to be that of the sample run
	
	abs_units_van_range=[-40,40] integral range for absolute vanadium data
	
	"""	
	#available keywords
	#abs_units_van_range
	global reducer, rm_zero,inst_name,van_mass,bleed_switch,rate,pixels
	print 'DGreduce run for ',inst_name,'run number ',sample_run
	print 'Output will be in absolute units of mb/str/mev/fu'
	try:
		n,r=lhs('both')
		wksp_out=r[0]
	except:
		if sample_run == 0:
			#deal with the current run being parsed as 0 rather than 00000
			sample_run='00000'
		wksp_out=str(sample_run)+'.spe'
	
	start_time=time.time()
	
	if sample_run=='00000' and mtd.workspaceExists(inst_name+'00000.raw')==True:
		print 'Deleteing previous instance of temp data'
		DeleteWorkspace(inst_name+'00000.raw')
	
	if kwargs.has_key('norm_method'):
		reducer.normalise_method = kwargs.get('norm_method')
		print 'Setting normalisation method to ', kwargs.get('norm_method')
	else:
		reducer.normalise_method = 'monitor-1'
		
	if kwargs.has_key('mask_run'):
		mask_run = kwargs.get('mask_run')
		print 'Using run ', kwargs.get('mask_run'),' for diag'
	else:
		mask_run=sample_run
	
	if kwargs.has_key('background'):
		reducer.background = kwargs.get('background')
		print 'Setting background option to ', kwargs.get('background')
	else:
		reducer.background = False
	
	if kwargs.has_key('fixei'):
		reducer.fix_ei = kwargs.get('fixei')
		print 'Setting fixei to ', kwargs.get('fixei')
	else:
		reducer.fix_ei = False
	
	if kwargs.has_key('save_format'):
		reducer.save_formats = kwargs.get('save_format')
		print 'Setting save format to ', kwargs.get('save_format')
	else:
		reducer.save_formats = ['.spe']
	#Set parameters for the run
	
	if kwargs.has_key('detector_van_range'):
		reducer.wb_integr_range = kwargs.get('detector_van_range')
		print 'Setting detector van int range to ', kwargs.get('detector_van_range')
	else:
		reducer.wb_integr_range=[20,100]
	
	#######DIAG###########
	if kwargs.has_key('bkgd_range'):
		background_range = kwargs.get('bkgd_range')
		print 'Setting background intergration to ', kwargs.get('bkgd_range')
	else:
		background_range=[15000,19000]
	
	if kwargs.has_key('tiny'):
		tinyval = kwargs.get('tiny')
		print 'Setting tiny ratelimit to ', kwargs.get('tiny')
	else:
		tinyval=1e-10
		
	if kwargs.has_key('large'):
		largeval = kwargs.get('large')
		print 'Setting large limit to ', kwargs.get('large')
	else:
		largeval=1e10
	
	if kwargs.has_key('diag_remove_zero'):
		sampzero = kwargs.get('diag_remove_zero')
		print 'Setting diag to reject zero backgrounds '
	else:
		sampzero =False
		
	if kwargs.has_key('diag_van_median_rate_limit_hi'):
		vanouthi = kwargs.get('diag_van_median_rate_limit_hi')
		print 'Setting diag_van_median_rate_limit_hi to ', kwargs.get('diag_van_median_rate_limit_hi')
	else:
		vanouthi=100
	
	if kwargs.has_key('diag_van_median_rate_limit_lo'):
		vanoutlo = kwargs.get('diag_van_median_rate_limit_lo')
		print 'Setting diag_van_median_rate_limit_lo to ', kwargs.get('diag_van_median_rate_limit_lo')
	else:
		vanoutlo=0.01
		
	if kwargs.has_key('diag_van_median_sigma_lo'):
		vanlo = kwargs.get('diag_van_median_sigma_lo')
		print 'Setting diag_van_median_sigma_lo to ', kwargs.get('diag_van_median_sigma_lo')
	else:
		vanlo=0.1
		
	if kwargs.has_key('diag_van_median_sigma_hi'):
		vanhi = kwargs.get('diag_van_median_sigma_hi')
		print 'Setting diag_van_median_sigma_hi to ', kwargs.get('diag_van_median_sigma_hi')
	else:
		vanhi=1.5
		
	if kwargs.has_key('diag_van_median_sigma'):
		vansig = kwargs.get('diag_van_median_sigma')
		print 'Setting diag_van_median_sigma to ', kwargs.get('diag_van_median_sigma')
	else:
		vansig=0.0
		
	if kwargs.has_key('diag_samp_median_sigma_lo'):
		samplo = kwargs.get('diag_samp_median_sigma_lo')
		print 'Setting diag_samp_median_sigma_lo to ', kwargs.get('diag_samp_median_sigma_lo')
	else:
		samplo=0.0
		
	if kwargs.has_key('diag_samp_median_sigma_hi'):
		samphi = kwargs.get('diag_samp_median_sigma_hi')
		print 'Setting diag_samp_median_sigma_hi to ', kwargs.get('diag_samp_median_sigma_hi')
	else:
		samphi=2.0
		
	if kwargs.has_key('diag_samp_median_sigma'):
		sampsig = kwargs.get('diag_samp_median_sigma')
		print 'Setting diag_samp_median_sigma to ', kwargs.get('diag_samp_median_sigma')
	else:
		sampsig=3.0
	
	if kwargs.has_key('bleed'):
		bleed_switch = kwargs.get('bleed')
		print 'Setting bleed ', kwargs.get('bleed')
	else:
		print 'bleed set to default'
	#####diad end########
	
	if kwargs.has_key('det_cal_file'):
		reducer.det_cal_file = kwargs.get('det_cal_file')
		reducer.relocate_dets = True
		print 'Setting detector calibration file to ', kwargs.get('det_cal_file')
	else:
		print 'Setting detector calibration to detector block info from ', sample_run
		reducer.det_cal_file =None
		reducer.relocate_dets = False
	
	if mtd.workspaceExists(str(sample_run))==True and kwargs.has_key('det_cal_file')==False:
		print 'For data input type: workspace detector calibration must be specified'
		print 'use Keyword det_cal_file with a valid detctor file or run number'
		return
	
	
	if kwargs.has_key('one2one'):
		reducer.map_file =None
		print 'one2one selected'
		
	else:
		reducer.map_file = map_file+'.map'
	
	if kwargs.has_key('hardmaskPlus'):
		HardMaskFile = kwargs.get('hardmaskPlus')
		print 'Use hardmask from ', HardMaskFile
		#hardMaskSpec=common.load_mask(HardMaskFile)
		#MaskDetectors(Workspace='masking',SpectraList=hardMaskSpec)
	else:
		HardMaskFile=None
		
	reducer.energy_bins = rebin
	#monovan info
	fileName, fileExtension = os.path.splitext(monovan_mapfile)
	if (not fileExtension):
		monovan_mapfile=monovan_mapfile+'.map'
	reducer.abs_map_file =monovan_mapfile 

	if kwargs.has_key('abs_units_van_range'):
		reducer.monovan_integr_range = kwargs.get('abs_units_van_range')
		print 'Setting absolute units vanadium integratiOn range to ', kwargs.get('abs_units_van_range')
	else:
		reducer.monovan_integr_range=[-40,40]
	
	#reducer.van_rmm =50.94
	reducer.van_mass=van_mass
	#sample info
	reducer.sample_mass=samp_mass
	reducer.sample_rmm =samp_rmm
	
	print 'output will be normalised to', reducer.normalise_method
	if (numpy.size(sample_run)) > 1 and kwargs.has_key('sum') and kwargs.get('sum')==True:
		#this sums the runs together before passing the summed file to the rest of the reduction
		#this circumvents the inbuilt method of summing which fails to sum the files for diag
		
		sumfilename=str(sample_run[0])+'sum'
		accum=sum_files(sumfilename, sample_run)
		#the D.E.C. tries to be too clever so we have to fool it into thinking the raw file is already exists as a workpsace
		RenameWorkspace(accum,inst_name+str(sample_run[0])+'.raw')
		sample_run=sample_run[0]
	
	if kwargs.has_key('hardmaskOnly'):
		hardmask = kwargs.get('hardmaskOnly')
		print 'Use hardmask from ', hardmask
		masking=hardmask
	else:
	
		masking = reducer.diagnose(wb_run, 
			sample=mask_run,
			second_white = None,
			tiny=tinyval, 
			huge=largeval, 
			van_out_lo=vanoutlo,
			van_out_hi=vanouthi,
			van_lo=vanlo,
			van_hi=vanhi,
			van_sig=vansig,
			samp_zero=sampzero,
			samp_lo=samplo,
			samp_hi=samphi,
			samp_sig=sampsig,
			bkgd_range=background_range, 
			variation=1.1,
			print_results=True,
			bleed_test=bleed_switch,
			bleed_maxrate=rate,
			bleed_pixels=pixels,
			hard_mask=HardMaskFile)
		
		masking2 = reducer.diagnose(wb_mono, 
			sample=mono_van,
			second_white = None,
			tiny=tinyval, 
			huge=largeval, 
			van_out_lo=vanoutlo,
			van_out_hi=vanouthi,
			van_lo=vanlo,
			van_hi=vanhi,
			van_sig=vansig,
			samp_zero=sampzero,
			samp_lo=samplo,
			samp_hi=samphi,
			samp_sig=sampsig,
			bkgd_range=background_range, 
			variation=1.1,
			print_results=True,
			bleed_test=bleed_switch,
			bleed_maxrate=rate,
			bleed_pixels=pixels,
			hard_mask=HardMaskFile)
		
		total_mask=masking+masking2
	
	
	if kwargs.has_key('use_sam_msk_on_monovan') and kwargs.get('use_sam_msk_on_monovan')==True:
		print 'applying sample run mask to mono van'
		reducer.spectra_masks=masking
	else:
		reducer.spectra_masks=total_mask
		
	fail_list=get_failed_spectra_list('total_mask')
	
	
	print 'Diag found ', len(fail_list),'bad spectra'
		
	
	#Run the conversion
	deltaE_wkspace = reducer.convert_to_energy(sample_run, ei_guess, wb_run, mono_van,ei_guess,wb_mono)
	end_time=time.time()
	results_name=str(sample_run)+'.spe'
	ei= (deltaE_wkspace.getSampleDetails().getLogData("Ei").value)
	
	if mtd.workspaceExists('_wksp.spe-white')==True:
		DeleteWorkspace('_wksp.spe-white')
	
	
	print 'Incident energy found ',ei,' meV'
	print 'Elapsed time =',end_time-start_time, 's'
	#get the name that convert to energy will use
	
	if mtd.workspaceExists(results_name)==False:
		RenameWorkspace(deltaE_wkspace,results_name)
	RenameWorkspace(results_name,wksp_out)
	
	return mtd[wksp_out]

def abs_units(wb_run,sample_run,mono_van,wb_mono,samp_rmm,samp_mass,ei_guess,rebin,map_file,monovan_mapfile,**kwargs):
	"""	
	dgreduce.abs_units(wb_run,sample_run,mono_van,wb_mono,samp_rmm,samp_mass,ei_guess,rebin,map_file,monovan_mapfile,keyword arguments)
	
	Available keywords
	norm_method =[monitor-1],[monitor-2][Current]
	background  =False , True
	fixei 		=False , True
	save_format	=['.spe'],['.nxspe'],'none'
	detector_van_range			=[20,40] in mev
	
	bkgd_range	=[15000,19000]	:integration range for background tests
	
	second_white	- If provided an additional set of tests is performed on this. (default = None)
	hard_mask  		- A file specifying those spectra that should be masked without testing (default=None)
	tiny        	- Minimum threshold for acceptance (default = 1e-10)
	large        	- Maximum threshold for acceptance (default = 1e10)
	bkgd_range 		- A list of two numbers indicating the background range (default=instrument defaults)
	diag_van_median_rate_limit_lo  	- Lower bound defining outliers as fraction of median value (default = 0.01)
	diag_van_median_rate_limit_hi 	- Upper bound defining outliers as fraction of median value (default = 100.)
	diag_van_median_sigma_lo      	- Fraction of median to consider counting low for the white beam diag (default = 0.1)
	diag_van_median_sigma_hi      	- Fraction of median to consider counting high for the white beam diag (default = 1.5)
	diag_van_sig  - Error criterion as a multiple of error bar i.e. to fail the test, the magnitude of the
			  		difference with respect to the median value must also exceed this number of error bars (default=0.0)
	diag_remove_zero 				- If true then zeroes in the vanadium data will count as failed (default = True)
	diag_samp_samp_median_sigma_lo  - Fraction of median to consider counting low for the white beam diag (default = 0)
	diag_samp_samp_median_sigma_hi  - Fraction of median to consider counting high for the white beam diag (default = 2.0)
	diag_samp_sig  					- Error criterion as a multiple of error bar i.e. to fail the test, the magnitude of the"
			  						  difference with respect to the median value must also exceed this number of error bars (default=3.3)
	variation  		-The number of medians the ratio of the first/second white beam can deviate from
			   		the average by (default=1.1)
	bleed_test 		- If true then the CreatePSDBleedMask algorithm is run
	bleed_maxrate 	- If the bleed test is on then this is the maximum framerate allowed in a tube
	bleed_pixels 	- If the bleed test is on then this is the number of pixels ignored within the
				 	bleed test diagnostic
	print_results - If True then the results are printed to the screen
	
	diag_remove_zero =True, False (default):Diag zero counts in background range
	
	bleed=True , turn bleed correction on and off on by default for Merlin and LET
	
	sum	=True,False(default) , sum multiple files

	det_cal_file= a valid detector block file and path or a raw file. Setting this
				  will use the detector calibraion from the specified file NOT the
				  input raw file
	mask_run = RunNumber to use for diag instead of the input run number
	
	one2one =True, False :Reduction will not use a mapping file
	
	hardmaskPlus=Filename :load a hardmarkfile and apply together with diag mask
	
	hardmaskOnly=Filename :load a hardmask and use as only mask
	
	use_sam_msk_on_monovan=False This will set the total mask to be that of the sample run
	
	abs_units_van_range=[-40,40] integral range for absolute vanadium data
	
	mono_correction_factor=float User specified correction factor for absolute units normalisation
	"""	
	#available keywords
	#abs_units_van_range
	global reducer, rm_zero,inst_name,van_mass,bleed_switch,rate,pixels
	print 'DGreduce run for ',inst_name,'run number ',sample_run
	print 'Output will be in absolute units of mb/str/mev/fu'

	#reducer.van_rmm =50.94
	reducer.van_mass=van_mass
	#sample info
	reducer.sample_mass=samp_mass
	reducer.sample_rmm =samp_rmm
	print 'Using vanadium mass: ',van_mass
	print '        sample mass: ',samp_mass    
	print '        sample_rmm : ',samp_rmm 
    # check if mono-vanadium is provided as multiple files list or just put in brackets ocasionally
	if isinstance(mono_van,list):
             if len(mono_van)>1:
                 raise IOError(' Can currently work only with single monovan file but list supplied')
             else:
                mono_van = mono_van[0];

    
	try:
		n,r=lhs('both')
		wksp_out=r[0]
	except:
		if sample_run == 0:
			#deal with the current run being parsed as 0 rather than 00000
			sample_run='00000'
		wksp_out=str(sample_run)+'.spe'
	
	start_time=time.time()
	
	if sample_run=='00000' and mtd.workspaceExists(inst_name+'00000.raw')==True:
		print 'Deleteing previous instance of temp data'
		DeleteWorkspace(inst_name+'00000.raw')
	
	if kwargs.has_key('norm_method'):
		reducer.normalise_method = kwargs.get('norm_method')
		print 'Setting normalisation method to ', kwargs.get('norm_method')
	else:
		reducer.normalise_method = 'monitor-1'
		
	if kwargs.has_key('mask_run'):
		mask_run = kwargs.get('mask_run')
		print 'Using run ', kwargs.get('mask_run'),' for diag'
	else:
		mask_run=sample_run
	
	if kwargs.has_key('background'):
		reducer.background = kwargs.get('background')
		print 'Setting background option to ', kwargs.get('background')
	else:
		reducer.background = False
	
	if kwargs.has_key('fixei'):
		reducer.fix_ei = kwargs.get('fixei')
		print 'Setting fixei to ', kwargs.get('fixei')
	else:
		reducer.fix_ei = False
	
	if kwargs.has_key('save_format'):
		reducer.save_formats = kwargs.get('save_format')
		print 'Setting save format to ', kwargs.get('save_format')
	else:
		reducer.save_formats = ['.spe']
	#Set parameters for the run
	
	if kwargs.has_key('detector_van_range'):
		reducer.wb_integr_range = kwargs.get('detector_van_range')
		print 'Setting detector van int range to ', kwargs.get('detector_van_range')
	else:
		reducer.wb_integr_range=[20,100]
	
	#######DIAG###########
	if kwargs.has_key('bkgd_range'):
		background_range = kwargs.get('bkgd_range')
		print 'Setting background intergration to ', kwargs.get('bkgd_range')
	else:
		background_range=[15000,19000]
	
	if kwargs.has_key('tiny'):
		tinyval = kwargs.get('tiny')
		print 'Setting tiny ratelimit to ', kwargs.get('tiny')
	else:
		tinyval=1e-10
		
	if kwargs.has_key('large'):
		largeval = kwargs.get('large')
		print 'Setting large limit to ', kwargs.get('large')
	else:
		largeval=1e10
	
	if kwargs.has_key('diag_remove_zero'):
		sampzero = kwargs.get('diag_remove_zero')
		print 'Setting diag to reject zero backgrounds '
	else:
		sampzero =False
		
	if kwargs.has_key('diag_van_median_rate_limit_hi'):
		vanouthi = kwargs.get('diag_van_median_rate_limit_hi')
		print 'Setting diag_van_median_rate_limit_hi to ', kwargs.get('diag_van_median_rate_limit_hi')
	else:
		vanouthi=100
	
	if kwargs.has_key('diag_van_median_rate_limit_lo'):
		vanoutlo = kwargs.get('diag_van_median_rate_limit_lo')
		print 'Setting diag_van_median_rate_limit_lo to ', kwargs.get('diag_van_median_rate_limit_lo')
	else:
		vanoutlo=0.01
		
	if kwargs.has_key('diag_van_median_sigma_lo'):
		vanlo = kwargs.get('diag_van_median_sigma_lo')
		print 'Setting diag_van_median_sigma_lo to ', kwargs.get('diag_van_median_sigma_lo')
	else:
		vanlo=0.1
		
	if kwargs.has_key('diag_van_median_sigma_hi'):
		vanhi = kwargs.get('diag_van_median_sigma_hi')
		print 'Setting diag_van_median_sigma_hi to ', kwargs.get('diag_van_median_sigma_hi')
	else:
		vanhi=1.5
		
	if kwargs.has_key('diag_van_median_sigma'):
		vansig = kwargs.get('diag_van_median_sigma')
		print 'Setting diag_van_median_sigma to ', kwargs.get('diag_van_median_sigma')
	else:
		vansig=0.0
		
	if kwargs.has_key('diag_samp_median_sigma_lo'):
		samplo = kwargs.get('diag_samp_median_sigma_lo')
		print 'Setting diag_samp_median_sigma_lo to ', kwargs.get('diag_samp_median_sigma_lo')
	else:
		samplo=0.0
		
	if kwargs.has_key('diag_samp_median_sigma_hi'):
		samphi = kwargs.get('diag_samp_median_sigma_hi')
		print 'Setting diag_samp_median_sigma_hi to ', kwargs.get('diag_samp_median_sigma_hi')
	else:
		samphi=2.0
		
	if kwargs.has_key('diag_samp_median_sigma'):
		sampsig = kwargs.get('diag_samp_median_sigma')
		print 'Setting diag_samp_median_sigma to ', kwargs.get('diag_samp_median_sigma')
	else:
		sampsig=3.0
	
	if kwargs.has_key('bleed'):
		bleed_switch = kwargs.get('bleed')
		print 'Setting bleed ', kwargs.get('bleed')
	else:
		print 'bleed set to default'
	#####diad end########
    
	
	if kwargs.has_key('det_cal_file'):
		reducer.det_cal_file = kwargs.get('det_cal_file')
		reducer.relocate_dets = True
		print 'Setting detector calibration file to ', kwargs.get('det_cal_file')
	else:
		print 'Setting detector calibration to detector block info from ', sample_run
		reducer.det_cal_file =None
		reducer.relocate_dets = False
	
	if mtd.workspaceExists(str(sample_run))==True and kwargs.has_key('det_cal_file')==False:
		print 'For data input type: workspace detector calibration must be specified'
		print 'use Keyword det_cal_file with a valid detctor file or run number'
		return
	
	
	if kwargs.has_key('one2one'):
		reducer.map_file =None
		map_file = ""
		print 'one2one selected'
		
	else:
		fileName, fileExtension = os.path.splitext(map_file)
		if (not fileExtension):
			map_file =  map_file+'.map'
		reducer.map_file = map_file;
	
	if kwargs.has_key('hardmaskPlus'):
		HardMaskFile = kwargs.get('hardmaskPlus')
		print 'Use hardmask from ', HardMaskFile
		#hardMaskSpec=common.load_mask(HardMaskFile)
		#MaskDetectors(Workspace='masking',SpectraList=hardMaskSpec)
	else:
		HardMaskFile=None
		
	reducer.energy_bins = rebin
	#monovan info
	fileName, fileExtension = os.path.splitext(monovan_mapfile)
	if (not fileExtension):
		monovan_mapfile=monovan_mapfile+'.map'
	reducer.abs_map_file =monovan_mapfile 

	if kwargs.has_key('abs_units_van_range'):
		reducer.monovan_integr_range = kwargs.get('abs_units_van_range')
		print 'Setting absolute units vanadium integration range to: ', kwargs.get('abs_units_van_range')
	else:
		reducer.monovan_integr_range=[-40,40]

	
	
	print 'output will be normalised to', reducer.normalise_method
	if (numpy.size(sample_run)) > 1 and kwargs.has_key('sum') and kwargs.get('sum')==True:
		#this sums the runs together before passing the summed file to the rest of the reduction
		#this circumvents the inbuilt method of summing which fails to sum the files for diag
		
		sumfilename=str(sample_run[0])+'sum'
		accum=sum_files(sumfilename, sample_run)
		#the D.E.C. tries to be too clever so we have to fool it into thinking the raw file is already exists as a workpsace
		RenameWorkspace(accum,inst_name+str(sample_run[0])+'.raw')
		sample_run=sample_run[0]
	
	if kwargs.has_key('hardmaskOnly'):
		if (kwargs.get('hardmaskOnly')):		
			totalmask = kwargs.get('hardmaskOnly')
			print 'Using hardmask from ', totalmask
			#next stable version can replace this with loadmask algoritum
			specs=diag_load_mask(totalmask)
		else:
			specs=""
		  
		CloneWorkspace(InputWorkspace=sample_run,OutputWorkspace='mask_wksp')
		MaskDetectors(Workspace='mask_wksp',SpectraList=specs)
		masking =mtd['mask_wksp']
	else:
		print '########### Run diagnose for sample run ##############'
		masking = reducer.diagnose(wb_run, 
			sample=mask_run,
			second_white = None,
			tiny=tinyval, 
			huge=largeval, 
			van_out_lo=vanoutlo,
			van_out_hi=vanouthi,
			van_lo=vanlo,
			van_hi=vanhi,
			van_sig=vansig,
			samp_zero=sampzero,
			samp_lo=samplo,
			samp_hi=samphi,
			samp_sig=sampsig,
			bkgd_range=background_range, 
			variation=1.1,
			print_results=True,
			bleed_test=bleed_switch,
			bleed_maxrate=rate,
			bleed_pixels=pixels,
			hard_mask=HardMaskFile)
            
	
    
	if kwargs.has_key('use_sam_msk_on_monovan') and kwargs.get('use_sam_msk_on_monovan')==True:
		print 'applying sample run mask to mono van'
		reducer.spectra_masks=masking
		fail_list=get_failed_spectra_list(masking)		
	else:
		print '########### Run diagnose for monochromatic vanadium run ##############'
		masking2 = reducer.diagnose(wb_mono, 
			sample=mono_van,
			second_white = None,
			tiny=tinyval, 
			huge=largeval, 
			van_out_lo=vanoutlo,
			van_out_hi=vanouthi,
			van_lo=vanlo,
			van_hi=vanhi,
			van_sig=vansig,
			samp_zero=sampzero,
			samp_lo=samplo,
			samp_hi=samphi,
			samp_sig=sampsig,
			bkgd_range=background_range, 
			variation=1.1,
			print_results=True,
			bleed_test=bleed_switch,
			bleed_maxrate=rate,
			bleed_pixels=pixels,
			hard_mask=HardMaskFile)
				
		total_mask=masking+masking2			
		reducer.spectra_masks=total_mask		    
		fail_list=get_failed_spectra_list('total_mask')
	
	
	print 'Diag found ', len(fail_list),'bad spectra '
	
	
	
	#Run the conversion first on the sample
	deltaE_wkspace_sample = reducer.convert_to_energy(sample_run, ei_guess, wb_run)

	
	if kwargs.has_key('mono_correction_factor'):
		absnorm_factor=kwargs.get('mono_correction_factor')
		print 'Using supplied correction factor for absolute units'
	else:
		print '##### Evaluate the integral from the monovan run and calculate the correction factor ######'
		print '      Using absolute units vanadion integration range : ', reducer.monovan_integr_range                
        #now on the mono_vanadium run swap the mapping file
		reducer.map_file = monovan_mapfile	
		deltaE_wkspace_monovan = reducer.convert_to_energy(mono_van, ei_guess, wb_mono)
        
		absnorm_factor = getAbsNormalizationFactor(deltaE_wkspace_monovan.getName(),str(reducer.monovan_integr_range[0]),str(reducer.monovan_integr_range[1]))
		
	print 'Absolute correction factor =',absnorm_factor
	CreateSingleValuedWorkspace(OutputWorkspace='AbsFactor',DataValue=absnorm_factor)
	end_time=time.time()
	results_name=str(sample_run)+'.spe'
	ei= (deltaE_wkspace_sample.getSampleDetails().getLogData("Ei").value)
	
	if mtd.workspaceExists('_wksp.spe-white')==True:
		DeleteWorkspace('_wksp.spe-white')
	
	
	print 'Incident energy found for sample run ',ei,' meV'
	print 'Incident energy found for mono vanadium run ',ei,' meV'
	print 'Elapsed time =',end_time-start_time, 's'
	#get the name that convert to energy will use
	
	if mtd.workspaceExists(results_name)==False:
		RenameWorkspace(deltaE_wkspace_sample,results_name)
	RenameWorkspace(results_name,wksp_out)
	Divide(LHSWorkspace=wksp_out,RHSWorkspace='AbsFactor',OutputWorkspace=wksp_out)
	DeleteWorkspace('AbsFactor')
	return mtd[wksp_out]

def chunk(wb_run,sample_run,ei_guess,rebin,mapingfile,nchunk,**kwargs):
	"""
	chunk(wb_run,sample_run,ei_guess,rebin,mapfile,**kwargs)

	dgreduce.chunk(1000,10001,80,'-10,.1,70','mari_res', additional keywords as required)
	
	dgreduce.arb_units(1000,10001,80,'-10,.1,70','mari_res',fixei=True)
	available keywords
	norm_method =[monitor-1],[monitor-2][uamph]
	background=False , True
	fixei =False , True
	save_format=['.spe'],['.nxspe']
	bkgd_range=[15000,19000]
	detector_van_range=[20,40] in mev
	diag_sigma=3.0
	sum=True , sum multiple files
	bleed=True , turn bleed correction on and off on by default for Merlin and LET
	det_cal_file= a valid detector block file and path or a raw file. Setting this
				  will use the detector calibraion from the specified file NOT the
				  input raw file
	"""
	global reducer,rm_zero,inst_name,van_mass,bleed_switch,rate,pixels
	print 'DGreduce run for ',inst_name,'run number ',sample_run
	try:
		n,r=lhs('both')
		wksp_out=r[0]
	except:
		if sample_run == 0:
			#deal with the current run being parsed as 0 rather than 00000
			sample_run='00000'
		wksp_out=inst_name+str(sample_run)+'.spe'
		if kwargs.has_key('sum') and kwargs.get('sum')==True:
			wksp_out=inst_name+str(sample_run[0])+'sum'+'.spe'
		
	start_time=time.time()
	
	if sample_run=='00000' and mtd.workspaceExists(inst_name+'00000.raw')==True:
		print 'Deleteing previous instance of temp data'
		DeleteWorkspace(inst_name+'00000.raw')
	
	
	reducer.energy_bins = rebin
	
	mon_list1=reducer.ei_mon_spectra
	mon_list2=reducer.mon1_norm_spec
	mon_list1.append(mon_list2)
	#mon_list1.sort()
	print 'Monitors for this chunk are: ',mon_list1
	# monitors for merlin[69634,69638]
	
	if inst_name == 'MER':
		#number of spectrums per instrument and where the detectors start (i.e. 5 for mari but 1 for merlin)
		numspec=69632
		spectrum_start=1
	if inst_name == 'MAP':
		#number of spectrums per instrument and where the detectors start (i.e. 5 for mari but 1 for merlin)
		numspec=41472
		spectrum_start=1
	
	if kwargs.has_key('det_cal_file'):
		cal_file = kwargs.get('det_cal_file')	
	else:
		print 'Setting detector calibration to detector block info from ', sample_run
		
		reducer.det_cal_file =None
		reducer.relocate_dets = False
	nums=range(spectrum_start,numspec,nchunk)
	output_wkspName=wksp_out
	for i in nums:
		print '=========================================================================='
		print  'start spectra for this chunk',i
		chunk=range(i,i+nchunk)
		endIndex=nchunk-1
		if i+nchunk > numspec:
			chunk=range(i,numspec+1)
			endIndex=len(chunk)-1
		print 'end spectra for this chunk ', i+endIndex
		
		speclist=mon_list1+chunk
		#print speclist
		LoadRaw(Filename=wb_run,OutputWorkspace="wb_wksp",LoadLogFiles="0",SpectrumList=speclist)
	
		LoadRaw(Filename=sample_run,OutputWorkspace="run_wksp",LoadLogFiles="0",SpectrumList=speclist)
	
		tmp=arb_units("wb_wksp","run_wksp",ei_guess,rebin,'none_for_this_run_type',one2one=True,bleed=False,**kwargs)
		
		
		DeleteWorkspace("wb_wksp")
		DeleteWorkspace("run_wksp")
		#DeleteWorkspace("_wksp.spe")
		#DeleteWorkspace("_wksp.spe-white")
		
		if i == spectrum_start:
			#crop the workspace to remove the monitors, the workpsace seems sorted on specnumber so this is ok for instruments where the monitors are at the end of the 
			# spectrum list
			CropWorkspace(tmp,wksp_out,StartWorkSpaceIndex=0,EndWorkSpaceIndex=endIndex)
		else:
			CropWorkspace(tmp,tmp,StartWorkSpaceIndex=0,EndWorkSpaceIndex=endIndex)
			ConjoinWorkspaces(wksp_out,tmp,CheckOverlapping='0')
		print int(((float(i+endIndex))/float(numspec))*100),'% complete'
		print '===============================================================================' 
	
	GroupDetectors(InputWorkspace=output_wkspName,OutputWorkspace=output_wkspName,MapFile=mapingfile)

	
	
	print 'Elapsed time =',time.time()-start_time, 's'
	return mtd[wksp_out]


def get_failed_spectra_list(diag_workspace):
    """Compile a list of spectra numbers that are marked as
    masked in the given workspace

    Input:

     diag_workspace  -  A workspace containing masking
    """
    if type(diag_workspace) == str:
        diag_workspace = mtd[diag_workspace]

    if hasattr(diag_workspace, "getAxis") == False:
        raise ValueError("Invalid input to get_failed_spectra_list. "
                         "A workspace handle or name is expected")
        
    spectra_axis = diag_workspace.getAxis(1)
    failed_spectra = []
    for i in range(diag_workspace.getNumberHistograms()):
	try:
            det = diag_workspace.getDetector(i)
	except RuntimeError:
            continue
	if det.isMasked():
            failed_spectra.append(spectra_axis.spectraNumber(i))

    return failed_spectra

def sum_files(accumulator, files):
	if type(files) == list:
		tmp_suffix = '_plus_tmp'

		for filename in files:
			print 'Summing run ',filename,' to workspace ',accumulator
			temp = common.load_run(filename, force=False)
			if mtd.workspaceExists(accumulator)==False:
				#check for existance of output workpsace if false clone and zero
				print 'create output file'
				CloneWorkspace(temp,accumulator)
				CreateSingleValuedWorkspace(OutputWorkspace="tmp",DataValue="0",ErrorValue="0")
				Multiply(LHSWorkspace=accumulator,RHSWorkspace="tmp",OutputWorkspace=accumulator)
			Plus(accumulator, temp, accumulator)
		return accumulator


def diag_load_mask(hard_mask):
    """
    Load a hard mask file and return the
    list of spectra numbers it contains as a string

    Each line of the file specifies spectra to be masked by either specifying a range
    using a hypen or a single spectra number using a space as a delimiter between fields i.e.

    48897 - 49152
    50000
    60100-60105
    """
   
    mask_file = open(hard_mask)
    spectra_list = ""
    for line in mask_file:
        numbers = line.split()
        if len(numbers) == 0:
            continue
        # Any non-numeric character at the start of the line marks a comment, 
        # check the first character of the first word
        if not numbers[0][0].isdigit():
            continue
        num_cols = len(numbers)
        remainder = num_cols % 3
        group_end = num_cols - remainder
        # Jump in steps of 3 where there are whole blocks
        for index in range(0, group_end, 3):
            # Can either have a range specified with a "-" or single numbers
            n_i = numbers[index]
            n_ip1 = numbers[index+1]
            n_ip2 = numbers[index+2]
            # If there is a dash it will have to be the middle value
            if n_ip1 == '-':
                spectra_list += n_i + '-' + n_ip2
            else:
                spectra_list += n_i + "," + n_ip1 + ',' + n_ip2
            spectra_list += ","
        # Now deal with the remainder
        for index in range(group_end,num_cols):
            spectra_list += numbers[index] + ","
            
    if len(spectra_list) < 1:
        mantid.sendLogMessage('Only comment lines found in mask file ' + hard_mask)
        return ''
    # Return everything after the very first comma we added in the line above
    return spectra_list.rstrip(',')
    
def getAbsNormalizationFactor(deltaE_wkspace,min,max):
    """
    get absolute normalization factor for monochromatic vanadium
    Inputs:
    @param: deltaE_wkspace  -- the name (string) of monovan workspace, converted to energy
    @param: min             -- the string representing the minimal energy to integrate the spectra
    @param: max             -- the string representing the maximal energy to integrate the spectra    
    
    @returns the value of monovan absolute normalization factor. 
             deletes monovan workspace (deltaE_wkspace) if abs norm factor was calculated successfully
    """
    global reducer,van_mass
    Integration(InputWorkspace=deltaE_wkspace,OutputWorkspace='van_int',RangeLower=min,RangeUpper=max,IncludePartialBins='1')
    input_ws = mtd[deltaE_wkspace]
    ei_monovan = input_ws.getSampleDetails().getLogData("Ei").value
    data_ws=mtd['van_int']
    nhist = data_ws.getNumberHistograms()
   #print nhist

    average_value = 0.0
    weight_sum = 0.0
    ic=0;
    izerc=0;
    for i in range(nhist):
        try:
            det = data_ws.getDetector(i)
        except Exception:
            continue
        if det.isMasked():
            continue

        signal = data_ws.readY(i)[0]
        error = data_ws.readE(i)[0]
        
        if signal != signal:	#ignore NaN
            continue
        if (error<=0):          # ignore Inf (0 in error are probably 0 in sign
            izerc+=1
            continue
            
        weight = 1.0/error
        average_value += signal * weight
        weight_sum += weight
        ic += 1
        #print 'signal value =' ,signal
        #print 'error value =' ,error        
        #print 'average ',average_value        
   #---------------- Loop finished
   
    integral_monovan=average_value / weight_sum
    van_multiplier = (float(reducer.van_rmm)/float(van_mass))
    absnorm_factor = integral_monovan * van_multiplier
    #print 'Monovan integral :' ,integral_monovan        
    
    if ei_monovan >= 200.0:
        xsection = 421.0
    else:
        xsection = 400.0 + (ei_monovan/10.0)

    absnorm_factor /= xsection
    sample_multiplier = (float(reducer.sample_mass)/float(reducer.sample_rmm))
    absnorm_factor= absnorm_factor *sample_multiplier
    
    if (absnorm_factor !=absnorm_factor)|(izerc!=0):    # It is an error, print diagnostics:
        if (absnorm_factor !=absnorm_factor):
            print '-----------> Absolute normalization factor is NaN <--------------------------------------------------'
        else:
            print '-----------> Warning, Monovanadium has zero spectra <------------------------------------------------'		
        print '-----------> Processing workspace: ',deltaE_wkspace
        print '-----------> Monovan Integration range : min=',min,' max=',max
        print '-----------> Calculated: ',ic,' spectra with average value: ',average_value, 'and total weight: ',weight_sum
        print '-----------> Dropped: ',izerc,' empty spectra'
        print '-----------> Van multiplier: ',van_multiplier,'  sample multiplier: ',sample_multiplier, 'and xsection: ',xsection		
        print '-----------------------------------------------------------------------------------------------------'	
    else:
        DeleteWorkspace(deltaE_wkspace)
    DeleteWorkspace(data_ws)
    return absnorm_factor
    
