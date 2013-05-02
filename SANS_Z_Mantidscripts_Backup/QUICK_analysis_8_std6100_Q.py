###########################################################################
#  14/8/9, 17/05/10 attempt to get neutrons PER UAMPHR vs  Q
#  no normalisation to monitor, no transmission correction
#  compare SANS2d 914 & LOQ 53335
# 27/4/10 NOTE line 186 normalises to current
#  NOTE edit beam centre coords and det offsets etc. lines ca 74 & 121
# #########################################################################
#  4/8/9 was import LOQFunctions  now
import SANSUtility_RKH
import math

# set this up for cycle 09/2, for which the spectra are now by rows, right to left from top right, and start at 9 after 8 monitors !
# i.e. still not as we really want them.  
#  ( see version 7 for cycle 09/1, spectra by columns from bottom right going upwards)
#path="c:/mantidinstall/data/"
# thus far SANS2d transmission files are same as for sans runs, but will change before too long
SANS2d_trans_inst_def_file="c:/mantidinstall/instrument/SANS2d_definition.xml"
#path="//isis/inst$/NDXSANS2D/Instrument/data/cycle_09_5/"
path="c:/mantidinstall//data/"
#path = '/home/dmn58364/Desktop/SANS2D_Data/'
# get whole raw file once only
#instr_dir = '/home/dmn58364/Mantid/trunk/Test/Instrument/'
#instr_dir = '/home/dmn58364/Mantid/trunk/Test/Instrument/'
#
# WHY is there no ctrl/F find function in the Mantid Python window ?
# NEED assemble these string properly, with correct number of leading zeros etc
#                                    S                            T
#  normal                       769                            772      SDS rot
#   -1000 +100  -12deg      774                            773        cell
#                                                                    771  DB
#  -860  +260 -10deg         801                            804      SDS
#                                    802                            805       cell
#                                                                     806
#    tk49                         914                            909
#   can                            915                            910
#                                                                     908
run_integer = 6100
prefix = 'SANS2D'
run_no = '0000'+str(run_integer)
extRaw = '.raw'
ext = '.raw'
# choose SANS2d detector HERE ++++++++++++++++++++++++++++++++++++++++++
detbank = 'rear-detector'
# This indicates whether a 1D or a 2D analysis is performed
correction_type = '1D'
#
# First the worksapces
scatter_sample = str(run_integer)+'_sans_raw'
scatter_can = ''
trans_sample = ''
trans_can = ''
direct_sample = ''
#
# set up TRANSMISSION RUNS
# sample trans run  CHECK if name is OK here
trans_sample_integer = 6096
trans_sample = str(trans_sample_integer) + '_trans_raw'
# direct beam trans run
direct_sample_integer = 6095
direct_sample = str(direct_sample_integer)+'_trans_raw'
#
#  instr_name and specmin (see below) determine which beam line and which detector bank - need a better way to do this 
#  as most people won't know the spectrum numbers
#
# LOG files for SANS2d ( on ndxsans2d//c/data/  ) will have these encoder readings  ============================================
# or will all this lot be found in NEXUS file ?
#     selog.front_det_z.value   etc
#     selog.rear_det_z.value
Front_Det_z = 2500.
Front_Det_x = -1050.
Front_Det_Rot = -21.
Rear_Det_z = 4000.
# Rear_Det_x  will be needed to calc relative X translation of front detector 
Rear_Det_x = 120.
#
# MASK file stuff ==========================================================
# correction terms to SANS2d encoders - store in MASK file
Rear_Det_z_corr = 58.
Rear_Det_x_corr = -16
Front_Det_z_corr = 47.
Front_Det_x_corr = -44.
Front_Det_y_corr = -20.
Front_Det_Rot_corr = 0.0
#
# Now the mask string, CORRECT expects one, else stops with an error !
# HOW include more CHECK
#  NEED  MASKF v0>v2   etc for front detector
# BEWARE spaces in these strings cause errors !
MaskString = 'h0,h127,v0,v127'
MaskStringRear = 'h0,h190>h191,v0,v191'
MaskStringFront = 'h0,h190>h191,h141>h160,v0,v191,h130>h167+v0>v6'
#  coords of SANS2d rear or LOQ det beam centre, from MASK file (though check from spectrum number for SANS2d below) 
# LOQ xcentre = 324/1000.
#        ycentre = 329/1000.
#	 backmon_start = '31000'
#	 backmon_end = '39000'
#	direct_beam_file = path + 'DIRECT.041'
#   SANS2d direct beam file
direct_beam_file = 'C:/mantidinstall/data/DIRECT_RUN524.dat'
#
# Default Analysis values  1.5 to 15.0 likely OK in direct beam, but try 2-14 for now
wav1 = 2.2
wav2 = 10.
rmin = 41./1000.
# 12/8/9 approx correction from quick fit to TK49 run 914 less can 915
rescale = 0.429*100.0
# rmax > 10m is ignored in CORRECT below
rmax = 20000./1000.
# check appropriate dwav for 0.75msec on current runs  T=0.253D.Lam     D~23m lam=0.13ang, 
# reduce Q2 for longer distances
dwav = 0.125
q1 = 0.006
q2 = 1.0
dq = -.08
qxy2 = 0.4
dqxy = 0.008
#
# SANS2d run 519 (rear at x=+100)  0.0842, -.1965 ( 0.15m LESS WITH NEW 21/7/9 INSTR DESCRIPTION)
#  same for run 626 but now spectrum 25848 due to new spectrum table at start if cycle 9/02
#xcentre = 0.0842
# from mask_101D
xcentre = 0.1167
# SANS2d run 242 ( rear at x=+200)
#xcentre = 0.1863
ycentre = -0.1703
#  NOTE these may need to vary with chopper phase and detector distance !
backmon_start = '88000'
backmon_end = '98000'
# sample dependent corrections to go in MASK file:
# SANS2d changer samples 37mm before shutter, which is 81 mm after notional sample, so add 44mm
# how do we get in the sample-detector distance for LOQ ?  This gets set zero below after first use
# to avoid repeating it.
# rotating rack, run 769, is 139 before shutter, so subtract 58 mm
Sample_z_corr = +28.
#
# hardware, rotation radius of front det
Front_Det_Radius = 306.
# how many monitors before the rear detector on SANS2d
monstart = 8
if run_integer < 568:
	monstart = 0
#
# SANS2d numbers in default instrument description, detectors 4m from sample and front det at x=1.1m
# BETTER to get the absolute z value from the file itself ? - but HOW ?
front_det_default_sd_m = 4.0
front_det_default_x_m = 1.1
rear_det_default_sd_m = 4.0
#

# The instrument dependent information
instr_name = prefix
# Which bank are we looking at
#detbank = ''
xshift = ''
yshift = ''
#
# sort out the monitor spectrum numbers
if instr_name == 'SANS2D':
#	The UDET hardcoding in CalculateTransmission was in the C++ code itself. 
#I've now changed that to add properties where you can pass in the appropriate detector id's. The property names are "IncidentBeamMonitor" 
#and "TransmissionMonitor", with the defaults being the LOQ values so nothing should be broken that wasn't already!

	# assume SANS2d, spectra are as per DAE, subtract one later
	monitorspectrum = 2
	# for transmissions
	transmonitorspectrum = 2
	transdetspectrum = 3
	usePreviousTrans = ' '
	# these are detector numbers, for calculateTransmission, which so far are always 2 & 3 (even for cycle 09/1) - may need 4 one day for beam stop transmission det
	UdetMonitor = 2
	UdetTrans = 3
	if run_integer < 568:
		monitorspectrum = 73730
	if direct_sample_integer < 568:
		transmonitorspectrum = 73730
		transdetspectrum = 73731
	dimension=192
 #========  load up the SANS2d transmission runs
	if trans_sample <> ' ':
		LoadRaw(path + prefix + '0000'+str(trans_sample_integer) + ext, trans_sample, SpectrumMin=transmonitorspectrum,SpectrumMax=transdetspectrum)
	if direct_sample <> ' ':
		LoadRaw(path + prefix + '0000'+str(direct_sample_integer) + ext, direct_sample, SpectrumMin=transmonitorspectrum,SpectrumMax=transdetspectrum)
#LoadISISNexus(path + prefix + run_no + ext, run_no.lstrip('0') + '_sans_raw')
#
# ========= load the SANS raw file, 
# WHY with raw file on local hard drive the time taken to load it varies hugely - is it fighting for memory ?
LoadRaw(path + prefix + run_no + extRaw, run_no.lstrip('0') + '_sans_raw')
#
NormaliseByCurrent(run_no.lstrip('0') + '_sans_raw',run_no.lstrip('0') + '_sans_raw')
#
# Sort out the input. The tags get replaced by information in the LOQ GUI
#
# pick spectrum for centre, enter number below
wksp = mantid.getMatrixWorkspace(scatter_sample)
instrument = wksp.getInstrument()
# run 242         det = wksp.getDetector(25401)
# NEED to pass this spectrum number into the CORRECT loop below ...
# run 519         det = wksp.getDetector(21561)
# run 626   25848
det = wksp.getDetector(11397)
# then check the coords for the shift 
out=str(det.getPos())
# NEED to read X&Y coords from string "out" here 
print "det 11397 is  "+out 
#det = wksp.getDetector(73728)
#out=str(det.getPos())
#print "det 73728 is  "+out 
# WHY, with run 519, centre at 21561,  having run around a couple of times does PICK in show intstrument now give a beam centre pixel as detector  as ~ 73728 ?
# what's more it returns coords [ 0, 0, 7.217]     whilst det 21561 still returns [0,  0,  25.339 ]
#
print "centre_rear = ( "+str(xcentre)+" , "+str(ycentre)+" )"
# debug stuff
# x coord of central pixel of front detector from instrument description
# NOTE this spectrum number will change when we go to proper wiring table
#  i = str(192*192*1.5)  =55296 less 1 plus 8 = 55303
# WHY can't I program this here ?  Will only take exact integer, no expressions !
#det = wksp.getDetector(55303)
#out=str(det.getPos())
#print "centre det 55296-1+8 is  "+out 
#det = wksp.getDetector(55207)
#out=str(det.getPos())
#print "edge det 55200-1+8 is  "+out 
#det = wksp.getDetector(36872)
#out=str(det.getPos())
#print "corner det 192^2 +1-1+8 is  "+out 
#det = wksp.getDetector(8)
#out=str(det.getPos())
#print "corner det 1-1+8 is  "+out 
#det = wksp.getDetector(36959)
#out=str(det.getPos())
#print "det 36960-1 is  "+out 
#

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
if detbank == 'front-detector':
	specmin = 192*192 + 1 + monstart
	specmax = 192*192*2 + monstart
if detbank == 'rear-detector':
	specmin = 1 + monstart
	specmax = 192*192 +monstart
#
# NOW SET UP THE BEAM LINE   22/7/9 CHANGE TO DO BOTH DETS ON SANS2D ONCE ONLY
# (ELSE HAD TO BE VERY CAREFUL LOOPING BACK AROUND TO DO THIS ONLY ONCE  FOR FRONT & ONCE FOR REAR DETECTOR)
# 
if instr_name == 'LOQ':
	monitorspectrum = 2
	dimension = 128
	# detbank names have to match those in instrument descriptions
	if specmin < 16387:
		detbank = 'main-detector-bank'
	else:
		detbank = 'HAB'	
	xshift = str( (317.5/1000.0) - xcentre)
	yshift = str( (317.5/1000.0) - ycentre)
	#direct_beam_file = path + 'DIRECT.041'
	direct_beam_file = 'C:/mantidinstall/data/DIRECT.041'
if instr_name == 'SANS2D':
	# corrected 19/6/9 was wrong way round for current SANS2d runs
	det = wksp.getDetector(36872)
	out=str(det.getPos())
	print "corner det 192^2 +1-1+8 is  "+out 
	det = wksp.getDetector(8)
	out=str(det.getPos())
	print "corner det 1-1+8 is  "+out 
	# =========================   front det =================================
	detbank_now ='front-detector'
	rotateDet = str( -Front_Det_Rot -  Front_Det_Rot_corr)
	print "FRONT rotateDet ="+rotateDet
	# was going to rotate and translate in CORRECT below, but for now will do here.
	# The actual roation is at a point 306mm behind the active plane of the detector.
	# then translate to correct relative position to rear detector, note due to packing plate under rear det is 20mm lower !
	# this for now found by using dead reckoning from rear detector using encoder values and a correction.
	# note SANS2d x encoders are -ve towards POLREF, but MANTID is +ve ,  "angle" is usually negative
	RotRadians = math.pi*(Front_Det_Rot + Front_Det_Rot_corr)/180.
	# this roates about the internal axis system of the component, the translation below corrects for the actual roation centre.
	RotateInstrumentComponent(scatter_sample,detbank_now,X="0.",Y="1.0",Z="0.",Angle=rotateDet)
	wksp = mantid.getMatrixWorkspace(scatter_sample)
	instrument = wksp.getInstrument()
	det = wksp.getDetector(55303)
	out=str(det.getPos())
	print "centre det 55296-1+8 is  "+out 
	det = wksp.getDetector(55207)
	out=str(det.getPos())
	print "edge det 55200-1+8 is  "+out 
	det = wksp.getDetector(36872)
	out=str(det.getPos())
	print "corner det 192^2 +1-1+8 is  "+out 
	det = wksp.getDetector(8)
	out=str(det.getPos())
	print "corner det 1-1+8 is  "+out 
	#		if run_integer < 676:
	# TEMPORARY TRANSLATION  for the flipped X&Y in cycle 09/1 & start of 09/2, due to 150mm shift in Y that is now in X !  NB  rotradians is NEGATIVE
	# 21/7/9 new instr description should not need this
	#			print "temp translate"
	#			MoveInstrumentComponent(scatter_sample, detbank, X = str(0.150*(1.0-cos(RotRadians))), Y = "0.0", Z = str(0.150*sin(-RotRadians)), RelativePosition="1")
	wksp = mantid.getMatrixWorkspace(scatter_sample)
	instrument = wksp.getInstrument()
	det = wksp.getDetector(55303)
	out=str(det.getPos())
	print "centre det 55296-1+8 is  "+out 
	det = wksp.getDetector(55207)
	out=str(det.getPos())
	print "edge det 55200-1+8 is  "+out 
	det = wksp.getDetector(36872)
	out=str(det.getPos())
	print "corner det 192^2 +1-1+8 is  "+out 
	det = wksp.getDetector(8)
	out=str(det.getPos())
	print "corner det 1-1+8 is  "+out 
#	print RotRadians
#	print sin(RotRadians)
	xshift = str( (Rear_Det_x + Rear_Det_x_corr - Front_Det_x - Front_Det_x_corr + Front_Det_Radius*math.sin(RotRadians ) )/1000. - front_det_default_x_m -xcentre)
	yshift = str( Front_Det_y_corr /1000.  - ycentre)
	# default in instrument description is 23.281m - 4.000m from sample at 19,281m !
	zshift = str( (Front_Det_z + Front_Det_z_corr + Front_Det_Radius*(1 - math.cos(RotRadians )) )/1000.  -front_det_default_sd_m)
	print "front translate"
	print "shifts = " + xshift + " , " + yshift + " , " + zshift
	MoveInstrumentComponent(scatter_sample, detbank_now, X = xshift, Y = yshift, Z = zshift, RelativePosition="1")
# =========== ================   rear det ===============================	
	detbank_now ='rear-detector'
	# this for now found by picking a spectrum for beam centre on show instrument, see above, will need to iterate to fine tune this
	xshift = str( -xcentre)
	yshift = str( -ycentre)
	zshift = str((Rear_Det_z + Rear_Det_z_corr)/1000. - 4.000)
	print "rear translate"
	print "shifts = " + xshift + " , " + yshift + " , " + zshift
	MoveInstrumentComponent(scatter_sample, detbank_now, X = xshift, Y = yshift, Z = zshift, RelativePosition="1")
	# more debug stuff ....
	wksp = mantid.getMatrixWorkspace(scatter_sample)
	instrument = wksp.getInstrument()
	det = wksp.getDetector(55303)
	out=str(det.getPos())
	print "centre det 55296-1+8 is  "+out 
	det = wksp.getDetector(55207)
	out=str(det.getPos())
	print "edge det 55200-1+8 is  "+out 
	det = wksp.getDetector(36872)
	out=str(det.getPos())
	print "corner det 192^2 +1-1+8 is  "+out 
	det = wksp.getDetector(8)
	out=str(det.getPos())
	print "corner det 1-1+8 is  "+out 
#
# ==================                move SANS2d sample 
# Component positions  - tested, +0.5m moves peak to higher Q
	MoveInstrumentComponent(scatter_sample, 'some-sample-holder', Z = str(Sample_z_corr/1000.), RelativePosition="1")	
# NEED to pull in the absolute sample position from instrument description and use that ?	
# meanwhile zero the correction here to avoid keep repeating it !
#	Sample_z_corr =0.0
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# all instruments


### Main correction function ###
visitcount = 0
def Correct(sampleWS, transWS, resultWS, suffix):
	'''Performs the data reduction steps'''
	# mop up any debug left overs from previous time through, incrementing visitcount here not allowed ?
#	if visitcount > 0:
#		mantid.deleteWorkspace(monitorWS)
#	        mantid.deleteWorkspace(trans_tmp_out)
#              mantid.deleteWorkspace(direct_tmp_out)
        # Set up the workspaces to work with
	#monitorWS = "Monitor"
	# Get the monitor ( StartSpectrum is off by one with cropworkspace)
	#CropWorkspace(sampleWS, OutputWorkspace = monitorWS, StartWorkspaceIndex=(monitorspectrum - 1), EndWorkspaceIndex=(monitorspectrum - 1))
	#if instr_name == 'LOQ':
		#RemoveBins(monitorWS, monitorWS, '19900', '20500', Interpolation="Linear")
	#FlatBackground(monitorWS, monitorWS, '0', backmon_start, backmon_end)

	# Get the bank we are looking at
	CropWorkspace(sampleWS, resultWS, StartWorkspaceIndex=(specmin - 1),EndWorkspaceIndex=(specmax - 1))
	
	# Mask the corners and beam stop if radius parameters are given:
	# CHECK   -  BETTER ONLY USE FOR REAR DET ?   - after rotations for FRONT det ????
	#  CHECK results log gives "error parsing" messages related to this ???
	print 'rmin = '+str(rmin)
	print 'rmax = '+str(rmax)
	if rmin > 0.0: 
		SANSUtility_RKH.MaskInsideCylinder(resultWS, rmin)

	if rmax > 0.0 and rmax < 10.0 :
		SANSUtility_RKH.MaskOutsideCylinder(resultWS, rmax)
	
	# Mask other requested spectra
	if instr_name == 'LOQ':
		speclist = SANSUtility_RKH.ConvertToSpecList(MaskString, specmin, dimension)
	if instr_name == 'SANS2D':
		if detbank == 'front-detector':
			print MaskStringFront
			speclist = SANSUtility_RKH.ConvertToSpecList(MaskStringFront, specmin, dimension)
		if detbank == 'rear-detector':
			print MaskStringRear
			speclist = SANSUtility_RKH.ConvertToSpecList(MaskStringRear, specmin, dimension)
	#		
	SANSUtility_RKH.MaskBySpecNumber(resultWS, speclist)
	#
	# Convert all of the files to wavelength and rebin
	# ConvertUnits does have a rebin option, but it's crude. In particular it rebins on linear scale.
	wavbin = str(wav1) + "," + str(dwav) + "," + str(wav2)

	#ConvertUnits(monitorWS, monitorWS, "Wavelength")
        #Rebin(monitorWS, monitorWS,wavbin)
	ConvertUnits(resultWS,resultWS,"Wavelength")
	Rebin(resultWS,resultWS,wavbin)
	      
	# Correct by incident beam monitor
	# At this point need to fork off workspace name to keep a workspace containing raw counts
        # NORMALLY move to tmpWS here, then back to resultWS at end, try equate them instead here!
        # tmpWS = "temporary_workspace"
        tmpWS = str(resultWS)
	print "tmpWS = "+tmpWS
	#Divide(resultWS, monitorWS, tmpWS)
	#mantid.deleteWorkspace(monitorWS)
 	# Correct by transmission
 	if instr_name == 'skipLOQ' and transWS and direct_sample:
                # Change the instrument definition to the correct one in the LOQ case 
                LoadInstrument(transWS, instr_dir + "/LOQ_trans_Definition.xml")
                LoadInstrument(direct_sample, instr_dir + "/LOQ_trans_Definition.xml")
		trans_tmp_out = LOQFunctions.SetupTransmissionData(transWS, instr_name, '1,2', backmon_start, backmon_end, wavbin) 
		direct_tmp_out = LOQFunctions.SetupTransmissionData(direct_sample, instr_name, '1,2', backmon_start, backmon_end, wavbin)              
                trans_ws = transWS.split('_')[0] + "_trans_" + suffix
                CalculateTransmission(trans_tmp_out,direct_tmp_out,trans_ws)
                Divide(tmpWS, trans_ws, tmpWS)
 	if instr_name == 'SANS2D'and transWS and direct_sample:
		trans_ws = transWS.split('_')[0] + "_trans_" + suffix
		if usePreviousTrans == ' ':
			print "calculate SANS2d trans"
			# Change the instrument definition to the correct one in the LOQ case 
			LoadInstrument(transWS, SANS2d_trans_inst_def_file)
			LoadInstrument(direct_sample,SANS2d_trans_inst_def_file)
			# I think '0,1' is a spectrum list, within the transWS which only has the monitor & transm detector
			trans_tmp_out = SANSUtility_RKH.SetupTransmissionData(transWS, instr_name, '0,1', backmon_start, backmon_end, wavbin) 
			direct_tmp_out = SANSUtility_RKH.SetupTransmissionData(direct_sample, instr_name, '0,1', backmon_start, backmon_end, wavbin)              
			
			# note this uses detector numbers, which for monitors on SANS2d are always 2 & 3 ( even for cycle 09/1 ), though might be 2 & 4 when we get a beam
			# stop transmission detector.
			CalculateTransmission(trans_tmp_out,direct_tmp_out,trans_ws,UdetMonitor,UdetTrans,wav1,wav2,True)
			mantid.deleteWorkspace(trans_tmp_out)
			mantid.deleteWorkspace(direct_tmp_out)
		else:
			print "using previous trans"
			Rebin(usePreviousTrans,trans_ws,wavbin)
                # Apply the correction
                Divide(tmpWS, trans_ws, tmpWS)
	#transmission, efficiency & direct beam are all simple functions of wavelength, so could multiply them together and do all at once,
	# though keeping track of the errors is fun
	# Correct for efficiency
	#NEED DIRECT BEAM FOR SANS2D (IF SAME FOR BOTH DETECTORS, THEN COULD TRY REDUCE BOTH SIMULTANEOULSY, BUT
	#WOULD FIRST NEED MASK COMMANDS TO IDENTIFY WHICH DETECTOR THEY APPLY TO.
	#CorrectToFile(tmpWS, direct_beam_file, tmpWS, "Wavelength", "Divide")

	# Correct for sample/Can volume
	#SANSUtility_RKH.ScaleByVolume(tmpWS, rescale);

	# Steps differ here depending on type
	if correction_type == '1D':
		# Convert to Q
		#ConvertUnits(tmpWS,tmpWS,"MomentumTransfer")
		ConvertUnits(resultWS,resultWS,"MomentumTransfer")
		
		# Need to mark the workspace as a distribution at this point to get next rebin right
		#ws = mantid.getMatrixWorkspace(tmpWS)
		#ws.isDistribution(True)
		# Calculate the solid angle corrections
		#NOTE STILL NEED TO ALLOW ALTERNATIVE USE OF FLAT_CELL (flood source data) FOR LOQ
		#SolidAngle(tmpWS,"solidangle")

		# Rebin to desired Q bins
		q_bins = str(q1) + "," + str(dq) + "," + str(q2)
		Rebin(resultWS, resultWS, q_bins)
		#Rebin(tmpWS, tmpWS, q_bins)
		#RebinPreserveValue("solidangle", "solidangle", q_bins)

		# Sum all spectra
		SumSpectra(resultWS,resultWS)
		#SumSpectra(tmpWS,tmpWS)
		#SumSpectra("solidangle","solidangle")
	
		# Correct for solidangle
		#Divide(tmpWS,"solidangle",tmpWS)
		#mantid.deleteWorkspace("solidangle")

		# Now put back the fractional error from the raw count workspace into the result
		#PoissonErrors(tmpWS, resultWS, resultWS)
		
	else:
		# Run 2D algorithm
#  HOW get 2d to work for SANS2d ?
		Qxy(tmpWS, resultWS, str(qxy2), str(dqxy))
		
# NEED really to crop off all leading & trailin nan,  this function sets them to zero 
	ReplaceSpecialValues(InputWorkspace=resultWS,OutputWorkspace=resultWS,NaNValue="0",InfinityValue="0")
# IS there a functionto check if a workspace exists ? as deleteWorkspace throws an error if it does not exist ....
# just do
#	if tmpWS:
#	mantid.deleteWorkspace(tmpWS)
#	visitcount=visitcount+1
	return

############################# End of Correct function ##################################


############################# Script begins here  ##################################

bank=''
if detbank == 'front-detector':
	bank = 'f'
if detbank == 'HAB':
	bank = 'h'
final_workspace = scatter_sample.split('_')[0] + bank + '_' + correction_type
# this calls the function above, set up lambda loop here !
# the final "sample" is only added as a suffix to the transmission workspace
Correct(scatter_sample, trans_sample, final_workspace, "sample")
#
# WHY does result log on screen give a lot of red "626 Object not found in ADS" ???  or 626_2 when doing 626_2_4 ???
#
# WHY is final value of transmission ZERO on qti plots but fortunately NOT zero in the show data ?
# issue with point mode data ?
# WHY do  _2_4 etc transmisions in nnn_trans_part NOT have the same values as nnn_trans_sample ?? -
#  presume it is doing a straight line fit, which over a small range will not give same answer as over a wide range ???  
print "STOP now, but you can execute more code to loop over wavelengths etc see below"
STOP
#
# BETTER to calc the transmission once only over the full wavelength range, then rebin for the others.
usePreviousTrans = trans_sample.split('_')[0] + "_trans_" + "sample"
print "usePreviousTrans = "+str(usePreviousTrans)
wav1=2.0
wav2=4.0
final_workspace = scatter_sample.split('_')[0] + bank + "_2_4"
Correct(scatter_sample, trans_sample, final_workspace, "part")
#
wav1=4.0
wav2=6.0
final_workspace = scatter_sample.split('_')[0] + bank + "_4_6"
Correct(scatter_sample, trans_sample, final_workspace, "part")
wav1=6.0
wav2=8.0
final_workspace = scatter_sample.split('_')[0] + bank + "_6_8"
Correct(scatter_sample, trans_sample, final_workspace, "part")
wav1=8.0
wav2=10.0
final_workspace = scatter_sample.split('_')[0] + bank + "_8_10"
Correct(scatter_sample, trans_sample, final_workspace, "part")
wav1=10.0
wav2=12.0
final_workspace = scatter_sample.split('_')[0] + bank + "_10_12"
Correct(scatter_sample, trans_sample, final_workspace, "part")
wav1=12.0
wav2=14.0
final_workspace = scatter_sample.split('_')[0] + bank + "_12_14"
Correct(scatter_sample, trans_sample, final_workspace, "part")
#
# re-set the defaults for next time through
wav1 = 2.2
wav2 = 10.
usePreviousTrans = ' '


wlist=["_1D","_2_4","_4_6","_6_8","_8_10","_10_12","_12_14"]
run=scatter_sample.split('_')[0] + bank
i=0
#  WHY do first two plosts come out in Horizontal steps not Lines ?
#  NOTE neither horizontal (block to right of data) nor vertical steps (block to left of data) is a proper histogram which would have the "point" in the centre of the block.
#  ALSO when qtiplot uses symbols or dots on graphs it seems to put the top left corner of the symbol at the data point, not the centre of the symbol !
for wext in wlist:
	wsp=run+wext
	if wext==wlist[0]:
		plot1=plotSpectrum(wsp,0)
		layer=plot1.activeLayer()
		layer.setTitle(run)
	else:
		mergePlots(plot1,plotSpectrum(wsp,0))
	layer.setCurveTitle(i,wsp)
	i=i+1
# HOW do you stop elegantly ?
STOP
# go round again with the other detector bank
MaskStringRear = 'h0,h190>h191,v0,v191,v64>v127+h64>h191,v128>v191+h0>h127'
MaskStringRear = 'h0,h190>h191,v0,v191'
MaskStringRear = 'h0,h190>h191,v0,v191,h131>h140'
detbank = "rear-detector"
detbank = "front-detector"
correction_type = '2D'
bank='_mask'
#
if detbank == 'front-detector':
	specmin = 192*192 + 1 + monstart
	specmax = 192*192*2 + monstart
if detbank == 'rear-detector':
	specmin = 1 + monstart
	specmax = 192*192 +monstart
bank=''
if detbank == 'front-detector':
	bank = 'f'
if detbank == 'HAB':
	bank = 'h'
usePreviousTrans = " "
wav1=2.0
wav2=14.0
final_workspace = scatter_sample.split('_')[0] + bank + '_' + correction_type
Correct(scatter_sample, trans_sample, final_workspace, "sample")
print final_workspace
# NEED then to loop over the CAN run, do the subtractions and write everything out.	
rmin=0.038
print rmin

# NEED v2 of SaveRKH to write out a series of 1d workspaces, with proper titles, into the same file (i.e append)
#
if scatter_can:
        Correct(scatter_can, trans_can, "tmp_can_output", "can")
        Minus(final_workspace, "tmp_can_output", final_workspace)
        mantid.deleteWorkspace("tmp_can_output")
print "end"

import SANSadd
SANSadd.add_runs(path, [774,780])
wksp = '769_sans_raw'
CropWorkspace(wksp, OutputWorkspace = 'rear', StartWorkspaceIndex= 8, EndWorkspaceIndex=(192*192 + 7))
CropWorkspace(wksp, OutputWorkspace = 'front', StartWorkspaceIndex= (192*192+8), EndWorkspaceIndex=(192*192*2 + 7))
# help(SumRowColumn)
# SumRowColumn(InputWorkspace, OutputWorkspace, Orientation, XMin=-9999, XMax=-9999, HVMin=-9999, HVMax=-9999)
# Orientation is "D_H" or "D_V" at present, Xmin/max is for TIMES, HVmin/max for rows or cols 
SumRowColumn('rear', 'rear_H', 'D_H', XMin=-9999, XMax=-9999, HVMin=-9999, HVMax=-9999)
SumRowColumn('rear', 'rear_V', 'D_V', XMin=-9999, XMax=-9999, HVMin=-9999, HVMax=-9999)
SumRowColumn('front', 'front_H', 'D_H', XMin=-9999, XMax=-9999, HVMin=-9999, HVMax=-9999)
SumRowColumn('front', 'front_V', 'D_V', XMin=-9999, XMax=-9999, HVMin=-9999, HVMax=-9999)
for i in range (52,62):
	SumRowColumn('rear', 'rear_H'+str(i), 'D_H', HVMin=i, HVMax=i)
SumRowColumn('rear', 'rear_V78', 'D_V', HVMin=78, HVMax=78)
SumRowColumn('front', 'front_H61', 'D_H', HVMin=61, HVMax=61)
#
# run the can through the whole process, rear & front dets, then subtract below
sam = '914'
can = '915'
Minus(sam+'_1D', can+'_1D', sam+'_sc_1d')
Minus(sam+'f_1D', can+'f_1D', sam+'f_sc_1d')
#SaveRKH(InputWorkspace="...",Filename="...",Append="0")
SaveRKH(sam+'_sc_1d',sam+'.Q',Append='0')
SaveRKH(sam+'f_sc_1d',sam+'.Q',Append='1')
#
#  SC wav for REAR only
sam = '914'
can = '915'
wlist=["_1D","_2_4","_4_6","_6_8","_8_10","_10_12","_12_14"]
i=0
for wext in wlist:
	wsp=sam+'-'+can+wext
	Minus(sam+wext, can+wext, wsp)
	if wext==wlist[0]:
		plot1=plotSpectrum(wsp,0)
		layer=plot1.activeLayer()
		layer.setTitle(wsp)
		SaveRKH(wsp,sam+'.Q',Append='0')
	else:
		mergePlots(plot1,plotSpectrum(wsp,0))
		SaveRKH(wsp,sam+'.Q',Append='1')
	layer.setCurveTitle(i,wsp)
	i=i+1
#

speclist = SANSUtility_RKH.ConvertToSpecList('v2+h0>h3', specmin, dimension)
print speclist
speclist = SANSUtility_RKH.ConvertToSpecList('h1', specmin, dimension)
print speclist
#SANSUtility_RKH.MaskBySpecNumber(resultWS, speclist)
#
wksp="SANS2D00000653"
CropWorkspace(wksp, OutputWorkspace = "D_T", StartWorkspaceIndex=(8+1-1), EndWorkspaceIndex=str(8+192*192-1))
# sum all into a single spectrum  this does the COLETTE, DISPLAY TIME or D T
SumSpectra("D_T","D_T")
Integration("D_T","D_T_integral",10,100000)	
CropWorkspace(wksp, OutputWorkspace = "D_Tf", StartWorkspaceIndex=(8+192*192-1), EndWorkspaceIndex=str(8+192*192*2-1))
# sum all into a single spectrum  this does the COLETTE, DISPLAY TIME or D T
SumSpectra("D_Tf","D_Tf")
Integration("D_Tf","D_Tf_integral",10,100000)	

import SANSUtility_RKH
# D_T for nine blocks
#wksp="SANS2D00000859"
wksp="SANS2D00000899"
spec1=8+1
dimdet=192
tag="low"
#print speclist
i=0
namelist=["A TL","D TC","G TR","B ML","E MC","H MR","C BL","F BC","I BR"]
masklist=["h0>h127,v64>v191","h0>h127,v0>v63,v128>v191","h0>h127,v0>v127","h0>h63,h128>h191,v64>v191","h0>h63,h128>h191,v0>v63,v128>v191","h0>h63,h128>h191,v0>v127","h63>h191,v64>v191","h63>h191,v0>v63,v128>v191","h63>h191,v0>v127"]
for name in namelist:
	print "i= "+str(i)
	list = SANSUtility_RKH.ConvertToSpecList(masklist[i], spec1, dimdet)
	CropWorkspace(wksp, OutputWorkspace = "D_T "+name+tag, StartWorkspaceIndex=(spec1-1), EndWorkspaceIndex=str(spec1+192*192-2))
	SANSUtility_RKH.MaskBySpecNumber("D_T "+name+tag, list)
	SumSpectra("D_T "+name+tag,"D_T "+name+tag)
	i=i+1
print "done"
for wext in namelist:
	wsp="D_T "+wext
	if wext==namelist[0]:
		plot1=plotSpectrum(wsp+tag,0)
		layer=plot1.activeLayer()
		layer.setTitle("D_T modules")
	else:
		mergePlots(plot1,plotSpectrum(wsp+tag,0))
	layer.setCurveTitle(i,name)
#



