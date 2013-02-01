######################################################################
#
# This is a partial copy from LeBailFitScript.py
#
# Python Script as Step 1 of Le Bail Fitting to
# 1. Load file
# 2. Create LeBailFitInput
# 3. Fit Peaks
#
# Step 1:   Load data, model (HKL list) and starting instrument parameters values,
#           and do an initial LeBailFit/calculation to see how much the starting
#           values are off;
#
######################################################################
from Calibration_ImportInformation import *

global numsteps 

numsteps = 1000

#--------------  Definition of Global Variables ---------------------
bankid = 0

datafilename = "" 
hklfilename = ""
irffilename = ""

datawsname = ""
instrparamwsname = ""
braggpeakparamwsname = ""

# Range for Le Bail Fit of all peaks 
startx = -1
endx =  -1

# Range for fitting single peaks for step 1~3 
tofmin_singlepeaks = -1
tofmax_singlepeaks = -1

# Background related
backgroundtype = "Polynomial"
backgroundorder = 6
bkgdtablewsname = ""
bkgdwsname = ""
bkgdfilename = ""
usrbkgdpoints = ''

latticesize = 4.1568899999999998
#--------------------------------------------------------------------

def setupGlobals(infofilename):
    """ Set up globals values
    """

    global datafilename, hklfilename, irffilename  
    global datawsname, instrparamwsname, braggpeakparamwsname
    global bkgdtablewsname, bkgdwsname, bkgdfilename, backgroundorder
    global tofmin_singlepeaks, tofmax_singlepeaks, startx, endx
    global bankid, latticesize
    global usrbkgdpoints

    bankid, calibDict = importCalibrationInformation(infofilename)
    bankid = int(bankid)

    datafilename = calibDict["DataFileDir"] + calibDict[bankid]["DataFileName"] 
    hklfilename  = calibDict["HKLFileDir"]  + calibDict[bankid]["HKLFileName"]
    irffilename  = calibDict["IrfFileDir"]  + calibDict[bankid]["IrfFileName"]
    
    startx = float(calibDict[bankid]["LeBailFitMinTOF"])
    endx   = float(calibDict[bankid]["LeBailFitMaxTOF"])

    # Name of workspaces
    datawsname = calibDict[bankid]["DataWorkspace"]
    instrparamwsname     = "Bank%sInstrumentParameterTable" % (bankid)
    braggpeakparamwsname = 'BraggPeakParameterTable'

    # Background related
    usrbkgdpoints   = calibDict[bankid]["UserSpecifiedBkgdPts"]
    bkgdwsname      = datawsname+"_Background"
    backgroundtype  = calibDict["BackgroundType"]
    backgroundorder = int(calibDict["BackgroundOrder"])
    bkgdfilename    = calibDict["WorkingDir"] + datawsname + "_Parameters.bak'"
    bkgdwsname      = datawsname + "_Background"
    bkgdtablewsname = datawsname + "_Background_Parameters"

    # Other constants
    latticesize   = calibDict[bankid]["LatticeSize"]

    return


def doStep4(numsteps):
    """ Use Monte Carlo random/drunken walk for solution
            FitRegion                   = '8500,49000',
            BackgroundParametersWorkspace   ='PG3_10808_Background_Parameters',
    """ 
    global datawsname 
    global instrparamwsname
    global braggpeakparamwsname 
    global startx, endx
    global backgroundorder, bkgdtablewsname

    outputwsname = '%s_MC_%s' % (datawsname, numsteps)
    
    print instrparamwsname    
    print braggpeakparamwsname
    print startx, endx

    LeBailFit(
            InputWorkspace                  = datawsname,
            OutputWorkspace                 = outputwsname,
            InputParameterWorkspace         = 'Bank%dInstrumentParameterTable1_0' % (bankid),
            OutputParameterWorkspace        = 'Bank%dInstrumentParameterTable_MC' % (bankid),
            InputHKLWorkspace               = 'BraggPeakParameterTable1',
            OutputPeaksWorkspace            = 'BraggPeakParameterTable3',
            FitRegion                       = '%f, %f' % (startx, endx),
            Function                        = 'MonteCarlo', 
            NumberMinimizeSteps             = numsteps, 
            BackgroundFunctionOrder         = backgroundorder,
            BackgroundParametersWorkspace   = bkgdtablewsname,
            UseInputPeakHeights             = False, 
            PeakRadius                      ='8',
            Minimizer                       = 'Levenberg-Marquardt',
            Damping                         = '0.90000000000000004',
            FitGeometryParameter            = True, 
            RandomSeed                      = numsteps, 
            AnnealingTemperature            = 10.0,
    	    DrunkenWalk                     = True)

    return


def main(argv):
    """ Main
    """
    global numsteps
    if numsteps <= 0:
	    numsteps = 200
	    print "Using default number of steps 200"

    setupGlobals("/home/wzz/Projects/MantidTests/LeBailFit/FinalExam/Calibration_Information.txt")

    doStep4(numsteps)

    return

if __name__=="__main__":
    main([])
