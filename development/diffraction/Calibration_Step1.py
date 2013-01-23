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

#|-------------------------------------------------------------------------------------------
#|Step 2.1: Fit single peaks for TOF_h, Alpha, Beta, and Sigma;
#|
#|Step 2.2: Plot parameters (TOF_H, Alpha, Beta and Sigma) against d-spacing;
#|
#|Step 2.3: Remove the peaks with bad fit;
#|
#|Step 2.4: Refine instrument geometry parameters; 
#| 
#|          It is possible to loop back to Step 2 to include more peaks with single peaks
#|-------------------------------------------------------------------------------------------
#
# Step 5: Do Le Bail Fit from previous result to see whether the peak parameters are 
#         close enough for Le Bail Fit
#
# !Step 6: Fit Alpha, Beta and apply the change to 
#
# !Step 7: Save the result files to HDD for Step 6, LeBailFit in random walk;
# 
# Step 8: Save the result files to HDD for Step 6, LeBailFit in random walk;
#
# Step 7: Not in this script;
#
# Step 8: Process Monte Carlo results from Step 6
#
#
######################################################################
from Calibration_ImportInformation import *

#--------------  Definition of Global Variables ---------------------
bankid = 0

datafilename = "" 
hklfilename = ""
irffilename = ""

# montecarlofilename = ""
# expirffilename = ""

datawsname = ""
instrparamwsname = ""
braggpeakparamwsname = ""

# outdataws1name = ""

minpeakheight = 0.001

# Range for Le Bail Fit of all peaks 
startx = -1
endx =  -1
# Range for fitting single peaks for step 1~3 
tofmin_singlepeaks = -1
tofmax_singlepeaks = -1

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
    global bkgdtablewsname, bkgdwsname, bkgdfilename
    global tofmin_singlepeaks, tofmax_singlepeaks, startx, endx
    global bankid, latticesize
    global usrbkgdpoints

    bankid, calibDict = importCalibrationInformation(infofilename)
    bankid = int(bankid)

    datafilename = calibDict["DataFileDir"] + calibDict[bankid]["DataFileName"] 
    hklfilename  = calibDict["HKLFileDir"]  + calibDict[bankid]["HKLFileName"]
    irffilename  = calibDict["IrfFileDir"]  + calibDict[bankid]["IrfFileName"]
    
    datawsname = calibDict[bankid]["DataWorkspace"]
    
    startx = float(calibDict[bankid]["LeBailFitMinTOF"])
    endx   = float(calibDict[bankid]["LeBailFitMaxTOF"])


    # expirffilename = r'/home/wzz/Projects/MantidTests/LeBailFit/PeakPositionFitting/UnitTest_NewBank1/Bank1.irf'
    # montecarlofilename = '/home/wzz/Projects/MantidTests/LeBailFit/PeakPositionFitting/UnitTest_NewBank1/Bank1InstrumentMC.dat'
    # outdataws1name = "PG3_10808_FittedSinglePeaks"
    
    
    # tofmin_singlepeaks = 15000.
    # tofmax_singlepeaks = 50000.
    
    # bkgdtablewsname = "PG3_10808_Background_Parameters"
    # bkgdwsname = "PG3_10808_Background"
    # bkgdfilename = '/home/wzz/Projects/MantidTests/LeBailFit/PeakPositionFitting/UnitTest_NewBank1/PG3_10808_Background_Parameters.bak'
    
    latticesize   = calibDict[bankid]["LatticeSize"]
    usrbkgdpoints = calibDict[bankid]["UserSpecifiedBkgdPts"]

    instrparamwsname     = "Bank%sInstrumentParameterTable" % (bankid)
    braggpeakparamwsname = 'BraggPeakParameterTable'

    return

#------------------------------------------------------------------------------
# File Processing
#------------------------------------------------------------------------------

def importParameterBoundaryFile(paramfilename):
    """ Import Instrumental Geometry Parameter File Range
    """
    try:
        infile = open(paramfilename, "r")
    except IOError:
	 print "Unable to open file %s" % (paramfilename)
	 raise IOError("Unable to open parameter boundary file %s" % (paramfilename))
    lines = infile.readlines()
    infile.close()

    # Parse
    paramdict = {}
    for line in lines:
        line = line.strip()
        if len(line) == 0:
            continue
        elif line[0] == '#':
            continue
        else:
            terms = line.split()
            name = terms[0]
            value = float(terms[1])
            parmin = float(terms[2])
            parmax = float(terms[3])
            stepsize = float(terms[4])
            
            paramdict[name] = [value, parmin, parmax, stepsize]
        # ENDIF
    # ENDFOR

    return paramdict


#------------------------------------------------------------------------------
# TableWorkspace Processing
#------------------------------------------------------------------------------

def getValueFromBestMCResultTable(tablews, parname, index):
    """ Get the number [index] best MC result from MC table tablews for parname
    """
    # 1. Read table workspace
    colnames = tablews.getColumnNames()
    colindex = -1
    for i in xrange(len(colnames)):
        if colnames[i] == parname:
            colindex = i
            break
    # ENDFOR

    # 2. Return if parameter is not found
    if colindex < 0:
        raise NotImplementedError("Parameter (columne) %s is not found in TableWorkspace %s" %
                (parname, str(tablews)))

    if index < 0 or index >= tablews.rowCount():
        raise NotImplementedError("MC result table contains %d best results.  Index %d is out of range." % 
                (tablews.rowCount(), index))

    # 3. Get value
    parvalue = tablews.cell(index, colindex)

    return parvalue

def getValueFromTable(tablews, parname):
    """ Get a value from a TableWorkspace such as parameter tableworkspace
    Requirement: Must have 2 columns named 'Name' and 'Value'
    """
    # 1. Get column names
    colnames = tablews.getColumnNames()
    colnamedict = {}
    for i in xrange(len(colnames)):
        colnamedict[colnames[i]] = i
    if colnamedict.has_key("Name") and colnamedict.has_key("Value"):
        pass
    else:
        raise NotImplementedError("TableWorkspace %s has either column 'Name' or 'Value' or both" % (str(tablews)))

    # 2. Find row
    numrows = tablews.rowCount()
    icolname = colnamedict["Name"]
    irow = -1
    for i in xrange(numrows):
        if tablews.cell(i, icolname) == parname:
            irow = i
            break
    # ENDFOR
    if irow < 0:
        raise NotImplementedError("TableWorkspace %s does not have parameter named %s" % 
                str(tablews), parname)

    # 3. Get value
    icol = colnamedict["Value"]
    parvalue = tablews.cell(irow, icol)

    return parvalue
   

def updateLeBailTable(instrumenttablews, paramdict):
    """ Update the instrumental geometry parameters MC to table workspace
    """ 
    wsname = instrumenttablews.name()

    for parname in paramdict.keys():
        parinfo =paramdict[parname]
        parval = parinfo[0]
        parmin = parinfo[1]
        parmax = parinfo[2]
        parstp = parinfo[3]

        UpdatePeakParameterTableValue(InputWorkspace=wsname,Column='Min',
                ParameterNames=parname,NewFloatValue=parmin)
        
        UpdatePeakParameterTableValue(InputWorkspace=wsname,Column='Max',
                ParameterNames=parname,NewFloatValue=parmax)
        
        UpdatePeakParameterTableValue(InputWorkspace=wsname,Column='StepSize',
                ParameterNames=parname,NewFloatValue=parstp)
    # ENDFOR 

    return


def parseRefineGeometryMCResultWorkspace(tablews):
    """ Parse the refinement result workspac
    """
    numrows = tablews.rowCount()
    numcols  = tablews.columnCount()
    colnames = tablews.getColumnNames()

    parameterdict = {}

    for irow in xrange(numrows):
        ppdict = {}
        for icol in xrange(numcols):
            colname = colnames[icol]
            value = tablews.cell(irow, icol)
            ppdict[colname] = value
        # ENDFOR Column
        chi2 = ppdict["Chi2"]

        parameterdict[chi2] = ppdict
    # ENDFOR

    return parameterdict

def updateInstrumentParameterValue(tablews, paramdict):
    """ Update the value to an instrument parameter table
    """
    paramnames = paramdict.keys()
    for parname in paramnames: 
        parvalue = paramdict[parname]
	# print "%s = %f" % (parname, parvalue)
	if parname.count("Chi2") == 0:
	    # Only update parameters nothing to do with chi2
            UpdatePeakParameterTableValue(InputWorkspace=tablews,
		Column='Value',
                ParameterNames=[parname],
		NewFloatValue=parvalue)

    return

def exportBackground(ofilename, bkgdtablews, backgroundtype, startx, endx):
    """ Export background background 
    """
    wbuf = ""
    wbuf += "Type,\t%s\n" % (backgroundtype)
    wbuf += "StartX,\t%f\n" % (startx)
    wbuf += "EndX,\t%f\n" % (endx)

    # Parse table
    colnames = bkgdtablews.getColumnNames()
    colnamedict = {}
    for icol in xrange(len(colnames)):
	colname = colnames[icol]
	colnamedict[colname] = icol
    iname = colnamedict["Name"]
    ivalue = colnamedict["Value"]

    numrows = bkgdtablews.rowCount()
    for irow in xrange(numrows):
	parname = bkgdtablews.cell(irow, iname)
	parvalue = bkgdtablews.cell(irow, ivalue)
	wbuf +="%s,\t%.10e\n" % (parname, parvalue)

    # Write
    ofile = open(ofilename, "w")
    ofile.write(wbuf)
    ofile.close()

    return
	
def importBackground(bkgdfilename):
    """ Import background file name.  

    Return: Parameter TableWorkspace, Background Type
    """
    # 1. Import
    ifile = open(bkgdfilename, "r")
    lines = ifile.readlines()
    ifile.close()

    # 2. Parse to dictionary
    pardict = {}
    for line in lines:
        line = line.strip()
        if len(line) == 0:
    	# Empty Line
    	    continue

        terms = line.split(",")
        if len(terms) != 2:
    	    print "Warning! Line [%s] is not expected." % (line)
    	    continue

        parname = terms[0].strip()
        parvaluestr = terms[1].strip()

        pardict[parname] = parvaluestr

    # ENDFOR

    # 3. Set up
    bkgdtype = pardict.pop("Type")
    startx = float(pardict.pop("StartX"))
    endx = float(pardict.pop("EndX"))

    # 4. Create table workspace
    tablews = CreateEmptyTableWorkspace(OutputWorkspace="LoadedBackgroundParameters")

    tablews.addColumn("str", "Name")
    tablews.addColumn("double", "Value")

    for parname in sorted(pardict.keys()):
        parvalue = float(pardict[parname])
        tablews.addRow([parname, parvalue])

    return (tablews, bkgdtype, startx, endx)


def TransformTableWorkspace(tablewsname, outputwsnames, xaxisname, yaxisnames):
    """ Transform the content in a table worskpace with specified format
    to several Workspace2D.
    """
    # 1. Validate input arguments
    tablews = mtd[str(tablewsname)]
    if tablews is None:
        print "Error!  Unable to locate TableWorkspace %s" % (str(tablewsname))
        return

    if len(outputwsnames) != len(yaxisnames):
        print "Error!  Input output workspace names and y-axis names are not in pair"
        return

    # 2. Information for table workspace
    numrows =  tablews.rowCount()
    colnames = tablews.getColumnNames()
    colnamedict = {}
    for i in xrange(len(colnames)):
        colname = colnames[i]
        colnamedict[colname] = i

    # 3. Create workspace2Ds
    # a) X-axis
    if colnamedict.has_key(xaxisname) is False:
        print "Error! Input X-axis name %s is not found in the table workspace columns" % (xaxisnamea)
        return
    else:
        ixcol = colnamedict[xaxisname]
    # ENDIF
    xarray = []
    for irow in xrange(numrows):
        x = tablews.cell(irow, ixcol)
        xarray.append(x)

    for iy in xrange(len(yaxisnames)):
        # b) Build arrays
        yaxisname = yaxisnames[iy]
        outputwsname = outputwsnames[iy]

        if colnamedict.has_key(yaxisname) is False:
            print "Error: Input Y-axis name %s is not found in the table workspace columns" % (yaxisnamea)
            continue

        iycol = colnamedict[yaxisname]

        yarray = []
        earray = []
        for irow in xrange(numrows):
            y = tablews.cell(irow, iycol)
            yarray.append(y)
	    earray.append(1.0)
        # ENDFOR Table Row

        # c) Create output workspace
        CreateWorkspace(OutputWorkspace=outputwsname, DataX=xarray, DataY=yarray, DataE=earray, 
                NSpec=1, UnitX="dSpacing")
    # ENDFOR Y-axis name

    return

    

def doStep1(): 
    """ Step 1: Load data, parameters and do an initial/test Le Bail Fit/calculation 
    """
    global datafilename, hklfilename, irffilename
    global datawsname, instrparamwsname, braggpeakparamwsname
    global outdataws1name
    global tofmin_singlepeaks, tofmax_singlepeaks
    global bankid
    global latticesize
    
    # 1. Load File 
    print datawsname
    LoadAscii(Filename=datafilename, OutputWorkspace=datawsname, Unit='TOF')
    
    # 2. Create Input Tables
    CreateLeBailFitInput(ReflectionsFile=hklfilename, 
        FullprofParameterFile=irffilename, 
        LatticeConstant=float(latticesize), 
        InstrumentParameterWorkspace=instrparamwsname+"1", 
        BraggPeakParameterWorkspace=braggpeakparamwsname+"1",
        Bank=bankid)
        
    # 3. Use Le Bail Fit 's calculation feature to generate the starting peak parameters (for fitting) values (can be off)
    LeBailFit(InputWorkspace=datawsname,
            OutputWorkspace="PatternWildGuess",
            InputParameterWorkspace=instrparamwsname+"1",
            OutputParameterWorkspace='TempTable',
            InputHKLWorkspace=braggpeakparamwsname+"1",
            OutputPeaksWorkspace=braggpeakparamwsname+"2",
            Function='Calculation',
            BackgroundType='Chebyshev',
            BackgroundFunctionOrder='1',
            BackgroundParameters='0,0,0',
            UseInputPeakHeights='0',
            PeakRadius='8',
            Minimizer='Levenberg-Marquardt')
    
    return


def main(argv):
    """ Main
    """
    setupGlobals("/home/wzz/Projects/MantidTests/LeBailFit/FinalExam/Calibration_Information.txt")
    
    print "Le Bail Fit Calibration Instrument... Step 1"
    doStep1()

    return


if __name__=="__main__":
    """ starter
    """
    main(["LeBailFitScript"])
