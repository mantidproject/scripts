######################################################################
#
# This is the formal final script!
#
# Python Script as Step 1 of Le Bail Fitting to
# 1. Load file
# 2. Create LeBailFitInput
# 3. Fit Peaks
#
# Step 1: Load data, model (HKL list) and starting instrument parameters values,
#           and do an initial LeBailFit/calculation to see how much the starting
#           values are off;
# 
# Step 2: Fit single peaks for TOF_h, Alpha, Beta, and Sigma;
# 
# Step 3: Plot parameters (TOF_H, Alpha, Beta and Sigma) against d-spacing;
# 
# Step 3.5:  Remove the peaks with bad fit;
# 
# Step 4: Refine instrument geometry parameters; 
#  
#         It is possible to loop back to Step 2 to include more peaks with single peaks
#
# Step 5: Do Le Bail Fit from previous result to see whether the peak parameters are 
#         close enough for Le Bail Fit
#
# Step 6: Save the result files to HDD for Step 6, LeBailFit in random walk;
#
# Step 7: Not in this script;
#
# Step 8: Process Monte Carlo results from Step 6
#
#
######################################################################

#--------------  BankID and Step To Execute (User Control) ----------
usrbankid = 2
usrstep = 3

#--------------  Definition of Global Variables ---------------------
bankid = 0

datafilename = "" 
hklfilename = ""
irffilename = ""
montecarlofilename = ""

expirffilename = ""

datawsname = ""
instrparamwsname = ""
braggpeakparamwsname = ""

outdataws1name = ""

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


def setupGlobals(usrbankid):
    """ Set up globals values
    """
    global datafilename, hklfilename, irffilename  
    global datawsname, instrparamwsname, braggpeakparamwsname
    global outdataws1name , montecarlofilename
    global bkgdtablewsname, bkgdwsname, bkgdfilename
    global expirffilename
    global tofmin_singlepeaks, tofmax_singlepeaks, startx, endx
    global bankid, latticesize
    global usrbkgdpoints

    bankid = usrbankid

    if bankid == 1:
	# Bank 1
	datafilename = '/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August/FERNS-LaB6/PG3_10808-1.dat'
	hklfilename = "/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August/Reflections/LB4844b1.hkl"
	irffilename = r'/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August/PeakProfiles/2011B_HR60b1.irf'

	expirffilename = r'/home/wzz/Projects/MantidTests/LeBailFit/PeakPositionFitting/UnitTest_NewBank1/Bank1.irf'

	montecarlofilename = '/home/wzz/Projects/MantidTests/LeBailFit/PeakPositionFitting/UnitTest_NewBank1/Bank1InstrumentMC.dat'

	datawsname = "PG3_10808"
    
	outdataws1name = "PG3_10808_FittedSinglePeaks"

        startx = 5053.       
        endx =  49387.

        tofmin_singlepeaks = 15000.
        tofmax_singlepeaks = 50000.

        bkgdtablewsname = "PG3_10808_Background_Parameters"
        bkgdwsname = "PG3_10808_Background"
        bkgdfilename = '/home/wzz/Projects/MantidTests/LeBailFit/PeakPositionFitting/UnitTest_NewBank1/PG3_10808_Background_Parameters.bak'

        latticesize = 4.1568899999999998

        usrbkgdpoints = '5243,8910,11165,12153,13731,15060,16511,17767,19650,21874,23167,24519,36000,44282,49000'   
    
    elif bankid == 2:
	# Bank2
	datafilename = '/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August/FERNS-LaB6/PG3_10809-2.dat'
	hklfilename = "/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August/Reflections/LB4853b2.hkl"
	irffilename = r'/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August/PeakProfiles/2011B_HR60b2.irf'

	expirffilename = "/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August_Bank2/Bank2.irf"

	montecarlofilename = '/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August_Bank2/Bank2InstrumentMC.dat'

	datawsname = "PG3_10809"
    
	outdataws1name = "PG3_10809_FittedSinglePeaks"

        startx = 7000.0
        endx =  70700.0

        tofmin_singlepeaks = 17322.
        tofmax_singlepeaks = 70700.


        bkgdtablewsname = "PG3_10809_Background_Parameters"
        bkgdwsname = "PG3_10808_Background"
        bkgdfilename = '/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August_Bani2/PG3_10808_Background_Parameters.bak'

        latticesize = 4.1568899999999998

        usrbkgdpoints = '5243,8910,11165,12153,13731,15060,16511,17767,19650,21874,23167,24519,36000,44282,49000'   

    else:
	# Bank ?
	raise NotImplementedError("To be implemented ASAP for Bank %d!" % (bankid))


    instrparamwsname = "Bank%sInstrumentParameterTable" % (bankid)
    braggpeakparamwsname = 'BraggPeakParameterTable'


    return

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

    
""" Input Information Required """
# ReflectionFile = ""
# FullprofParameterFile = ""
# 
# LatticeConstant = 1.0
# 
# MinTOF = 15000
# MaxTOF = 50000

def Step1():
    """ Step 1: Load data, parameters and do an initial/test Le Bail Fit/calculation
    """
    global datafilename, hklfilename, irffilename
    global datawsname, instrparamwsname, braggpeakparamwsname
    global outdataws1name
    global tofmin_singlepeaks, tofmax_singlepeaks
    global bankid
    global latticesize

    # 1. Load File 
    LoadAscii(Filename=datafilename, 
	OutputWorkspace=datawsname, 
	Unit='TOF')

    # 2. Create Input Tables
    CreateLeBailFitInput(ReflectionsFile=hklfilename, 
        FullprofParameterFile=irffilename, 
        LatticeConstant=latticesize, 
        InstrumentParameterWorkspace=instrparamwsname+"1", 
        BraggPeakParameterWorkspace=braggpeakparamwsname+"1",
        Bank=bankid)
	
    # 3. Use Le Bail Fit 's calculation feature to generate the starting peak parameters (for fitting) values (can be off)
    LeBailFit(InputWorkspace=datawsname,
            OutputWorkspace="CalculatedPattern1",
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

def Step2():
    """ Step 2: Fit single peaks
    """
    global datafilename, hklfilename, irffilename
    global datawsname, instrparamwsname, braggpeakparamwsname
    global outdataws1name
    global tofmin_singlepeaks, tofmax_singlepeaks
    global bankid

    # 1. Use LeBailFit/Calculation to create the initial peak parameters
    # print "Fit Range: %f to %f" % (tofmin_singlepeaks, tofmax_singlepeaks)
    # print "InputWorkspace = %s" % (datawsname)
    # print "BraggPeakParameterWorkspace = %s" % (braggpeakparamwsname+"2")
    # print "InstrumentParameterWorkspace = %s" % (instrparamwsname+"1")
    # print "OutputWorkspace = %s" % (outdataws1name)
    FitPowderDiffPeaks(InputWorkspace=datawsname,
        OutputWorkspace=outdataws1name,
        BraggPeakParameterWorkspace=braggpeakparamwsname+"2",
        InstrumentParameterWorkspace=instrparamwsname+"1",
        OutputBraggPeakParameterWorkspace=braggpeakparamwsname+"3",
        MinTOF=tofmin_singlepeaks,
	MaxTOF=tofmax_singlepeaks)

    print "\n............................... and you need to delete some peak by simple observation.\n"

    return
   
#def Step1_5():
#    """ Step 1.5: Best via observation! 
#    """ 
#    SelectPowderDiffPeaks(BraggPeakParameterWorkspace=braggpeakparamwsname+"3", 
#            ZscoreWorkspace='ZscoreTable',
#            OutputBraggPeakParameterWorkspace=braggpeakparamwsname+"4",
#            MinimumPeakHeight=minpeakheight, 
#            ZscoreFilter='Alpha, 3.0, Beta, 3.0') 


def Step3():
    """ Convert the resultant table workspace of peak parameters and convert several workspace2D
    """
    global braggpeakparamwsname

    tablewsname = braggpeakparamwsname+"3"
    TransformTableWorkspace(tablewsname, ["PeakCentres", "PeakAlphas", "PeakBetas", "PeakSigmas"], 
            "d_h", ["TOF_h", "Alpha", "Beta", "Sigma"])

    return

def Step4():
    """ Step 4: Refine the instrumental peak profile
    """ 
    # 1. Remove peaks with parameters away from model
    global minpeakheight,  montecarlofilename
    global bankid

    instrumenttablews = mtd["Bank%sInstrumentParameterTable1"%(bankid)]

    # 2. Import the parameter boundary file
    paramdict = importParameterBoundaryFile( montecarlofilename)

    # 3. Update the instrumental parameter workspace
    updateLeBailTable(instrumenttablews, paramdict) 

    # 4. Refine parameters
    RefinePowderInstrumentParameters(BraggPeakParameterWorkspace=braggpeakparamwsname+"4",
	InstrumentParameterWorkspace='Bank%dInstrumentParameterTable1'%(bankid),
        OutputWorkspace='Bank%dPeakPositions'%(bankid),
        OutputInstrumentParameterWorkspace='Bank%dInstrumentParameterTable2'%(bankid),
        OutputBestResultsWorkspace='BestMCResult1', 
        RefinementAlgorithm = "MonteCarlo",
        RandomWalkSteps='1000',
        ParametersToFit='Dtt1,Dtt1t,Dtt2t,Zerot,Width',
        NumberBestFitRecorded='10',
        MonteCarloRandomSeed='0')

    print "... ... Consider Loop Back to Step 2 To Increase Range of Single Peaks"

    return


def Step5():
    """ Step 5: to check whether the refinement result is good for low-d region
    FunctionName : name=Chebyshev,EndX=1,StartX=-1,n=6,A0=0.657699,A1=3.68433e-05,
                   A2=-6.48189e-09,A3=2.91945e-13,A4=-5.89298e-18,A5=5.53119e-23,A6=-1.95804e-28'
    """
    global startx, endx
    global datawsname, instrparamwsname, braggpeakparamwsname
    global bkgdtablewsname, bkgdwsname, backgroundtype, backgroundorder
    global bankid
    global usrbkgdpoints

    # 1. Process background 
    bkgdwsname = datawsname+"_Background"
    ProcessBackground(InputWorkspace=datawsname, 
            OutputWorkspace=bkgdwsname, 
            Options='SelectBackgroundPoints',
            LowerBound=startx, 
            UpperBound=endx, 
            BackgroundType=backgroundtype,
            BackgroundPoints=usrbkgdpoints,
            NoiseTolerance='0.10000000000000001')

    functionstr = "name=%s,n=%d" % (backgroundtype, backgroundorder)
    for iborder in xrange(backgroundorder+1):
        functionstr = "%s,A%d=%.5f" % (functionstr, iborder, 0.0)
    print "Background function: %s" % (functionstr)

    Fit(Function='name=Polynomial,n=6,A0=0.657699,A1=3.68433e-05,A2=-1.29638e-08,A3=1.16778e-12,A4=-4.71439e-17,A5=8.8499e-22,A6=-6.26573e-27',
            InputWorkspace=bkgdwsname,
            MaxIterations='1000',
            Minimizer='Levenberg-MarquardtMD',
            CreateOutput='1',
            Output=bkgdwsname,
            StartX=startx,
            EndX=endx)

    # 2. Do calculation
    paramdict = parseRefineGeometryMCResultWorkspace(mtd["BestMCResult1"]) 

    index = 0
    for chi2 in sorted(paramdict.keys()):
        # 1. Update parameter values
        updateInstrumentParameterValue(mtd["Bank%sInstrumentParameterTable1"%(bankid)], paramdict[chi2])

        LeBailFit(InputWorkspace=datawsname,
            OutputWorkspace='CalculatedPattern%d'%(index),
            InputParameterWorkspace='Bank%dInstrumentParameterTable1'%(bankid),
            OutputParameterWorkspace='Bank%dInstrumentParameterTable1_%d'%(bankid, index),
            InputHKLWorkspace='BraggPeakParameterTable1',
            OutputPeaksWorkspace='BraggPeakParameterTable2_%d'%(index),
            Function='Calculation',
            BackgroundType='Polynomial',
            BackgroundFunctionOrder='6',
            #BackgroundParameters='0.657699, 3.68433e-05, -1.29e-08, 1.168e-12, -4.714e-17, 8.8499e-22, -6.266e-27', 
            BackgroundParametersWorkspace=bkgdtablewsname, 
            UseInputPeakHeights='0',
            PeakRadius='8',
            Minimizer='Levenberg-Marquardt')

        print "Pattern %s:  Instrumental Geometry Chi^2 = %.5f" % ('CalculatedPattern%d'%(index), chi2)
        index += 1
    # ENDFOR

    return


def Step4():
    """ Step 4: Save result to disk for LeBailFit/Random Walk launched in parallel
    Launch a multiple-threaded MC random walk for the best Le Bail Fit result
    """
    # TODO  How to convert the information such that a Monte Carlo can be performed ... ...
    global bkgdfilename
    global bkgdtablewsname, backgroundtype, startx, endx
    global expirffilename

    # 1. Background
    bkgdtablews = mtd[bkgdtablewsname]
    exportBackground(bkgdfilename, bkgdtablews, backgroundtype, startx, endx)

    # 2. Save .irf file
    paramdict = parseRefineGeometryMCResultWorkspace(mtd["BestMCResult1"]) 

    index = 0
    for chi2 in sorted(paramdict.keys()):
        # 1. Update parameter values
        updateInstrumentParameterValue(mtd["Bank1InstrumentParameterTable1"], paramdict[chi2])
	irffilenameX = expirffilename.split(".irf")[0]+"_Stage1_"+str(index)+".irf"
        print "Export MC Result To File %s" % (irffilenameX)
        SaveFullprofResolution(InputWorkspace="Bank1InstrumentParameterTable1",
                OutputFile=irffilenameX,
                Bank=1)
	index += 1
    # ENDFOR

    return

def Step5():
    """ Step 5: Test Load background
    """
    global bkgdfilename
    
    results = importBackground(bkgdfilename)
    tablews = results[0]
    bkgdtype = results[1]
    startx = results[2]
    endx = results[3]
    print "Import Background Parameter Workspace %s.  Type = %s.  Range = %f to %f" % (tablews.name(), bkgdtype, startx, endx)

    return


def Step6():
    """ Step 6:  Process LeBailFit/Random Walk Result
    """
    print "Step6: Process LeBailFit/random walk result. "
    print "       Call LeBailScript_Step6ProcessMCResults.py"
    return


def main(argv):
    """ Main
    """
    global datafilename, datawsname

    if len(argv) < 3: 
        print "%s [Bank ID]  [Step]" % (argv[0])
	return

    bankid = int(argv[1]) 
    step = int(argv[2])

    setupGlobals(bankid)

    if step == 1: 
	print datafilename
        Step1()

    elif step == 2:
        Step2()

    elif step == 3:
        Step3()

    elif step == 4:
        Step4()

    elif step == 5:
        Step5()

    elif step == 6:
        Step6()

    else:
        print "Step %s is not defined." % (step)

    return


if __name__=="__main__":
    global minpeakheight 
    global usrbankid, usrstep
   
    """ Information Required For Step 2 """
    minpeakheight = 0.31  
   
    
    main(["LeBailFitScript", usrbankid, usrstep])
