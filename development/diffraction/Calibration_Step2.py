################################################################################
#
# Iteratively fit single peaks of diffraction pattern, then use the 
# the fitted peak parameters to fit powder instrument parameters. 
#
# It will be migrated to a Mantid algorithm if appropriate
#
################################################################################
from Calibration_ImportInformation import *

class FitPowderInstrumentParameters:
    """ Class to implement the algorithm
    """
    def __init__(self, inputdict=None):
        """ Initialization
        """
        self._inputdict = inputdict

        return

    def init(self):
        """ Define input properties
        """

        return


    def getProperty(self, paramname):
        """ Get property. THIS IS A FAKE METHOD
        """
        return (self._inputdict[paramname])


    def execute(self):
        """ Main execution body
        """
        # 1. Get inputs
        datawsname = self.getProperty("InputWorkspace")
        instwsname = self.getProperty("InstrumentParameterWorkspace")
        peakwsname = self.getProperty("BraggPeakParameterWorkspace")
        outdatawsname = self.getProperty("OutputWorkspace")
        tofmin_singlepeaks = self.getProperty("TofMin")
        tofmax_singlepeaks = self.getProperty("TofMax")
	self._bankid = self.getProperty("BankID")
	self._monteCarloFile = self.getProperty("MCParmetersFile")
        fitpeakcentreonly = self.getProperty("FitPeakCentreOnly")
	
	self._instwsname = instwsname
	self._peakwsname = peakwsname

        # 2. Fit powder diffraction peaks 
        FitPowderDiffPeaks(InputWorkspace=datawsname, 
                OutputWorkspace=outdatawsname,
                BraggPeakParameterWorkspace=peakwsname, 
                InstrumentParameterWorkspace=instwsname, 
                OutputBraggPeakParameterWorkspace=peakwsname+"_0",
                MinTOF=tofmin_singlepeaks, 
                MaxTOF=tofmax_singlepeaks)
	""" FIXME: Extend the functionality
                ConfidentInInput=confident, 
                NumPeaksToFit=numpeaks)
	"""

        tablews = mtd[instwsname]
        self._sig0 = getValueFromTable(tablews, "Sig0")
        self._sig1 = getValueFromTable(tablews, "Sig1")
        self._sig2 = getValueFromTable(tablews, "Sig2") 

        self._alph0 = getValueFromTable(tablews, "Alph0")
        self._alph1 = getValueFromTable(tablews, "Alph1")
        self._alph0t = getValueFromTable(tablews, "Alph0t")
        self._alph1t = getValueFromTable(tablews, "Alph1t")

        self._beta0 = getValueFromTable(tablews, "Beta0")
        self._beta1 = getValueFromTable(tablews, "Beta1")
        self._beta0t = getValueFromTable(tablews, "Beta0t")
        self._beta1t = getValueFromTable(tablews, "Beta1t")

        # 4. Fit Parameters in loop
        continueloop = True
	iteration = 1
        maxloop = 3
        while continueloop is True:
            self.doFit(iteration, fitpeakcentreonly)
	    iteration += 1
            if iteration > maxloop:
                break
        # ENDWHILE

        # 5. Apply fitted values to 
        CloneWorkspace(InputWorkspace=instwsname, OutputWorkspace=instwsname+"_01")
	self.applyFitResult(instwsname+"_01")

        return 
    
    def doFit(self, iteration, fitpeakcentreonly=False):
        """ Do fit once 
        """ 
        # 1. Transform to workspace2D for future fitting : previous iteration's output
        tablewsname = self._peakwsname+"_"+str(iteration-1)
        if fitpeakcentreonly: 
            # Fit peak centre only! 
            raise NotImplementedError("Not mature!")
            transformTableWorkspace(tablewsname, ["PeakCentres"], "d_h", ["TOF_h"])
        else: 
            transformTableWorkspace(tablewsname, ["PeakCentres", "PeakAlphas", "PeakBetas", "PeakSigmas"], 
                "d_h", ["TOF_h", "Alpha", "Beta", "Sigma"])

        # 2. Fit Sigma 
        sig0 = self._sig0
        sig1 = self._sig1
        sig2 = self._sig2
        sigmafunction = 'name=ThermalNeutronBk2BkExpSigma,Sig0=%f,Sig1=%f,Sig2=%f' % (sig0, sig1, sig2)
        
        # a) Evaluate for sigma about the screening 
        Fit(Function=sigmafunction, InputWorkspace='PeakSigmas', 
                MaxIterations='0', Minimizer='Levenberg-Marquardt', 
                CreateOutput='1', Output='PeakSigmaEvalulation')

        # b) Fit
        Fit(Function=sigmafunction, InputWorkspace='PeakSigmas', 
                MaxIterations='1000', Minimizer='Levenberg-Marquardt', 
                CreateOutput='1', Output='PeakSigmaLM')

        # c) Zscore
        CalculateZscore(InputWorkspace="PeakSigmaLM_Workspace", 
                OutputWorkspace="PeakSigmaLM_Zscore",
                WorkspaceIndex=2)

        # d) Remove peaks with 'bad' fit 
        self.removeBadPeaks(tablewsname, "PeakSigmaLM_Zscore", 0, "d_h", 4.0, tolerance=0.0001)

        # e) Transform again
        transformTableWorkspace(tablewsname, ["PeakCentres", "PeakAlphas", "PeakBetas", "PeakSigmas"], 
                "d_h", ["TOF_h", "Alpha", "Beta", "Sigma"])

        # 3. Fit geometry parameters
        instrumenttablews = mtd["Bank%sInstrumentParameterTable1"%(self._bankid)]

        # a) Import the parameter boundary file
        paramdict = importParameterBoundaryFile(self._monteCarloFile)

        # b) Update the instrumental parameter workspace
        updateLeBailTable(instrumenttablews, paramdict) 

        # c) Refine parameters
        RefinePowderInstrumentParameters(BraggPeakParameterWorkspace=tablewsname,
            InstrumentParameterWorkspace='Bank%dInstrumentParameterTable1'%(self._bankid),
            OutputWorkspace='Bank%dPeakPositions'%(self._bankid),
            OutputInstrumentParameterWorkspace='Bank%dInstrumentParameterTable2'%(self._bankid),
            OutputBestResultsWorkspace='BestMCResult1', 
            RefinementAlgorithm = "MonteCarlo",
            RandomWalkSteps='1000',
            ParametersToFit='Dtt1,Dtt1t,Dtt2t,Zerot,Width',
            NumberBestFitRecorded='10',
            MonteCarloRandomSeed='0',
	    StandardError='InvertedPeakHeight')

        # d) Zscore
        fitws = mtd["Bank%dPeakPositions"%(self._bankid)]
        CalculateZscore(InputWorkspace=fitws, 
                OutputWorkspace="%s_Zscore" % (str(fitws)),
                WorkspaceIndex=2)

        # e) Remove peaks with 'bad' fit 
        self.removeBadPeaks(tablewsname, "%s_Zscore"%(str(fitws)), 0, "d_h", 3.0, tolerance=0.0001)

        # f) Transform again
        transformTableWorkspace(tablewsname, ["PeakCentres", "PeakAlphas", "PeakBetas", "PeakSigmas"], 
                "d_h", ["TOF_h", "Alpha", "Beta", "Sigma"])

        tablews = mtd["BestMCResult1"]
        # FIXME Use the BEST MC result as default
        # TODO  Shall make it more flexible to use n-th BEST Monte Carlo result
        width = getValueFromBestMCResultTable(tablews, "Width", 0)
        tcross = getValueFromBestMCResultTable(tablews, "Tcross", 0)

        # 4. Fit Alpha
        alph0  = self._alph0 
        alph1  = self._alph1 
        alph0t = self._alph0t
        alph1t = self._alph1t
        alphafunction = 'name=ThermalNeutronBk2BkExpAlpha,Width=%f,Tcross=%f,Alph0=%f,Alph1=%f,Alph0t=%f,Alph1t=%f,ties=(Width=%f,Tcross=%f,Alph1t=%f)' % (width, 
                tcross, alph0, alph1, alph0t, alph1t, width, tcross, alph1t)

        # a) Alpha evalulation 
        Fit(Function=alphafunction,InputWorkspace='PeakAlphas',MaxIterations='0', 
                CreateOutput='1',Output='PeakAlphaEvaluation')

        # b) Alpha fit 
        Fit(Function=alphafunction,InputWorkspace='PeakAlphas', 
                MaxIterations='1000',CreateOutput='1',Output='PeakAlphaLM')

        # c) Zscore
        CalculateZscore(InputWorkspace="PeakAlphaLM_Workspace", 
                OutputWorkspace="PeakAlphaLM_Zscore",
                WorkspaceIndex=2)

        # d) Remove peaks with 'bad' fit 
        self.removeBadPeaks(tablewsname, "PeakAlphaLM_Zscore", 0, "d_h", 4.0, tolerance=0.0001)

        # e) Transform again
        transformTableWorkspace(tablewsname, ["PeakCentres", "PeakAlphas", "PeakBetas", "PeakSigmas"], 
                "d_h", ["TOF_h", "Alpha", "Beta", "Sigma"])

        # 5. Fit Beta
        beta0  = self._beta0 
        beta1  = self._beta1 
        beta0t = self._beta0t
        beta1t = self._beta1t
        betafunction = 'name=ThermalNeutronBk2BkExpBeta,Width=%f,Tcross=%f,Beta0=%f,Beta1=%f,Beta0t=%f,Beta1t=%f,ties=(Width=%f,Tcross=%f)' % (width, 
                tcross, beta0, beta1, beta0t, beta1t, width, tcross) 

        # a) Evalulate Beta 
        Fit(Function=betafunction,InputWorkspace='PeakBetas', 
                MaxIterations='0',CreateOutput='1',Output='PeakBetaEvaluation')

        # b) Fit Beta 
        Fit(Function=betafunction,InputWorkspace='PeakBetas', 
                MaxIterations='1000',CreateOutput='1',Output='PeakBetaLM')

        # c) Zscore
        CalculateZscore(InputWorkspace="PeakBetaLM_Workspace", 
                OutputWorkspace="PeakBetaLM_Zscore",
                WorkspaceIndex=2)

        # d) Remove peaks with 'bad' fit 
        self.removeBadPeaks(tablewsname, "PeakBetaLM_Zscore", 0, "d_h", 4.0, tolerance=0.0001)

        # e) Transform again
        transformTableWorkspace(tablewsname, ["PeakCentres", "PeakAlphas", "PeakBetas", "PeakSigmas"], 
                "d_h", ["TOF_h", "Alpha", "Beta", "Sigma"])

        # 6. Create a new table workspace
        tablewsname = self._peakwsname+"_"+str(iteration-1)
        newtablewsname = self._peakwsname+"_"+str(iteration)
        SelectPowderDiffPeaks(BraggPeakParameterWorkspace=tablewsname, 
                ZscoreWorkspace='ZscoreTable', 
                OutputBraggPeakParameterWorkspace=newtablewsname, 
                MinimumPeakHeight=0.01, 
                ZscoreFilter='Alpha, 10000000.0')

        return

    def applyFitResult(self, instrumentwsname):
	""" Apply the fit results to instrument parameter workspace
	It will go throw Peak_Parameter, Alpha_Parameter, Beta_Parameter workspaces to 
	apply the value there to instrument workspace
	"""
        # 1. Best MC results (for peak position fit)
        # FIXME As prototype, workspace name and index are fixed.  They should be flexible.
        paramws = mtd["BestMCResult1"]
        paramnames = paramws.getColumnNames()
        for i in xrange(len(paramnames)):
            parname = paramnames[i]
            parvalue = paramws.cell(0, i) 
            UpdatePeakParameterTableValue(InputWorkspace=instrumentwsname, Column="Value", 
                    ParameterNames=[parname], NewFloatValue=parvalue)

        # 2. For Sigma, Alpha and Beta
	paramwsnames = ["PeakAlphaLM_Parameters", "PeakBetaLM_Parameters",
	    "PeakSigmaLM_Parameters"]

	for paramwsname in paramwsnames:
	    paramws = mtd[paramwsname]
	    numrows = paramws.rowCount()
	    for i in xrange(numrows):
		parname = paramws.cell(i, 0)
		parvalue = paramws.cell(i, 1)
		UpdatePeakParameterTableValue(InputWorkspace=instrumentwsname, Column="Value", 
                        ParameterNames=[parname], NewFloatValue=parvalue)
	    #ENDFOR
	#ENDFOR

	return

    def removeBadPeaks(self, peakwsname, zscorewsname, wsindex, parname, maxzscore, tolerance):
        """ Remove peaks (HKL) of the parameters with bad zscore

        Argument:
         - peakwsname:   name of the workspace containing peak parameters
         - zscorewsname: name of the workspace containing zscores
         - wsindex     : workspace index in Zscore workspace as a reference filter
         - parname     : parameter name in peakworkspace for zscore's X value
         - maxzscore   : maximum z-score allowed
        """
        # 1. Collect the data points exceeding Zscore
        vecx = mtd[zscorewsname].readX(wsindex)
        vecy = mtd[zscorewsname].readY(wsindex)

        dataptstokill = []
        zscoretokill = []
	print "Select Peak From %s According To Zscore Workspace : %s" % (peakwsname, zscorewsname)
        for i in xrange(len(vecy)):
            #print "Z(%f) = %f" % (vecx[i], vecy[i])
            if vecy[i] > maxzscore:		
                dataptstokill.append(vecx[i])
                zscoretokill.append(vecy[i])
                #print "Data point %f has Zscore %f exeeding the maximum limit %f." % (
                #        vecx[i], vecy[i], maxzscore)
            # ENDIF
        # ENDFOR

        # 2. Remove peak
        rows = []
        paramdicts = parseTableWorkspace(peakwsname)
        for i in xrange(len(paramdicts)):
            paramdict = paramdicts[i]
            for j in xrange(len(dataptstokill)):
	        datapt = dataptstokill[j]
                if abs(paramdict[parname]-datapt) < tolerance:
                    rows.append(i)
                    print "Peak (%d, %d, %d) @ %s = %f as %d-th entry, having zscore = %f, will be removed." % (
                            paramdict["H"], paramdict["K"], paramdict["L"], parname,
                            paramdict[parname], i, zscoretokill[j])
                # ENDIF
            # ENDFOR
        # ENDFOR

        DeleteTableRows(TableWorkspace=peakwsname, Rows=rows)

        return

    #ENDDEFCLASS

#------------------------------------------------------------------------------
# TableWorkspace Processing
#------------------------------------------------------------------------------

def getValueFromBestMCResultTable(tablews, parname, index):
    """ Get the number [index] best MC result from MC table tablews for parname
    """
    #  1. Read table workspace
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


def parseTableWorkspace(tablewsname):
    """ Parse a table workspace to a list of dictionaries
    """
    dicts = []
    tablews = mtd[str(tablewsname)]

    colnames = tablews.getColumnNames()
    numcols = len(colnames)
    numrows = tablews.rowCount()

    for i in xrange(numrows):
        paramdict = {}

        for j in xrange(numcols):
            parname = colnames[j]
            parvalue = tablews.cell(i, j)
            paramdict[parname] = parvalue
        # ENDFOR

        dicts.append(paramdict)
    # ENDFOR

    return dicts

def transformTableWorkspace(tablewsname, outputwsnames, xaxisname, yaxisnames):
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
        print "Error! Input X-axis name %s is not found in the table workspace columns" % (xaxisname)
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
#------------------------------------------------------------------------------
# ASCII File Processing
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
# Main... ...
#------------------------------------------------------------------------------

def main(argv):
    """ Main
    "MCParmetersFile": "/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August_Bank1/Bank1InstrumentMC.dat"
    """
    infofilename = "/home/wzz/Projects/MantidTests/LeBailFit/FinalExam/Calibration_Information.txt"
    bankid, calibDict = importCalibrationInformation(infofilename)
    bankid = int(bankid)

    inputdict = { 
            "InputWorkspace": calibDict[bankid]["DataWorkspace"],
	    "InstrumentParameterWorkspace": "Bank%dInstrumentParameterTable1"%(bankid),
            "BraggPeakParameterWorkspace": "BraggPeakParameterTable2",
            "OutputWorkspace": "Bank%dFittedPeaks" % (bankid),
            "TofMin": calibDict[bankid]["SinglePeakFitMin"], 
            "TofMax": calibDict[bankid]["SinglePeakFitMax"], 
            "BankID": bankid,
            "MCParmetersFile": calibDict["WorkingDir"] + "Bank%dInstrumentMC.dat" % (bankid)
            }


    # 1. Set up input
    # if bankid == 1:
    #     inputdict = {
    #             "InputWorkspace": "PG3_10808",
    #             "InstrumentParameterWorkspace": "Bank1InstrumentParameterTable1",
    #             "BraggPeakParameterWorkspace": "BraggPeakParameterTable2",
    #             "OutputWorkspace": "Bank1FittedPeaks",
    #             "TofMin": 19516.0,
    #             "TofMax": 49225.0,
    #             "BankID": 1,
    #             "MCParmetersFile": "/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August_Bank1/Bank1InstrumentMC.dat"
    #             }
    # elif bankid == 2:
    #     inputdict = {
    #             "InputWorkspace": "PG3_10809",
    #             "InstrumentParameterWorkspace": "Bank2InstrumentParameterTable1",
    #             "BraggPeakParameterWorkspace": "BraggPeakParameterTable2",
    #             "OutputWorkspace": "Bank2FittedPeaks",
    #             "TofMin": 17000.0,
    #             "TofMax": 74000.0,
    #             "BankID": 2,
    #             "MCParmetersFile": "/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August_Bank2/Bank2InstrumentMC.dat"
    #             }

    # elif bankid == 3:
    #     inputdict = {
    #             "InputWorkspace": "PG3_10810",
    #             "InstrumentParameterWorkspace": "Bank3InstrumentParameterTable1",
    #             "BraggPeakParameterWorkspace": "BraggPeakParameterTable2",
    #             "OutputWorkspace": "Bank2FittedPeaks",
    #             "TofMin": 17000.0,
    #             "TofMax": 74000.0,
    #             "BankID": 2,
    #             "MCParmetersFile": "/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August_Bank3/Bank3InstrumentMC.dat"
    #             }

    inputdict["FitPeakCentreOnly"] = False

    fitalg = FitPowderInstrumentParameters(inputdict)

    fitalg.execute()

if __name__=="__main__":
    main([])
