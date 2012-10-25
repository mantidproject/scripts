################################################################################
#
#  Process the results/output files from Le Bail Fit/Random Walk
#  This is a script based on MantidPlot
#  It serves as Step6 in LeBailFit Script
#
################################################################################

"""----------------------   Definition of Global Variables  ----------------------------"""
bankid = 1
directory = "./Result02/"
seedrange = 400

datadict = {
        1: {"run": 10808},
        2: {"run": 10809},
        }

"""-------------------------------------------------------------------------------------"""

_PARAMETERSPOOL = ["Zero", "Dtt1", "Zerot", "Dtt1t", "Dtt2t", "Width", "Tcross",
        "LatticeConstant", 
        "Alph0", "Alph1", "Beta0", "Beta1",
        "Alph0t", "Alph1t", "Beta0t", "Beta1t",
        "Sig0", "Sig1", "Sig2",
        "Gam0", "Gam1", "Gam2"]
	
class Parameter:
    """ A class to hold the information of a parameter 
    """
    def __init__(self, name):
        """ Initialization
        """
        self.name = name
        self.value = 0.0
        self.newvalue = 0.0
        self.valuechange = False
        self.minvalue = -0.0
        self.maxvalue = -0.0
        self.stepsize = 0.0

	return

    def __str__(self):
        """ Formatted output
        """
        rs = "%s = %.5E.  Range: [%.5E, %.5E].  MC Step size = %.3E.  Proposed new value ? = %s: %.5E" % (self.name,
                self.value, self.minvalue, self.maxvalue, self.stepsize, str(self.valuechange), self.newvalue)

        return rs


class LeBailCalculation:
    """ Class to do LeBail fitting with Monte Carlo
    """
    def __init__(self, bankid, runnumber):
        """ Initialization
        """
        # Local variables for bank 1
        datafilename = "PG3_%d-%d.dat" % (runnumber, bankid)
        
        if bankid == 1:
            hklfilename = "LB4844b1.hkl"
            self.xmin = 14346.0
            self.xmax = 48547.0

        elif bankid == 2:
	    hklfilename = "LB4853b2.hkl"
	    self.xmin = 20000.0
            self.xmax = 64000.0

        irffilename = "2011B_HR60b%d.irf" % (bankid)
        irffilename = "2012_Ashfia_B%d.irf" % (bankid)

        self.bankid = bankid
       
        # Parameters constant to all banks 
        self.backgroundtype = "Chebyshev"
        self.backgroundorder = 7
        self.lattice = 4.15689  # A

        directory = "/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August/Bank1Test/"
        datadir = "FERNS-LaB6/"
        hkldir = "Reflections/"
        irfdir = "PeakProfiles/" 

        # Background
        backgroundfilename = "./PG3_%d_Background_Parameters.bak" % (runnumber)
        results = self.loadBackgroundFile(backgroundfilename)
        self.bkgdtablews = results[0]
        self.bkgdparameterwsname = self.bkgdtablews.name()
        self.backgroundtype = results[1]
        bkgdxmin = results[2]
        bkgdxmax = results[3]

        if self.backgroundtype == "Chebyshev" and (abs(bkgdxmin-self.xmin) > 1.0E-5
                or abs(bkgdxmax-self.xmax) > 1.0E-5):
            print "For Chebyshev background, Xmin/Xmax for the background function must be same as LeBailFit's range."
            raise NotImplementedError("Chebyshev background range incompatible to LeBailFit range.")

       
        # (Composed) user variables
        self.datafilename = "%s%s%s" % (directory, datadir, datafilename)
        self.hklfilename = "%s%s%s" % (directory, hkldir, hklfilename)
        self.irffilename = "%s%s%s" % (directory, irfdir, irffilename)
			
        # Random walk parameter file
	self.parameterMCfilename = "%s%s" % (directory, "randomwalk.txt")
        
        self.datawsname = "Data%d_Bank%d" %(runnumber, bankid)
        self.bkgdwsname = "Data%d_Bank%d_Bkgd" %(runnumber, bankid)
        self.bkgdparameterrootname = "Background%d"%(runnumber)
        
        self.bkgdfunction = "name=%s,n=%d" % (self.backgroundtype, self.backgroundorder)
        self.peakparamwsname = "LeBailParameterBank%d" % (bankid)
        self.hklwsname = "HKLBank%d" % (bankid)
        self.calwsname = "Calculated_Bank%d" % (bankid)
        self.peakswsname = "PeaksParameters_Bank%d" % (bankid)

        self.wsindex = 0

        # ----  Instance variables for parameters ----

        # parameternames and paramNameIndexDict are of a pair for look up
        self.parameternames = []
        self.paramNameIndexDict = {}

        self.parameterdict = {}     # current parameter value
        self.tempparameterdict = {}  # temp storage for calculated/fitted result

        index = 0
        for parname in _PARAMETERSPOOL:
            self.parameternames.append(parname)
            self.paramNameIndexDict[parname] = index

            self.parameterdict[parname] = Parameter(parname)  # Trace current parameters
            self.tempparameterdict[parname] = Parameter(parname)

            index += 1
        # ENDFOR	
		
        # Names of parameters
        self.parametersToFit = []

	self.curchi2 = -1.0

	# ---- MC records ----
	self.recordChi2 = []  # Current Chi2
	self.recordCase = []  # Case (cal/fit)
	self.recordTake = []  # Decision

        self.bestchi2 = 1.0E100
        self.bestchi2step = -1

        self.bestresults = {} # key = chi2, value = [parameterdict, step]
        self.bestchi2s = []
        self.numresults = 5

        self.numvisitstore = 0

        return


    def initLeBailFit(self):
        """ Initialize LeBail Fitting, including
        (1) loading data;
        (2) processing background (initial)
        (3) generating inputs for LeBailFit
        """
        # 1. Load data
        LoadAscii(Filename=self.datafilename, OutputWorkspace=self.datawsname,Unit='TOF')
        CropWorkspace(InputWorkspace=self.datawsname, OutputWorkspace=self.datawsname, XMin = self.xmin, XMax = self.xmax)
        
        dataws = mtd[self.datawsname]
        
        # 2. Calculate background roughly 
        # ProcessBackground(InputWorkspace=self.datawsname,OutputWorkspace=self.bkgdwsname,Options='SelectBackgroundPoints', \
        # 	BackgroundPoints=self.bkgdpoints,NoiseTolerance='0.05')
        # 	
        # # 3. Fit background
        # startx = dataws.dataX(self.wsindex)[0]
        # endx = dataws.dataX(self.wsindex)[-1]
        # Fit(Function=self.bkgdfunction, InputWorkspace=self.bkgdwsname,\
	# 	Minimizer='Levenberg-MarquardtMD',CreateOutput='1',Output=self.bkgdparameterrootname, 
        #         StartX=startx,EndX=endx)
        
        # 3. Create LeBailFit inputs
        CreateLeBailFitInput(ReflectionsFile=self.hklfilename, 
                FullprofParameterFile=self.irffilename,
        	Bank=self.bankid,
                LatticeConstant=self.lattice,
                InstrumentParameterWorkspace=self.peakparamwsname,
                BraggPeakParameterWorkspace=self.hklwsname)

        # 4. Set up Monte Carlo random walk parameters from file
        self.setupMCParameters(self.parameterMCfilename)

        # 5. Import (Le Bail) parameters to storage
        self.parseParameterWorkspace(self.peakparamwsname, self.parameterdict)

        return


    def reproduceLeBail(self, parameterdict):
        """ Do a single step of Le Bail Fit, including
        (1) calcualtion
        (2) fit

        Return
         - better? :  boolean to show new chi2 is better
         - case    :  integer to show (1) fit is better or (2) cal is better
         - wsname  :  name of the peak parameter worskpace to pass the value to next step
        """
        # 0. Check workspace availability 
        for wsname in [self.datawsname, self.peakparamwsname, self.hklwsname, self.bkgdparameterwsname]: 
            ws = mtd[wsname]
            if ws is None: 
                raise NotImplementedError("In doLeBailFit(), Workspace %s does not exist" % (wsname))
	#ENDFOR

	# 1. Update the parameters proposed in MC
        print "******************************************"
        for parindex in xrange( len(self.parameternames) ):
            parname = self.parameternames[parindex]
            if parameterdict.has_key(parname):
                newvalue = parameterdict[parname]
		print "*\tUpdate %s to %f" % (parname, newvalue)
                self.updatePeakParameterTableValue(InputWorkspace=self.peakparamwsname, OutputWorkspace=self.peakparamwsname, 
                        Column='Value', ParameterNames=[parname], NewFloatValue=newvalue)
            # ENDIF
        # ENDFOR
        print "******************************************"

        # 2. Calculate LeBail Patterns
        # a) Fix all parameters and to calculate chi^2
	self.updatePeakParameterTableValue(InputWorkspace=self.peakparamwsname, OutputWorkspace=self.peakparamwsname,
                Column='FitOrTie', NewStringValue='t')

        # b) Output workspace name
	calwsname = "Calculated_Bank%d" % (self.bankid)
	newcalparameterswsname = "CalculatedParameters_Bank%d" % (self.bankid)

        # c) Remove existing output workspace
        try:
            ws = mtd[calwsname]
            DeleteWorkspace(Workspace=ws)
        except KeyError:
            pass

        try:
            ws = mtd[newcalparameterswsname]
            DeleteWorkspace(Workspace=ws)
        except KeyError:
            pass


	# b) Calculate (fit by fixing all peak parameters)
        print "LeBailFit/Calculation: Input Parameter Workspace = %s; Output Parameter Workspace = %s;" % \
                (self.peakparamwsname, newcalparameterswsname)
        print "                       Output Data Workspace = %s;     Output Peaks Parameter Workspace = %s; " % \
	       (calwsname, self.peakswsname)

        # c) 
        try: 
            LeBailFit( 
                    InputWorkspace=mtd[self.datawsname],  
                    OutputWorkspace=calwsname,  
                    InputParameterWorkspace=mtd[self.peakparamwsname],  
                    OutputParameterWorkspace=newcalparameterswsname, 
                    InputHKLWorkspace=mtd[self.hklwsname], 
                    OutputPeaksWorkspace=self.peakswsname, 
		    BackgroundType=self.backgroundtype, 
                    BackgroundParametersWorkspace=mtd[self.bkgdparameterwsname], 
                    Function='LeBailFit',
                    UseInputPeakHeights='0',
                    PeakRadius='8')

            self.calwsname = calwsname

        except Exception:
            print "Catch an error from LeBailFit (calculating) and continue."
            self.newcalchi2 = 1.0E100
            pass

        return 


    def parseParameterWorkspace(self, parameterwsname, parameterdict):
	""" Parse parameter table workspace to a parameter dictionary
        Arguments
         - parameterwsname  : name of the table workspace
         - parameterdict    : dictionary (input/output)
	
	return: chi2 ... -1 if does not exist
	""" 
        # 1. Get workspace
        tablews = mtd[parameterwsname]
        if tablews is None:
            raise NotImplementedError("Input table workspace's name is incorrect!")

	# 2. Create a dictionary for flexible columns
        tmpcolnames = tablews.getColumnNames()
	colnames = []
        colnamedict = {} 
        for ic in xrange(len(tmpcolnames)):
            colnamedict[tmpcolnames[ic]] = ic
	    colnames.append(tmpcolnames[ic])
	   
        #    Check table workspace
        if colnames.count("Name") != 1:
            raise NotImplementedError("Input table workspace is not supported due to lacking of column Name.")
        if colnames.count("Value") != 1:
            raise NotImplementedError("Input table workspace is not supported due to lacking of column Value.")
        hasmin = colnamedict.has_key("Min")
        hasmax = colnamedict.has_key("Max")
        hasstep = colnamedict.has_key("StepSize")

        # 3. Parse
        chi2 = -1.0

        numrows = tablews.rowCount()
    
        for irow in xrange(numrows):
    	    parname = tablews.cell(irow, colnamedict["Name"])
            if parameterdict.has_key(parname):
                # Value
                parvalue = tablews.cell(irow, colnamedict["Value"]) 
                parameterdict[parname].value = parvalue

                # Min
                if hasmin is True:
                    minvalue = tablews.cell(irow, colnamedict["Min"])
                    parameterdict[parname].minvalue = minvalue
                
                # Max
                if hasmax is True:
                    maxvalue = tablews.cell(irow, colnamedict["Max"])
                    parameterdict[parname].maxvalue = maxvalue
                    
                # Step size
                if hasstep is True:
                    temp = tablews.cell(irow, colnamedict["StepSize"])
                    parameterdict[parname].stepsize = temp
            elif parname == "Chi2":
                chi2 = tablews.cell(irow, colnamedict["Value"])
		if self.curchi2 < 0.0:
			self.curchi2 = chi2
                    
        # ENDFOR

        return chi2

    def proposeNewValue(self, paramindex):
        """ use MC to generate a new value

        newvalue = currentvalue + r * stepsize.  r is random number (-0.5, 0.5)
	"""
	parname = self.parameternames[paramindex]
        parameter = self.parameterdict[parname]
        
        stepsize = parameter.stepsize

        randomvalue = random.random()-0.5

        newvalue = self.parameterdict[parname].value + randomvalue*stepsize

        # If outside of boundary, use periodic boundary
        if newvalue > parameter.maxvalue:
            newvalue = parameter.minvalue + (newvalue - parameter.maxvalue)
        elif newvalue < parameter.minvalue:
            newvalue = parameter.maxvalue + (newvalue - parameter.minvalue)

	return newvalue

    def setupMCParameters(self, parameterMCfilename):
        """ Read input column file (1) to update the Monte Carlo random walk parameters (2)
        to parameter table workspace

        (1) Format: name, min, max, stepsize 
        (2) Parameters include min, max, step size
        """
        # 1. Open file
        ifile = open(parameterMCfilename, "r")
        lines = ifile.readlines()
        ifile.close()

        # 2. Parse
        for line in lines:
            line = line.strip()
            if len(line) == 0 or line.startswith('#'):
                # a. Empty line or comment line
                continue

            terms = line.split()
            if len(terms) < 4:
                # b. Terms in line does not reach definition
                continue

            # c. Parse
            parname = terms[0]
            minvalue = float(terms[1])
            maxvalue = float(terms[2])
            stepsize = float(terms[3])
	
            # d. Update the workspace
            self.updatePeakParameterTableValue(InputWorkspace=self.peakparamwsname, OutputWorkspace=self.peakparamwsname, 
                    Column='Min', ParameterNames=[parname], NewFloatValue=minvalue)
            self.updatePeakParameterTableValue(InputWorkspace=self.peakparamwsname, OutputWorkspace=self.peakparamwsname, 
                    Column='Max', ParameterNames=[parname], NewFloatValue=maxvalue)
            self.updatePeakParameterTableValue(InputWorkspace=self.peakparamwsname, OutputWorkspace=self.peakparamwsname, 
                    Column='StepSize', ParameterNames=[parname], NewFloatValue=stepsize)

            # e. 

        # ENDFOR

        return
    #ENDDEF setupMCParameters(self, parameterMCfilename)
    
    def updatePeakParameterTableValue(self, InputWorkspace, Column, ParameterNames=None,
            NewFloatValue=0.0, NewStringValue="", OutputWorkspace=""):
        """ Update peak parameter table values... this will be a simplied verison of 
        the corresponding Pythong algorithm in order to share the 
        """

        # 0. Convert names
        tablews = mtd[str(InputWorkspace)]
        column = Column  # Column name
        if ParameterNames is None:
            ParameterNames = []

        # 1. Get necessary information from table workspace
        colnames = tablews.getColumnNames()
        colnamedict = {}
        for ic in xrange( len(colnames) ):
            colnamedict[colnames[ic]] = ic

        numrows = tablews.rowCount()
        parnamedict = {}
        parameternames = []
        for irow in xrange(numrows):
    	    parname = tablews.cell(irow, 0)
            parnamedict[parname] = irow
            parameternames.append(parname)
        # ENDFOR

        # 2. Parse for the cells (index) for new value
        paramindexlist = []
        for parname in ParameterNames:
            parindex = parnamedict[parname]
            paramindexlist.append(parindex)
        # ENDFOR
        
        if len(paramindexlist) == 0:
            # Default
            for ind in xrange(numrows):
                paramindexlist.append(ind)
        # ENDIF

        usestring = False
        usedouble = False
        if column in ["FitOrTie", "Name"]:
            usestring = True
        else:
            usedouble = True

        indexcolumn = colnamedict[column]

        # 3. Set new value
        for paramindex in paramindexlist:
            if usestring is True:
                tablews.setCell(paramindex, indexcolumn, NewStringValue)
            elif usedouble is True:
                tablews.setCell(paramindex, indexcolumn, NewFloatValue)
            # ENDIF
        # ENDFOR

        return


    def loadBackgroundFile(self, bkgdfilename):
	""" Import background file name.  

	Return: Parameter TableWorkspace, Background Type
	return (tablews, bkgdtype, startx, endx) 
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
        CreateEmptyTableWorkspace(OutputWorkspace="LoadedBackgroundParameters")
        tablews = mtd["LoadedBackgroundParameters"]

        tablews.addColumn("str", "Name")
        tablews.addColumn("double", "Value")

        for parname in sorted(pardict.keys()):
            parvalue = float(pardict[parname])
            tablews.addRow([parname, parvalue])

        return (tablews, bkgdtype, startx, endx)
  
#####################################################################################

def reproduceResults(rwparamdict, bankid, runnumber, numoutputs=1):
    """ Reproduce results
    """
    # FIXME Should be a proper way to do the rating for best N results
    seeds = sorted(rwparamdict.keys())
    seed = seeds[0]
    bestchi2 = rwparamdict[seed][0]["Chi2"]
    bestseed = seed

    for ic in xrange(1, len(seeds)):
        seed = seeds[ic]
        chi2 = rwparamdict[seed][0]["Chi2"]
        if chi2 < bestchi2:
            bestchi2 = chi2
            bestseed = seed

    print "Best Solution From Seed = %d.  Chi2 = %.5E" % (bestseed, bestchi2)

    bestseeds = []
    bestseeds.append(bestseed)

    # 2. Reproduce
    for seed in bestseeds: 
        
        alg = LeBailCalculation(bankid, runnumber)
        alg.initLeBailFit()

        alg.reproduceLeBail(rwparamdict[seed][0])


    return

#####################################################################################

def parseBestParameterFile(paramfname):
    """ Parse (best) parameter files
    """
    # 1. Return structure
    parameterdicts = []
    
    # 2. 
    ifile = open(paramfname, "r")
    lines = ifile.readlines()
    ifile.close()

    # 3. 
    for line in lines:
        line = line.strip()
        if len(line) == 0:
            continue

        terms = line.split()
        
        # construct output if necessary
        if len(parameterdicts) == 0:
            for i in xrange(1, len(terms)):
                parameterdicts.append({})
        
        # parse 
        parname = terms[0]
        for i in xrange(1, len(terms)):
            value = float(terms[i])
            parameterdicts[i-1][parname] = value
        # ENDFOR
    # ENDFOR

    return parameterdicts

def parseRandomWalkResult(infofiles, paramfiles):
    """ Parse the output files
    """
    # 0. Return data structure.  dictionary (seed) of lists (rank of chi2) of dicitonary (parameter/value)
    rwparamdict = {}

    # 1. Parse
    seeds = sorted(paramfiles.keys())
    for seed in seeds:
        paramfname = paramfiles[seed]
        parameterdicts = parseBestParameterFile(paramfname)
        rwparamdict[seed] = parameterdicts

    return rwparamdict


#####################################################################################

def main(argv):
    """ Main
    """
    global datadict

    if len(argv) < 4:
        print "Format %s [Bank ID]  [Seed Range]  [Directory]" % (argv[0])
        return
    
    # 1. Inputs
    bankid = int(argv[1])
    seedrange = int(argv[2])
    directory = argv[3]

    runnumber = datadict[bankid]["run"]

    seeds = range(seedrange)

    infofiles = {}
    paramfiles = {}

    for seed in seeds:
        infofname = "%soutput_seed%d.txt" % (directory, seed)
        paramfname = "%sbestparameters_seed%d.dat" % (directory, seed)
	print infofname, paramfname

        try:
            ifile = open(infofname, "r")
            ifile.close()
            infofiles[seed] = infofname
        except IOError:
            print "[Warning] Information file %s for random seed %d cannot be opened." % (infofname, seed)

        try:
            ifile = open(paramfname, "r")
            ifile.close()
            paramfiles[seed] = paramfname
        except IOError:
            print "[Warning] Best-parameter file %s for random seed %d cannot be opened." % (paramfname, seed)

    # ----- Parse Information and parameter files ----
    rwparamdict = parseRandomWalkResult(infofiles, paramfiles)

    # ----- Reproduce the result -----
    reproduceResults(rwparamdict, bankid, runnumber, numoutputs = 3)

    return


if __name__=="__main__":
    import sys

    global directory, seedrange, bankid
    
    main(["", bankid, seedrange, directory])
    
    
