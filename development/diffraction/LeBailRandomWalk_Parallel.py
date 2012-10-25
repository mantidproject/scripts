################################################################################
#  Using Random-Walk to find the local minima of Le Bail Fit
#  FOR PARALLEL (MULTIPLE THREADING) EXECUTION
#  Version 2.2
#     New Features
#     1. Add fixed multiple step refinement
#  Version 2.3
#     1. Use LocalChi2 other than Chi2 for best-chi2
#  Version 2.4 
#     1. Multiple step fitting.  Use MC to determine whether to accept the fitting result. 
#  Version 2.5
#     1. Script is converted to one step in LeBailFitScript
#     2. Instrumental geometry parameters and background are fitted as 
#     a. XXX.irf
#     b. YYY.bak
################################################################################
import random
import math
import os
import sys
import shutil
import time
sys.path.append("/home/wzz/Mantid/Code/debug/bin")

from MantidFramework import mtd
mtd.initialise()
from mantidsimple import *
from mantid.api import *

_PARAMETERSPOOL = ["Zero", "Dtt1", "Zerot", "Dtt1t", "Dtt2t", "Width", "Tcross",
        "LatticeConstant", 
        "Alph0", "Alph1", "Beta0", "Beta1",
        "Alph0t", "Alph1t", "Beta0t", "Beta1t",
        "Sig0", "Sig1", "Sig2",
        "Gam0", "Gam1", "Gam2"]

# _PARAMETERSFIT = {
#         0: ["Alph0", "Alph1", "Beta0", "Beta1"],
#         1: ["Alph0t", "Alph1t", "Beta0t", "Beta1t"],
#         2: ["Sig0", "Sig1", "Sig2"]
#         }

_PARAMETERSFIT = {
        0: ["Alph0", "Beta0",],
        1: ["Alph1", "Beta1"],
        2: ["Alph0", "Beta0", "Sig0"],
        3: ["Alph1", "Beta1", "Sig1"],
        4: ["Alph0t", "Beta0t", "Sig2"],
        5: ["Alph1t", "Beta1t"]
        }
_PARAMETERSFIT = {
        0: ["Alph0",    "Beta0"],
	1: ["Alph0t",  "Beta0t"],
	2: ["Alph1",    "Beta1"],
	3: ["Alph1t",   "Sig0"],
        4: ["Beta1t",   "Sig1"],
	5: ["Sig2"]
        }

_MAXNONUPDATE = 5
	
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
        self.chi2 = -1.0
        self.fit = False

	return

    def __str__(self):
        """ Formatted output
        """
        rs = "%s = %.5E.  Range: [%.5E, %.5E].  MC Step size = %.3E.  Proposed new value ? = %s: %.5E.  Fit? = %s.  Chi^2 = %.5E" % (self.name,
                self.value, self.minvalue, self.maxvalue, self.stepsize, str(self.valuechange), self.newvalue, self.fit, self.chi2)

        return rs


class LeBailMonteCarlo:
    """ Class to do LeBail fitting with Monte Carlo
    """
    def __init__(self, bankid, runnumber):
        """ Initialization
        (1) Set up the (hardcoded) parameters
        """
        self.bankid = bankid

        # Data File
        datafilename = "PG3_%d-%d.dat" % (runnumber, bankid)
        
        # Local variables for bank X
        if bankid == 1:
            hklfilename = "LB4844b1.hkl"
            self.bkgdpoints = '14392.0, 16850.0, 24894.0, 29079.0, 30657.9, 32873.0, 36704.0, 40309.0, 44553.0, 48460.0'
            self.xmin = 14346.0
            self.xmax = 48547.0

        elif bankid == 2:
	    hklfilename = "LB4853b2.hkl"
	    self.bkgdpoints = '20809., 23088., 24426., 28907., 34700., 51038., 60327., 64000.0 '
	    self.xmin = 20000.0
            self.xmax = 64000.0

        elif bankid == 3:
            raise NotImplementedError("Implement ASAP")

        elif bankid == 4:
            raise NotImplementedError("Implement ASAP")

        elif bankid == 5:
            raise NotImplementedError("Implement ASAP")

        elif bankid == 6:
            raise NotImplementedError("Implement ASAP")

        elif bankid == 7:
            raise NotImplementedError("Implement ASAP")

        # Instrument Files
        irffilename = "Bank%d_Stage1.irf" % (bankid)
        """
        irffilename = "2011B_HR60b%d.irf" % (bankid)
        irffilename = "2012_Ashfia_B%d.irf" % (bankid)
        """
      
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

        # Lattice parameter
        if bankid < 6: 
            self.lattice = 4.15689  # A
        else:
            raise NotImplementedError("Put the right value ASAP for bank = %d" % (bankid))

        # Parameters constant to all banks 
        directory = "/home/wzz/Projects/MantidTests/LeBailFit/PG3_2012August/"
        datadir = "FERNS-LaB6/"
        hkldir = "Reflections/"
       
        # irfdir = "PeakProfiles/" 
        # directory = "./"

        # (Composed) user variables
        self.datafilename = "%s%s%s" % (directory, datadir, datafilename)
        self.hklfilename = "%s%s%s" % (directory, hkldir, hklfilename)
        self.irffilename = "%s%s%s" % ("./", "", irffilename)
			
        # Random walk parameter file
	self.parameterMCfilename = "%s%s" % (directory, "randomwalk.txt")
        
        self.datawsname = "Data%d_Bank%d" %(runnumber, bankid)
        self.bkgdwsname = "Data%d_Bank%d_Bkgd" %(runnumber, bankid)
        self.bkgdparameterrootname = "Background%d"%(runnumber)
        
        self.peakparamwsname = "LeBailParameterBank%d" % (bankid)
        self.hklwsname = "HKLBank%d" % (bankid)
        self.calwsname = "Calculated_Bank%d" % (bankid)
        self.peakswsname = "PeaksParameters_Bank%d" % (bankid)
        # self.bkgdfunction = "name=%s,n=%d" % (self.backgroundtype, self.backgroundorder)

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

        # ---- MC Variables ----
        self.StepSizeFactor = 1.0

	# ---- MC records ----
	self.recordChi2 = []  # Current Chi2
	self.recordCase = []  # Case (cal/fit)
	self.recordTake = []  # Decision

	self.curchi2 = 1.0E10
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
        #ProcessBackground(InputWorkspace=self.datawsname,OutputWorkspace=self.bkgdwsname,Options='SelectBackgroundPoints', \
        #	BackgroundPoints=self.bkgdpoints,NoiseTolerance='0.05')
        	
        # 3. Fit background
        #startx = dataws.dataX(self.wsindex)[0]
        #endx = dataws.dataX(self.wsindex)[-1]
        #Fit(Function=self.bkgdfunction, InputWorkspace=self.bkgdwsname,\
	#	Minimizer='Levenberg-MarquardtMD',CreateOutput='1',Output=self.bkgdparameterrootname, 
        #        StartX=startx,EndX=endx)
        
        # 3. Create LeBailFit inputs
        CreateLeBailFitInput(ReflectionsFile=self.hklfilename, 
                FullprofParameterFile=self.irffilename,\
        	Bank=self.bankid,
                LatticeConstant=self.lattice,
                InstrumentParameterWorkspace=self.peakparamwsname,
                BraggPeakParameterWorkspace=self.hklwsname)

        # 4. Set up Monte Carlo random walk parameters from file
        self.setupMCParameters(self.parameterMCfilename)

        # 5. Import (Le Bail) parameters to storage
        self.parseParameterWorkspace(self.peakparamwsname, self.parameterdict)

        return


    def setParametersToFit(self, parnames):
        """ Set parameters to fit
        An algorithm to recognize parameter with blurred definition will be use such as
        (1) exactly same and case insensitive;
        (2) *partialname
        (3) partialname*
        (4) *partialname*
	(5) partialname?
        """

        for parnametofit in parnames:
            # for each input parameter name to fit

            parnametofit = parnametofit.lower()

            for parname in self.parameternames:
                parname_low = parname.lower()

                if parname_low == parnametofit:
                    # exactly match
                    self.parametersToFit.append(parname)

                elif parnametofit.startswith("*") and parnametofit.endswith("*"):
                    # *XXXX*
                    corestr = parnametofit.split("*")[1]
                    if parname_low.count(corestr) >= 1:
                        self.parametersToFit.append(parname)

                elif parnametofit.startswith("*"):
                    # *NNNNN
                    corestr = parnametofit.split("*")[1]
                    if parname_low.endswith(parnametofit):
                        self.parametersToFit.append(corestr)

                elif parnametofit.endswith("*"):
                    # NNNNN*
                    corestr = parnametofit.split("*")[0]
                    if parname_low.startswith(corestr):
                        self.parametersToFit.append(parname)

		elif parnametofit.endswith("?"):
		    # NNNNN?
		    corestr = parnametofit.split("?")[0]
		    if parname_low.startswith(corestr) and len(parname_low) == len(parnametofit)+1:
                        self.parametersToFit.append(parname)

                # ENDIFELSE
            # ENDFOR

        # ENDFOR

        return
	

    def doLeBailFit(self, fixedstep=False):
        """ Do a multiple-step Le Bail Fit, including
        (1) calcualtion
        (2) fit

        MC Parameters that will be set:
        1. newfitchi2
        2. newcalchi2

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
        print "SECTION 1:  Update Peak Parameters ****************************" 
        wbuf = ""
        numchange = 0
        for parindex in xrange( len(self.parameternames) ):
            parname = self.parameternames[parindex]
            if self.parameterdict[parname].valuechange is True:
                newvalue = self.parameterdict[parname].newvalue
                self.updatePeakParameterTableValue(InputWorkspace=self.peakparamwsname, OutputWorkspace=self.peakparamwsname, 
                        Column='Value', ParameterNames=[parname], NewFloatValue=newvalue)
		wbuf += "*\tUpdate %s to %.7f\n" % (parname, newvalue)
                numchange += 1
            # ENDIF
        # ENDFOR
        print wbuf
        print "******* %d Parameters Updated ******************************\n" % (numchange)
        time.sleep(0.1) 

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
            if ws is not None:
                DeleteWorkspace(Workspace=calwsname)
        except (KeyError, ValueError) as e:
            pass

        try:
            ws = mtd[newcalparameterswsname]
            if ws is not None: 
                DeleteWorkspace(Workspace=ws)
        except (KeyError, ValueError) as e:
            pass

        print "LeBailFit/Calculation: Input Parameter Workspace = %s; Output Parameter Workspace = %s;" % \
                (self.peakparamwsname, newcalparameterswsname)
        print "                       Output Data Workspace = %s;     Output Peaks Parameter Workspace = %s; " % \
	       (calwsname, self.peakswsname)

        # d) LeBail Fit/Calculation for Starting Chi2 (no parameter to fit = calculation)
        lbcalok = True
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
            print "[Error] Catch an error from LeBailFit (calculating).  Continue to fit. "
            lbcalok = False
            pass

        # f) Process outcome
        earlyreturn = False
        startchi2s = None
        if lbcalok is True:
            # Fit without error
            try:
	        physicalchi2, startchi2s = self.parseFitResult(newcalparameterswsname)
                newcalchi2 = startchi2s[1]
            except KeyError:
                physicalchi2 = False
                print "[Error] LeBailFit/Calcualtion error and does not generate output parameter TableWorkspace."
                pass

            if physicalchi2 is True: 
                # Good fit
                print "Chi^2 (Calculated) = %.4E" % (newcalchi2)
            else:
                # Diverged fit
                print "[Error] Chi^2 is not physical.  Proposed peak parameters are not refinable. Return!"
                earlyreturn = True
        else:
            # Fit with error
            earlyreturn = True
        # ENDIFELSE

        # g) Early return!
        if earlyreturn is True:
            # Early-error return
            better = False
            case = -1  # -1 for no way to select
            thisbestwsname = ""
            newcalchi = 1.0E10
            newfitchi = 1.0E10
            return (better, case, thisbestwsname, newcalchi2, newfitchi2)
        # ENDIF 

        # 3. Mutliple step Le Bail Fit
        chi2slist = [startchi2s]
        # a) Get number of steps for fitting 
        if fixedstep is True:
            numsteps = len(_PARAMETERSFIT)
        else:
            numsteps = 1

        lbfitresults = {}

        # b) Init: Input parameter workspace
        inputparamwsname = self.peakparamwsname

        # c) Fit in multiple steps
        allparstofit = [] 
        for istep in xrange(numsteps+1):
            # i.   Output workspace name
            fitwsname = "Bank%dFit_Step%d" % (self.bankid, istep) 
            newparameterwsname = "FittedParameters_Bank%d_%d" % (self.bankid, istep) 

            # ii.  Set specific set of parameters to fit 
            if istep < numsteps:
                # Get parameters from pre-setup
                if fixedstep is True:
                    parstofit = _PARAMETERSFIT[istep]
                else: 
                    parstofit = self.fitCandidates
                #      Add the parameters to list 'allparstofit'
                for parname in parstofit:
                    if allparstofit.count(parname) == 0:
                        allparstofit.append(parname)
            else:
                # Get parameters from previous fitted for the LAST 
                parstofit = allparstofit
            # ENDIFELSE

            # iii. Fix all parameters but set the chosed on to fit
            updategood1 = self.updatePeakParameterTableValue(InputWorkspace=inputparamwsname, Column='FitOrTie', 
                    ParameterNames=[], NewStringValue='t')            
            updategood2 = self.updatePeakParameterTableValue(InputWorkspace=inputparamwsname, Column='FitOrTie', 
                    ParameterNames=parstofit, NewStringValue='f')         
            #  update parameter dictionary
            for parname in self.parameterdict.keys():
                if parname in parstofit:
                    self.tempparameterdict[parname].fit = True
                    self.tempparameterdict[parname].chi2 = -1.0
                else: 
                    self.tempparameterdict[parname].fit = False 
                    self.tempparameterdict[parname].chi2 = 0.0
            # ENDFOR each parameter

            # iv.  (Possible) early error return
            if updategood1 is False or updategood2 is False:
                print "[Error] Failed to update peak parameter workspace %s @ Step %d.  Return!" % (inputparamwsname, istep)
                better = False
                case = -1  # -1 for no way to select
                thisbestwsname = ""
                return (better, case, thisbestwsname, 1.0E20, 1.0E20)

            print "SECTION 2 Fit Step %s ****************************************" % (istep)
            # print "LeBailFit/Fit: Updated Input Parameter Workspace = %s; Output Parameter Workspace = %s" % \
            #         (inputparamwsname, newparameterwsname)
            # print "               Output Data Workspace = %s;       Output Peaks Parameter Workspace = %s " % \
            #         (fitwsname, self.peakswsname)
            wbuf = "Fit Parameters: "
            for param in parstofit:
                wbuf += "%s\t\t" % (param)
            print wbuf
		    
            # v.   Delete pre-existing output workspace
            try:
                ws = mtd[newparameterwsname]
                if ws is not None: 
                    DeleteWorkspace(Workspace=ws)
            except (KeyError, RuntimeError) as e:
                pass

            try:
                ws = mtd[fitwsname]
                DeleteWorkspace(Workspace=ws)
            except (KeyError, RuntimeError) as e:
                pass
            
            # vi.   Do Le Bail Fit!  
            lbfitok = True
            minimizer = 'Levenberg-MarquardtMD'
            if istep == numsteps:
                minimizer = 'Simplex'
            try: 
                LeBailFit( 
                        InputWorkspace=mtd[self.datawsname],  
                        OutputWorkspace=fitwsname,  
			InputParameterWorkspace=inputparamwsname,  
                        OutputParameterWorkspace=newparameterwsname, 
			InputHKLWorkspace=mtd[self.hklwsname], 
                        OutputPeaksWorkspace=self.peakswsname, 
			BackgroundType=self.backgroundtype, 
                        BackgroundParametersWorkspace=mtd[self.bkgdparameterwsname], 
			Function='LeBailFit',
			UseInputPeakHeights='0',
                        PeakRadius='8',
                        Minimizer=minimizer) 
            except Exception:
                lbfitok = False
		print "Catch an error from LeBailFit (fitting) and continue for the rest"
                pass

            # vii.  Pass Result
            chi2s = None
            if lbfitok is True:
                try:
                    physicalchi2, chi2s = self.parseFitResult(newparameterwsname)
                except KeyError:
                    print "[Error] LeBailFit/Fit does not succeed. "
                    physicalchi2 = False
                    lbfitok = False
                    pass

            # viii. Process result
            if lbfitok is True: 
                # Update input parameter workspace
                newfitchi2 = chi2s[1]
                lbfitresults[istep] = (newparameterwsname, newfitchi2, fitwsname) 
                inputparamwsname = newparameterwsname

                wbuf = "[FIT] Chi^2 (Fit Step %d) = %.5E\n" % (istep, newfitchi2)
                for parname in sorted(self.tempparameterdict.keys()):
                    parameter = self.tempparameterdict[parname]
                    if parameter.fit is True: 
                        wbuf += "%s = %.5E, \t Error = %.5E;\n" % (parname, parameter.value, parameter.chi2)
                    else:
                        if abs(parameter.chi2) > 0.0:
                            raise NotImplementedError("Parameter %s is not for fit.  But its Chi^2 = %.5E" % (parname, parameter.chi2))
                    # ENDIFELSE
                # ENDFOR
                print wbuf

            else:
                # Fit error.  Discard the effort! 
                chi2slist.append(chi2s)
                pass 

            time.sleep(0.1)

        # ENDFOR: isteps 
        
        print "ENDSECTION2--------------------------------------------------------------\n\n"

	# 4. Process the fitting/calculation result
        thisbestchi2 = newcalchi2
        thisbestwsname = newcalparameterswsname
        case = 2

        for istep in lbfitresults.keys():
            stepresult = lbfitresults[istep]
            wsname = stepresult[0]
            tmpchi2 = stepresult[1]
            if math.isnan(tmpchi2) is False and tmpchi2 < thisbestchi2:
                thisbestchi2 = tmpchi2
                thisbestwsname = wsname
                case = 1
            # ENDIF
        # ENDFOR

        # 5. Determine the best fit-data-workspace
        self.fitwsname = ""
        bestfitchi2 = startchi2s[1]
        if math.isnan(bestfitchi2) is True:
            bestfitchi2 = 1.0E20
        beststep = -1

        for istep in xrange(numsteps+1):
            try:
                fitresult = lbfitresults[istep]
                tmpchi2 = fitresult[1]
        
                if math.isnan(tmpchi2) is False:
                    fitwsname = fitresult[2]

                    if self.fitwsname == "" or tmpchi2 < bestfitchi2:
                        # Init case of better case
                        self.fitwsname = fitwsname
                        bestfitchi2 = tmpchi2
                        beststep = istep
                    # ENDIF
                # ENDIF
            except KeyError:
                pass
        # ENDFOR is True
        bestfitchi2

        # If all NAN, then use the first one
        if self.fitwsname == "":
            self.fitwsname = lbfitresults[0][2]

        print "Best Fit Chi2 Data Workspace = %s   Local Chi^2 = %.5E  @ Step %d " % (self.fitwsname, bestfitchi2, beststep)

        if thisbestchi2 < self.curchi2:
            better = True
        else:
            better = False

        return (better, case, thisbestwsname, newcalchi2, bestfitchi2)

    #**********************  Parse Fitting Result  ****************************
    #
    def parseFitResult(self, parameterwsname):
        """ Parse fitting result to a data stucture for storing. 
        The results are stored in instance variables including
        - newchi2
        - calparameterdict/fitparameterdict

        Arguments:
         - fittype: calculation, fit

        Return 
         - meaningful chi2
        """
        parameters = self.tempparameterdict
        
        chi2s = self.parseParameterWorkspace(parameterwsname, parameters)
        physicalchi2 = True
        gslchi2 = chi2s[0]
        localchi2 = chi2s[1]
        fitchi2 = chi2s[2]

        if localchi2 < 0.0: 
            physicalchi2 = False
            print "TableWorkspace %s does not have row for Chi2 or Chi2 is meaningless." % \
                    (parameterwsname)
        else:
            pass
            # self.newchi2 = localchi2
            # if fittype == "calculation": 
            #     self.newcalchi2 = localchi2
            # elif fittype == "fit":
            #     self.newfitchi2 = localchi2
            # else:
            #     errmsg = "Fit type = %s is not supported. " % (fittype)
            #     raise NotImplementedError(errmsg)

        # print "[Sub] Parse Fit Result Type = '%s'   Chi2 = %f" % (fittype, localchi2)

        return (physicalchi2, chi2s)
    
    
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
        haschi2 = colnamedict.has_key("chi^2")

        # 3. Parse
        gslchi2 = -1.0
        fitchi2 = -1.0
        localchi2 = -1.0

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

                # Chi2
                if haschi2 is True:
                    temp = tablews.cell(irow, colnamedict["chi^2"])
                    parameterdict[parname].chi2 = temp

            elif parname == "Chi2":
                gslchi2 = tablews.cell(irow, colnamedict["Value"])

            elif parname == "FitChi2":
                fitchi2 = tablews.cell(irow, colnamedict["Value"])
                    
            elif parname == "LocalChi2":
                localchi2 = tablews.cell(irow, colnamedict["Value"])
                    
        # ENDFOR

        print "[FYI] Chi2:  GSL = %.5E,  Local = %.5E, GSL Fit = %.5E" % (gslchi2, localchi2, fitchi2)

        return (gslchi2, localchi2, fitchi2)

    #------------------------- END Parse Fitting Result -----------------------

    #********************** MC Codes  ****************************************
    #
    def randomWalk(self, seed, maxIteration=100, startiteration=1):
        """ Change parameters randomly on one random seed
        Inputs
         - seed:  random seed

        Return:
         - iteration:  number of interation of this random walk
                       -1  starting values are wrong! 
        """
        # 1. Work on the input configuration
        print "\n==========  MC Iteration %d of %d ===============\n" % (0, maxIteration)

        # a) Do a simple LeBail Fit 
        result = self.doLeBailFit(True)

        # b) Get output
        better = result[0]
        case = result[1]
        bestparwsname = result[2]
        newcalchi2 = result[3]
        newfitchi2 = result[4]

        if better is False:
            print "[Error] Starting proposed value is not physical.  Try a better a better starting value!"
            raise NotImplementedError("Starting values are rejected!")

        self.acceptNewConfiguration(bestparwsname)
        try:
            CloneWorkspace(InputWorkspace=self.calwsname, OutputWorkspace="StartCalculatedData")
        except Exception:
            pass
        try:
            CloneWorkspace(InputWorkspace=self.fitwsname, OutputWorkspace="StartFittedData")
        except Exception:
            pass
        if case == 1: 
            self.storeResult(self.parameterdict, 0, self.fitwsname)
        else:
            self.storeResult(self.parameterdict, 0, self.calwsname)
        # ENDIF

        startiteration += 1

        # 2. Start random number
        random.seed(seed) 
        self.lastupdate = 0
    
        # 3. Start MC random walk
	iteration = startiteration
        while iteration <= maxIteration:
            print "\n\n==========  MC Iteration %d of %d ===============\n" % (iteration, maxIteration)
            result = self.moveOneMonteCarloStep()
	    better = result[0]
	    case = result[1]
            parwsname = result[2]
            newcalchi2 = result[3]
            newfitchi2 = result[4]

            if better:
                # Better result.  Accept the change
                self.acceptNewConfiguration(parwsname)
		acceptnewvalues = True 
                if case == 1:
                    self.storeResult(self.tempparameterdict, iteration, self.fitwsname)
                else: 
                    self.storeResult(self.tempparameterdict, iteration, self.calwsname) 
                
		print "+++  Better Chi2 = %f of Case %d.  Accept the MC proposed step" % (self.curchi2, case)
    
            else:
                # Worse result.   Decide randomly 
                if case < 0:
                    accept = False
                else: 
                    mcresult = self.makeMCChoice(case, newfitchi2=newfitchi2, newcalchi2=newcalchi2)
                    accept = mcresult[0]

                if accept is True:
		    # Accept the proposed new value
                    self.acceptNewConfiguration(parwsname)
		    acceptnewvalues = True
		    print "+++  Worse Chi2 comparing to current Chi2 = %.4E of Case %d.  Accept the MC proposed step" % (self.curchi2, case)
		else:
		    # Reject the proposed new value
		    acceptnewvalues = False
		    print "+++  Worse Chi2 comparing to current Chi2 = %.4E of Case %d.  Reject the MC proposed step" % (self.curchi2, case)
                # ENDIF
            # ENDIFELSE

	    self.recordChi2.append(self.curchi2)
	    self.recordCase.append(case)
	    self.recordTake.append(acceptnewvalues)

            if acceptnewvalues is True:
                self.lastupdate = iteration

            # Y. Write temp file to update situation
            wbuf = "Current Iteration : %d    Last Update:  %d   Random Seed = %d\n" % (iteration, self.lastupdate, seed)
            maxi = min(11, len(self.recordChi2)+1)
            for i in xrange(1, maxi): 
                wbuf += "Iteration %d,  Chi2 = %.5E,  Case = %d, Taken = %s\n" % (iteration+1-i, self.recordChi2[-i], 
                        self.recordCase[-i], str(self.recordTake[-i]))
            ofile = open("realtimerecord.txt", "w")
            ofile.write(wbuf)
            ofile.close()
            
            # X. If no update for a certain amount of cycle, i.e., stuck in local min, 
            if iteration-self.lastupdate > _MAXNONUPDATE:
                # Increase step size by 2.0 and recount
                self.StepSizeFactor *= 2.0
                self.lastupdate = iteration
            elif iteration == self.lastupdate and self.StepSizeFactor > 1.0:
                # If there is change, then move step size factor back to 1.0 
                self.StepSizeFactor = 1.0
            
            iteration += 1
        # ENDWHILE

        # Process output
        print "\n\n==========  MC Result Until Iteration %d ===============\n" % (iteration)
	print "Walk %d Steps. Best Chi2 = %f @ Step %d" % (maxIteration, self.bestchi2, self.bestchi2step)
        for i in xrange(len(self.bestchi2s)):
            chi2 = self.bestchi2s[i]
            print "%d:  chi2 = %f    @ Step %d" % (i, chi2, self.bestresults[chi2][1])
        print "Visit storeResults() %d times." % (self.numvisitstore)

        # MC history file
        wbuf = ""
        for i in xrange(len(self.recordChi2)):
            wbuf += "%i\t\t%f\t\t%s\t%d\n" % (i, self.recordChi2[i], self.recordCase[i], self.recordTake[i])
        rfilename = "montecarlo_record_seed%s_iter%d.dat" % (seed, iteration)
        rfile = open(rfilename, "w")
        rfile.write(wbuf)
        rfile.close()

        # Parameter file
        bestchi2 = self.bestchi2s[0]
        wbuf = ""
        # Chi^2
        wbuf += "Chi2\t\t"
        for chi2 in sorted(self.bestresults.keys()):
            wbuf += "%f\t\t" % (chi2)
        wbuf += "\n"
        # Step
        wbuf += "Step\t\t"
        for chi2 in sorted(self.bestresults.keys()): 
            newtuple = self.bestresults[chi2]
            wbuf += "%d\t\t" % (newtuple[1])
        wbuf += "\n"
        # Parameters
        parnames = sorted(self.bestresults[bestchi2][0].keys())
        for parname in parnames:
            tempbuf = "%s\t\t" % (parname)
            for chi2 in sorted(self.bestresults.keys()):
		newtuple = self.bestresults[chi2]
		parvalue = newtuple[0][parname]
                # tempbuf += "%f\t\t" % (parvalue.value)
                tempbuf += "%f\t\t" % (parvalue)
            tempbuf += "\n"
            wbuf += tempbuf
        #ENDFOR
        pfilename = "bestparameters_seed%d.dat" % (seed)
        pfile = open(pfilename, "w")
        pfile.write(wbuf)
        pfile.close()
    
        return iteration


    
    def moveOneMonteCarloStep(self):
        """ Move one step in MC

        Return
         - better? :  boolean to show new chi2 is better
         - case    :  integer to show (1) fit is better or (2) cal is better
        """
        # 1. Determine how  many parameters to walk
        numparms = random.randint(1, len(self.parametersToFit))
    
        # 2. Determine parameter index to change value
        # print "Number of parmeters to have valued changed = %d out of %d" % (numparms, len(self.parametersToFit))
        paramindexes = []
        for ir in xrange(numparms):
            fitparindex = random.randint(0, len(self.parametersToFit)-1)
            fitparname = self.parametersToFit[fitparindex]
            parindex = self.paramNameIndexDict[fitparname]
            # print "  Parameter %s is selected for random walk.  Index = %d  (check name = %s)" % (fitparname, parindex, self.parameternames[parindex])
            if paramindexes.count(parindex) == 0:
                paramindexes.append(parindex)
            # ENDIF
        # ENDFOR
    
        # 3. Determine the new values
        for paramindex in xrange(len(self.parameterdict)):
            if paramindex in paramindexes:
                # choosed parameter
                newvalue = self.proposeNewValue(paramindex)
		parname = self.parameternames[paramindex]
                self.parameterdict[parname].valuechange = True
                self.parameterdict[parname].newvalue = newvalue
            else:
                # parameter ignored
		parname = self.parameternames[paramindex]
                self.parameterdict[parname].valuechange = False
            # ENDIF
        # ENDFOR

        # 4. Determine the parameters to fit
        self.fitCandidates = []

        fitsamesetpar = random.randint(0, 1)
        if fitsamesetpar == 0:
            # fit with same set of parmeters with value changes
            for paramindex in paramindexes:
                parname = self.parameternames[paramindex]
                self.fitCandidates.append(parname)
            # ENDFOR
        else:
            # re-select a set of parameters to fit
            numparms = random.randint(1, len(self.parameternames))
            for ir in xrange(numparms):
                parindex = random.randint(0, len(self.parameternames)-1)
                parname = self.parameternames[parindex]
                if self.fitCandidates.count(parname) == 0:
                    self.fitCandidates.append(parname)
                # ENDIF
            # ENDFOR
        # ENDIF
    
        # 5. Do LeBailFit
        result = self.doLeBailFit(True)
    
        return result
    
    
    def makeMCChoice(self, case, newfitchi2, newcalchi2):
        """ Make a MC choise to decide whether to take the change
        of parameters

        Return: (boolean, value of new chi2)
        """
	# 1. Survey
	if case == 1: 
            # fit is better
	    if newfitchi2 > newcalchi2:
	        raise NotImplementedError("Impossible to have fit-chi2 = %.5E > cal-chi2 = %.5E as case = 1" % (
		    newfitchi2, self.newcalchi2))
	    else:
		newchi2 = newfitchi2
	
	elif case == 2: 
            # cal is better 
            if newfitchi2 < newcalchi2: 
                raise NotImplementedError("Impossible to have fit-chi2 = %.5E < cal-chi2 = %.5E as case = 1" % 
                        (newfitchi2, newcalchi2))
	    else:
		newchi2 = newcalchi2

	else: 
            # undefined case
	    raise NotImplementedError("Case %d is not supported here.  Invalid argument!" % (case))

        # ENDIFELSE

	# 2. MC selection
        deltachi2 = newchi2 - self.curchi2
        ratio = math.exp(-1.*deltachi2/self.curchi2)
        if ratio > 1.0: 
            errmsg = "Ratio is larger than 1.  Chi2(old) = %f, Chi2(new) = %f!!!, Weird not ruled out in previous filters.  But won't affect result!" % \
                    (self.curchi2, self.newfitchi2) 
            raise NotImplementedError(errmsg)

        # 3. Generate a random number to make the decision
        choice = random.random()

        print "[MC] New Chi2 = %.5E.   Current Chi2 = %.5E.  Ratio = %.5E.  Choice = %.3E" % (newchi2, self.curchi2, ratio, choice)

        if choice < ratio:
            return (True, newchi2)
        else: 
            return (False, newchi2)
    
    
    def acceptNewConfiguration(self, paramwsname):
        """ Accept new configuration by
        (1) Apply the change of the parameters to self.parameterdict

        Setup 
        1. self.curchi2

	Arguments
	 - paramwsname:  the workspace to update
        """
        # 1. Check input workspace
        try:
            parws = mtd[paramwsname]
        except KeyError:
            errmsg = "Parameter (table) workspace %s does not exist. " % (paramwsname)
            raise NotImplementedError(errmsg)

        # 2. Parse from the workspace 
        newchi2s = self.parseParameterWorkspace(paramwsname, self.tempparameterdict)
        gslchi2 = newchi2s[0]
        localchi2 = newchi2s[1]

        if math.isnan(localchi2) is True:
            errmsg = "Client is not allowed to input chi2 equal to NAN (%.5E)" % (localchi2)
            raise NotImplementedError(errmsg)

        # 3. Update value
        msg = ""
        newline = False
        for parindex in xrange(len(self.parameternames)):
	    parname = self.parameternames[parindex]

            origvalue = self.parameterdict[parname].value
            newvalue = self.tempparameterdict[parname].value

            if abs(origvalue-newvalue) > 1.0E-5:
                # An updated value!
                self.parameterdict[parname].value = newvalue
                self.parameterdict[parname].newvalue = newvalue
                self.parameterdict[parname].valuechange = True
		
		msg += "%s = %.5E  (From %.5E, Delta = %.5E); " % (parname, newvalue, origvalue, newvalue-origvalue)
                if newline is True:
                    msg += "\n"
                newline = not newline
            # ENDIF 
        # ENDFOR 
        if newline is True:
            msg += "\n"
        print "SECTION 3 ******************************************"
        print "Fit/Calculate Result:\n%s" % (msg)
        print "****************************************************"

        # 4. Chi2
        self.curchi2 = localchi2
    
        return

    #************************  Book Keep Methods  **************************
    # 
    def storeResult(self, parametersdict, iteration, outdatawsname):
        """ Store result...
        """
        self.numvisitstore += 1

        # 1. The BEST chi2 
        if self.curchi2 < self.bestchi2: 
            self.bestchi2 = self.curchi2 
            self.bestchi2step = iteration
            CloneWorkspace(InputWorkspace=self.fitwsname, OutputWorkspace="BestLeBailData")

        # 2. Store first
        simpardict = {}
        for parname in parametersdict:
            parvalue = parametersdict[parname].value
            simpardict[parname] = parvalue
        
        self.bestresults[self.curchi2] = (simpardict, iteration)

        # 3. Remove the last if too many 
        self.bestchi2s = sorted(self.bestresults.keys()) 
        if len(self.bestchi2s) > self.numresults: 
            self.bestresults.pop(self.bestchi2s[-1])
            self.bestchi2s.pop()

        # ENDIF


        return


    def proposeNewValue(self, paramindex):
        """ use MC to generate a new value

        newvalue = currentvalue + r * stepsize.  r is random number (-0.5, 0.5)
	"""
	parname = self.parameternames[paramindex]
        parameter = self.parameterdict[parname]
        
        stepsize = parameter.stepsize*self.StepSizeFactor

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
        print "Parse Parameter Monte-Carlo File %s" % (parameterMCfilename)

        # 1. Open file
        ifile = open(parameterMCfilename, "r")
        lines = ifile.readlines()
        ifile.close()

        # 2. Parse
        wbuf = "Set Instrument Profile Parameters Monte Carlo Parameters\n"
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
            # UpdatePeakParameterTableValue(InputWorkspace=self.peakparamwsname,
            #         Column='Min', ParameterNames=[parname], NewFloatValue=minvalue)
            # UpdatePeakParameterTableValue(InputWorkspace=self.peakparamwsname,
            #         Column='Max', ParameterNames=[parname], NewFloatValue=maxvalue)
            # UpdatePeakParameterTableValue(InputWorkspace=self.peakparamwsname,
            #         Column='StepSize', ParameterNames=[parname], NewFloatValue=stepsize)

            wbuf += "%s:\t Range = %.5E\t,  %.5E\t\tStep = %f\n" % (parname, minvalue, maxvalue, stepsize)
        # ENDFOR
        print wbuf

        # print "-------------   Debug x338  ---------------"
        # tablews = mtd[self.peakparamwsname]
        # print "Number of rows = %d" % (tablews.rowCount())
        # for i in xrange(tablews.rowCount()):
        #     wbuf = ""
        #     for j in xrange(5):
        #         wbuf += "%s,  " % (str(tablews.cell(i, j)))
        #     print wbuf
        # print "-------------   Debug x338  ---------------"
        # raise NotImplementedError("Debug Stop x338")

        return
    #ENDDEF setupMCParameters(self, parameterMCfilename)

    def updatePeakParameterTableValue(self, InputWorkspace, Column, ParameterNames=None,
            NewFloatValue=0.0, NewStringValue="", OutputWorkspace=""):
        """ Update peak parameter table values... this will be a simplied verison of 
        the corresponding Pythong algorithm in order to share the 
        """

        # 0. Convert names
        try:
            tablews = mtd[str(InputWorkspace)]
        except KeyError:
            print "[Error] TableWorkspace %s does not exist.  " % (str(InputWorkspace))
            return False

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

        return True

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
  

#--------------------------------------------------------------------------
def main(argv):
    """ Main
    """
    """ ****************   Global Setup  ******************** """
    if len(argv) < 3:
        print "%s [Bank ID] [Random Seed] [Max Iteration]" % (argv[0])          
        return

    else:
        bankid = int(argv[1])
        randomseed = int(argv[2])
        if len(argv) >= 4:
            maxiteration = int(argv[3])
        else: 
            maxiteration = 20

    datadict = { 
            1: {"run": 10808},
            2: {"run": 10809},
            3: {"run": 10810}
            }

    runnumber = datadict[bankid]["run"]
    
    # bankid = 2
    # runnumber = 10809
    
    parameterstofit = ["alph*", "beta*", "sig*"]
    curiteration = 0
    lastiteration = 0
    
    """ ****************   Execution  ******************** """
    globalrecord = []
    
    # 1. Init LeBail
    mc = LeBailMonteCarlo(bankid, runnumber)
    mc.initLeBailFit() 

    # 2. Set up parameters to fit
    mc.setParametersToFit(parameterstofit)
    wbuf = "Parameters To Fit: "
    for parname in mc.parametersToFit:
        wbuf += " %s, " % (parname)
    print wbuf
   
    # 3. Random walk 
    if curiteration == 0:
        calinitconfig = True
    else:
        calinitconfig = False
    
    lastieration = curiteration
    curiteration = mc.randomWalk(randomseed, maxiteration, curiteration)
    
    if curiteration < 0:
        raise NotImplementedError("No random walk. Starting parameters are too bad. ")
    
    if lastieration == curiteration:
        raise NotImplementedError("Coding error!")
    
    bestchi2 = mc.bestchi2
    bestchi2step = mc.bestchi2step
    
    globalrecord.append( (bestchi2, bestchi2step, lastiteration, curiteration) )
    
    
    # print "---------------------------  Global Result  ------------------------"
    # index = 0
    # for tup in globalrecord:
    #     bc = tup[0]
    #     bs = tup[1]
    #     it0 = tup[2]
    #     itf = tup[3]
    #     print "Round %d:  Best Chi2 = %.5E   @ Step = %d   Iteration: from %d to %d" % (index, bc, bs, it0, itf)
    #     index += 1

    return

#--------------------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    argv = sys.argv
    
    """ defaults """
    bankid = 1
    randomseed = 0
    maxiteration = 2 
    """ end def defaults """
    
    if len(argv) == 1:
	    argv.append(bankid)
	    argv.append(randomseed)
	    argv.append(maxiteration)
    main(sys.argv)
