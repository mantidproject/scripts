#!/usr/bin/env python

import sys
import os
import imp
import stat
import shutil as shutil

from XMLparser import XMLparser
from mantid import *
from mantid.simpleapi import *
from MaskBTP import *
from MaskAngle import *
from numpy import *



def GetPathFromRunNumber(instrument,run):
    if instrument in ['ARCS','CNCS','SEQUOIA','HYSPEC']:
        #First try find_data
        finder=FileFinder.Instance()
            
        try:
            #path = finder.findRuns(str(run))[0] - ORIGINAL
            path = finder.findRuns(instrument+str(run))[0]
            return path
        except:
            pass

        instrumentsn = config.getFacility().instrument(instrument).shortName()
        f = os.popen('findnexus --event -i %s %s' % (instrumentsn,run), 'r')
        path = f.readline().strip()

        #check if there is an "ERROR" in the output string
        if path.find('ERROR') == -1 :
            return path
        else:
            raise ValueError("Event Nexus file not found")
    else:
        raise ValueError("Instrument not yet implemented")



class dgsreduction(object):
    def __init__(self, XMLfile=None):
        self.instrument=None
        self.filterbadpulses=None
        self.vanruns=None
        self.units=None
        self.vanmin=None
        self.vanmax=None
        self.processedfilename=None
        self.maskfilename=None
        self.mask=None
        self.normalizedcalibration=None
        self.ipts=None
        self.runs=None
        self.efixed=None
        self.t0=None
        self.calce=None
        self.eimon1=None
        self.eimon2=None
        self.emin=None
        self.emax=None
        self.loadmon=True
        self.ebin=None
        self.qstep=None
        self.kiokf=None
        self.tibg=None
        self.tibgstart=None
        self.tibstop=None
        self.grouping=None
        self.powderanglestep=None        
        self.goniometermotor=None
        self.goniometermotoroffset=None
        self.goniometermotoraxis=None
        self.goniometermotordirection=None
        self.save=None
        self.friendlyname=None
        self.calibrationtext=''
        self.datatext=''
        self.scantype = 'single'
        self.logvalue=''
        self.logvaluestep = None
        self.logvaluemin = None
        self.logvaluemax = None
        self.logconststep =True
        self.friendlynamelogs = None
        self.vanpath = ''
        self.datapath = ''


        if XMLfile!=None:
            if not(os.path.isfile(XMLfile)):
                raise IOError ("data text file "+ XMLfile+ " not found")
            self.RunFromXML(XMLfile)    


    def RunFromXML(self,filename):      
        #TODO check the xml against a schema
        parsed=XMLparser(filename)
        if parsed.root.tag !='dgsreduction':
            raise RuntimeError("This is not a dgsreduction xml file")
 
        if parsed.calibdict.has_key('instrument'):
            self.LoadInstrumentSettings(parsed.calibdict['instrument'])
        else:
            raise RuntimeError("No instrument defined in the XML file")
        self.LoadParameters(parsed.calibdict)

        if len(self.vanpath) > 1:
            config.appendDataSearchDir(self.vanpath)


        self.PerformCalibration()
        if parsed.datadicts!=[]:    
            for d in range(len(parsed.datadicts)-1):
                dgsi=dgsreduction()
                dgsi.LoadInstrumentSettings(parsed.calibdict['instrument'])
                dgsi.LoadParameters(parsed.calibdict)
                dgsi.PerformCalibration() 
                dgsi.LoadParameters(parsed.datadicts[d])
                if len(self.datapath) > 1:
                    config.appendDataSearchDir(self.datapath)
                dgsi.Execute()

            self.LoadParameters(parsed.datadicts[-1])
            self.Execute()


    #Where the data are actually reduced.
    def Execute(self):
        self.loadmon = self.calce
        self.resetEnergyToNone=False            
        if self.efixed==None:
            self.resetEnergyToNone=True
        if self.scantype == 'single':
            #load and filter bad pulses
            #load and add the runs together.
            #get the path for the first file
            path = GetPathFromRunNumber(self.instrument,self.runs[0])
            #load the file.
            Load(Filename=path,OutputWorkspace = 'data')
            print "Datafile "+path+" loaded."
            if self.loadmon:            
                LoadNexusMonitors(Filename=path,OutputWorkspace = 'monitorws')

            if self.filterbadpulses:
                FilterBadPulses(InputWorkspace = 'data', OutputWorkspace = 'data')

            self.datatext += "Loaded data run from "+path +"\n"

            if len(self.runs) > 1:
                for i in range(1,len(self.runs)):
                    path = GetPathFromRunNumber(self.instrument,self.runs[i])
                    Load(Filename=path,OutputWorkspace = 'datatemp')
                    print "Datafile "+path+" loaded."
                    if self.loadmon:                    
                        LoadNexusMonitors(Filename=path,OutputWorkspace = 'monitortemp')
                    if self.filterbadpulses:
                        FilterBadPulses(InputWorkspace = 'datatemp', OutputWorkspace = 'datatemp')
                    Plus(LHSWorkspace='data', RHSWorkspace = 'datatemp', OutputWorkspace='data')
                    if self.loadmon:                    
                        Plus(LHSWorkspace='monitorws', RHSWorkspace = 'monitortemp', OutputWorkspace='monitorws')
                    self.datatext += "Added data run from "+path +"\n"

            if self.filterbadpulses:
                self.datatext += "Bad pulses have been filterd from the data file(s).\n"

            #This is where the reduction is done.
            self.ProcessWorkspace('data')

        if self.scantype == 'step':
            #load each file and process individually, ONE summary file.
            for run in self.runs:
                #get the path for the first file
                path = GetPathFromRunNumber(self.instrument,run)
                #load the file.
                Load(Filename=path,OutputWorkspace = 'data')
                print "Datafile "+path+" loaded."
                if self.loadmon:
                    LoadNexusMonitors(Filename=path,OutputWorkspace = 'monitorws')
                if self.filterbadpulses:
                    FilterBadPulses(InputWorkspace = 'data', OutputWorkspace = 'data')

                self.datatext += "Loaded data run from "+path +"\n"
                if self.filterbadpulses:
                    self.datatext += "Bad pulses have been filterd from the data file(s).\n"

                #This is where the reduction is done.
                self.ProcessWorkspace('data')
                if self.resetEnergyToNone:
                    self.efixed=None

        if self.scantype == 'sweep':
            #load and filter bad pulses
            #load and add the runs together.
            #get the path for the first file
            path = GetPathFromRunNumber(self.instrument,self.runs[0])
            #load the file.
            Load(Filename=path,OutputWorkspace = 'data')
            print "Datafile "+path+" loaded."
            if self.loadmon:
                LoadNexusMonitors(Filename=path,OutputWorkspace = 'monitorws')
            if self.filterbadpulses:
                FilterBadPulses(InputWorkspace = 'data', OutputWorkspace = 'data')

            self.datatext += "Loaded data run from "+path +"\n"

            if len(self.runs) > 1:
                for i in range(1,len(self.runs)):
                    path = GetPathFromRunNumber(self.instrument,self.runs[i])
                    Load(Filename=path,OutputWorkspace = 'datatemp')
                    print "Datafile "+path+" loaded."
                    if self.loadmon:                    
                        LoadNexusMonitors(Filename=path,OutputWorkspace = 'monitortemp')
                    if self.filterbadpulses:
                        FilterBadPulses(InputWorkspace = 'datatemp', OutputWorkspace = 'datatemp')
                    Plus(LHSWorkspace='data', RHSWorkspace = 'datatemp', OutputWorkspace='data')
                    if self.loadmon:                    
                        Plus(LHSWorkspace='monitorws', RHSWorkspace = 'monitortemp', OutputWorkspace='monitorws')
                    self.datatext += "Added data run from "+path +"\n"

            if self.filterbadpulses:
                self.datatext += "Bad pulses have been filterd from the data file(s).\n"
            

            wsrun = mtd['data'].run()
            #Check if the workspace has the variable of interest
            
            if (self.logvalue == None or wsrun.hasProperty(self.logvalue)==False):
                raise ValueError("No log value given OR the given log value was not found in the file.")

            #need to split the data by an independt variable , some log value.
            #Create the array of logvalue BOUNDARIES
            if self.logvaluemin == None:
                self.logvaluemin= array(wsrun.getProperty(self.logvalue).value).min()
            if self.logvaluemax == None:
                self.logvaluemax= array(wsrun.getProperty(self.logvalue).value).max()
            if self.logvaluestep == None:
                self.logvaluestep = self.logvaluemax - self.logvaluemin

            bounds = arange(self.logvaluemin, self.logvaluemax+self.logvaluestep, self.logvaluestep)

            #Get the time correlation correct if you set the time correlation keyword.
            #To first approximation, set the time to zero for the first.

            for i in range(len(bounds)-1):
                FilterByLogValue(InputWorkspace="data",OutputWorkspace = 'dataslice', LogName= self.logvalue,MinimumValue=float(bounds[i]) ,MaximumValue = float(bounds[i+1]))
                dataslice=mtd['dataslice']
                if dataslice.getNumberEvents()>0:
                    values=array(dataslice.run().getProperty(self.logvalue).value)
                    self.datatext+= "Processing data for "+self.logvalue+" between "+str(bounds[i])+" and "+str(bounds[i+1])+", mean="+str(values.mean())+" std="+str(values.std())+"\n"
                    self.ProcessWorkspace('dataslice')                
                    if self.resetEnergyToNone:
                        self.efixed=None


    def ProcessWorkspace(self,datawsname):
        if self.efixed == None:
            if mtd[datawsname].run().hasProperty('EnergyRequest'):
                self.efixed =  float(mean(mtd[datawsname].run().getProperty('EnergyRequest').value))
            else:
                raise ValueError("no Efixed has been set, and not found in file.")

	    #get Ei, or use Ei
        if self.calce == False:
            efixed = self.efixed
            if self.t0==None:
                if self.instrument in ['HYSPEC','CNCS']:
                    t0 = self.t0fromei(efixed,self.instrument)
                else:
                    t0 = 0.0
            else:
                t0 = self.t0
            self.datatext += "User set value of incident energy, Ei="+str(efixed)+" meV, and t0="+str(t0)+" micro-seconds.\n"
       
	    #now deal with calculating the incident energy
        else:
            #check that the monitors are in memory
            try:
                mtd['monitorws']
            except:
                raise RuntimeError("monitor workspace not found")

            self.datatext += "Incident energy is calculated from monitor data.\n"

            mon1spec=mtd['monitorws'].getInstrument().getNumberParameter("ei-mon1-spec")[0]
            mon2spec=mtd['monitorws'].getInstrument().getNumberParameter("ei-mon2-spec")[0]

            
            alg=GetEi(InputWorkspace="monitorws",Monitor1Spec=int(mon1spec),Monitor2Spec=int(mon2spec),EnergyEstimate=self.efixed)	        

            efixed = float(alg[0])
            t0     = float(alg[3])		
	
            self.datatext += "Ei ="+str(efixed)+" meV, t0 ="+str(t0)+" microseconds\n"

        #Now adjust the data by the value of t0
        ChangeBinOffset(InputWorkspace=datawsname,OutputWorkspace=datawsname,Offset=-t0)

        #define Erange
        if self.emin == None:
            emin = -0.5*efixed
        else:
            emin = self.emin

        if self.emax == None:
            emax = efixed
        else:
            emax = self.emax

        if self.ebin == None:
            ebin = (emax-emin)/100.0
        else:
            ebin = self.ebin

        Erange = str(emin)+","+str(ebin)+","+str(emax)    


	    #Time-ind-bg subtraction.
        if (self.tibg):
            #check if tibmin and tibmax have been defined.
            tibmin = self.tibgstart
            tibmax = self.tibgstop
            if tibmin== None or tibmax == None:
                raise ValueError("Time independent background subtraction selected, but no limits set.")
            tibstep=tibmax-tibmin
            tibpar=str(tibmin)+","+str(tibstep)+","+str(tibmax)

            Rebin(InputWorkspace=datawsname,OutputWorkspace="background_origin_ws",Params=tibpar,PreserveEvents=False)
            ConvertUnits(InputWorkspace=datawsname,OutputWorkspace=datawsname,Target="DeltaE",EMode="Direct",Efixed=efixed)

	        #Do the Binning into energy bins
            Rebin(InputWorkspace=datawsname,OutputWorkspace=datawsname,Params=Erange,PreserveEvents=False)
            ConvertUnits(InputWorkspace=datawsname,OutputWorkspace=datawsname,Target="TOF",EMode="Direct",Efixed=efixed)

            ConvertToDistribution(Workspace=datawsname)
            FlatBackground(InputWorkspace="background_origin_ws",OutputWorkspace="background_ws",StartX=tibmin,EndX=tibmax,Mode="Mean",OutputMode="Return Background")
            ConvertToDistribution(Workspace="background_ws")
            Minus(LHSWorkspace=datawsname,RHSWorkspace="background_ws",OutputWorkspace=datawsname)
            ConvertFromDistribution(Workspace=datawsname)
            self.datatext  += "Time-independent background between "+str(tibmin)+" and "+str(tibmax)+" microseconds was subtracted.\n"
        else:
            self.datatext += "No time-independent background subtraction performed.\n"


        #normalize by charge
        NormaliseByCurrent(InputWorkspace=datawsname,OutputWorkspace=datawsname)
        w=mtd[datawsname]
        totalmuAhr = w.run().getProtonCharge()
        totalcoul  = totalmuAhr/1000*3.6
        self.datatext += "Data normalized by proton charge ("+str(totalmuAhr) + " micro-Ah), ("+str(totalcoul)+" C).\n"


        #detector wavelength sensitivity
        ConvertUnits(InputWorkspace=datawsname,OutputWorkspace=datawsname,Target="Wavelength",EMode="Direct",EFixed=efixed)
        He3TubeEfficiency(InputWorkspace=datawsname,OutputWorkspace=datawsname)												
        self.datatext += "Data corrected for He3 Tube Efficiency.\n"

        #Convert the data to units of energy
        ConvertUnits(InputWorkspace=datawsname,OutputWorkspace=datawsname,Target="DeltaE",EMode="Direct",EFixed=efixed)

        
	    #ki/kf
        if self.kiokf == True:
            CorrectKiKf(InputWorkspace=datawsname,OutputWorkspace=datawsname)
            self.datatext += "ki/kf factor has been applied to the data.\n" 

        #Rebinning the data if it was not already done.
        Rebin(InputWorkspace=datawsname,OutputWorkspace=datawsname,Params=Erange,PreserveEvents=False)
        self.datatext += "Data binned with emin="+str(emin)+", emax="+str(emax)+", ebin=" + str(ebin)+ " meV.\n"

        #normalize data by 1/(bin)
        #Convert to differential cross section by dividing by the energy bin width
        ConvertToDistribution(Workspace=datawsname)
        self.datatext += "Data converted to differential cross section by dividing by the energy bin width.\n"

	    #do the masking and calibration
        MaskDetectors(Workspace=datawsname,MaskedWorkspace="calibration")
        Divide(LHSWorkspace=datawsname,RHSWorkspace="calibration",OutputWorkspace=datawsname)
        self.datatext += "Data have been normalized and masked by the calibration file.\n"

	    #deal with angles
        [psiangle, angletext]= definegoniometer(self.goniometermotor, self.goniometermotoroffset, self.goniometermotordirection, self.goniometermotoraxis, datawsname)
        self.datatext += angletext

	    #grouping
        #in general grouping files will be stored in Mantid/instrument/Grouping
        #check if the requested grouping file is present
        #Two types of grouping, powder and pixel
        if ((self.grouping == 'powder') and (self.powderanglestep != None)):
            #do the powder work
            outputdir = os.path.abspath(os.curdir)
            mapping=createanglelist(datawsname,self.powderanglestep)
            GroupDetectors(InputWorkspace=datawsname,OutputWorkspace=datawsname,MapFile=outputdir+"/powdergroup.map",Behaviour="Sum")
            SolidAngle(InputWorkspace=datawsname,OutputWorkspace="sa")
            Divide(LHSWorkspace=datawsname,RHSWorkspace="sa",OutputWorkspace=datawsname)
            DeleteWorkspace(Workspace = "sa")
            self.datatext += "Detectors grouped by angle with a step of "+str(self.powderanglestep)+ ".\n"

        #case of powder grouping and NO valid angle step
        elif ((self.grouping == 'powder') and (self.powderanglestep == None)):
            raise ValueError("Powder grouping chosen, but anglestep is invalid.")

        elif (self.grouping != None and len(self.grouping)>0):
            #do the pixel grouping.
            #check if the grouping file is present

            #try looking in the current directory first
            grouppath = os.path.abspath(os.curdir)+'/'+self.instrument+"_"+self.grouping+"_grouping.xml"
            if (os.path.isfile(grouppath)):
                GroupDetectors(InputWorkspace=datawsname,OutputWorkspace=datawsname,MapFile=grouppath,Behaviour="Average")
 #               GroupDetectors(InputWorkspace=datawsname,OutputWorkspace=datawsname,MapFile=grouppath,Behaviour="Sum")
 #               SolidAngle(InputWorkspace=datawsname,OutputWorkspace="sa")
 #               Divide(LHSWorkspace=datawsname,RHSWorkspace="sa",OutputWorkspace=datawsname)
 #               DeleteWorkspace(Workspace = "sa")

                self.datatext += "Detectors grouped by averaging over "+self.grouping+" pixels.\n"
            else:
                grouppath=config["instrumentDefinition.directory"]+'Grouping/'+self.instrument+"_"+self.grouping+"_grouping.xml"
                if (os.path.isfile(grouppath)):
                    GroupDetectors(InputWorkspace=datawsname,OutputWorkspace=datawsname,MapFile=grouppath,Behaviour="Average")
#                    GroupDetectors(InputWorkspace=datawsname,OutputWorkspace=datawsname,MapFile=grouppath,Behaviour="Sum")
#                    SolidAngle(InputWorkspace=datawsname,OutputWorkspace="sa")
#                    Divide(LHSWorkspace=datawsname,RHSWorkspace="sa",OutputWorkspace=datawsname)
#                    DeleteWorkspace(Workspace = "sa")
                    self.datatext += "Detectors grouped by averaging over "+self.grouping+" pixels.\n"            
                else:
                    raise ValueError("Grouping file "+grouppath+" NOT FOUND.")


        #Now deal with saving files.
        if self.save != None:
            #create a friendly name
            friendlynamebase = self.CreateFriendlyFilename(datawsname)

            #Parse the filetypes to save
            if 'nxspe' in self.save:
                #if there is a powder maping file in the current directory, then use it with the .nxspe file
                if self.grouping == 'powder':
                    SaveNXSPE(Filename=friendlynamebase+".nxspe", InputWorkspace=datawsname, Efixed=str(efixed),Psi=str(psiangle), KiOverKfScaling=self.kiokf, ParFile=os.path.abspath(os.curdir) +"/powdergroup.par")
                else:
                    SaveNXSPE(Filename=friendlynamebase+".nxspe", InputWorkspace=datawsname, Efixed=str(efixed),Psi=str(psiangle), KiOverKfScaling=self.kiokf)
                self.datatext += "Data have been saved as a .nxspe file, FILENAME="+friendlynamebase+".nxspe.\n"
                #change permissions of the directory and file
                changepermissions(friendlynamebase+".nxspe")
 
            if 'nxs' in self.save:
                #save the nxs
                SaveNexus(Filename=friendlynamebase+".nxs", InputWorkspace=datawsname)
                self.datatext += "Data have been saved as a .nxs file, FILENAME="+friendlynamebase+".nxs.\n"
                #change permissions of the directory and file
                changepermissions(friendlynamebase+".nxs")

            #Save the diffraction pattern as I(Q) three column text format
            if 'iofq' in self.save:
                #check that qbining has been set
                if self.qstep == None:
                    #if not set then set it to 0.1 inv angstroms
                    qstep = 0.1
                else:
                    qstep = self.qstep

                #get the qmin and q max
                [qmin, qmax] = calqrangefromworkspace(datawsname)
                qbinparams = str(qmin)+','+str(qstep)+','+str(qmax)
                ws = mtd[datawsname]
                #get the energy binning.
                eminloc = ws.readX(0)[0]
                emaxloc = ws.readX(0)[-1]
                efullbin = (emaxloc - eminloc)*1.50
                ebinparams = str(eminloc)+","+str(efullbin)+","+str(emaxloc)
                SofQW(InputWorkspace=datawsname,OutputWorkspace='SofQWdata',QAxisBinning=qbinparams,Emode="Direct",Efixed=efixed)
                Transpose(InputWorkspace='SofQWdata',OutputWorkspace='SofQWdata')
                wsofq = Rebin2D(InputWorkspace='SofQWdata',Axis1Binning=qbinparams,Axis2Binning=ebinparams)
                SaveAscii(Filename=friendlynamebase+"_iofq.dat",InputWorkspace='wsofq')
                self.datatext += "Data have been saved as a iofq.dat file, FILENAME="+friendlynamebase+"_iofq.dat.\n" 
                #change permissions of the directory and file
                changepermissions(friendlynamebase+"_iofq.dat")


            if 'par' in self.save:
                if self.grouping == 'powder':
                    shutil.copy(os.path.abspath(os.curdir) +"/powdergroup.par",friendlynamebase+".par")
                else:
                    SavePAR(Filename=friendlynamebase+".par", InputWorkspace=datawsname)
                self.datatext += "PAR file has been saved as FILENAME="+friendlynamebase+".par.\n"
                #change permissions of the directory and file
                changepermissions(friendlynamebase+".par")

            if 'phx' in self.save:
                SavePHX(Filename=friendlynamebase+".phx",InputWorkspace=datawsname)
                self.datatext += "PHX file has been saved as FILENAME="+friendlynamebase+".phx.\n"

            if 'spe' in self.save:
                SaveSPE(Filename=friendlynamebase+".spe",InputWorkspace=datawsname)
                self.datatext += "SPE file has been saved as FILENAME="+friendlynamebase+".spe.\n"

            if 'iofe' in self.save:
                #get the q binning.
                #get the qmin and q max
                [qmin, qmax] = calqrangefromworkspace(datawsname)
                qfullbin = (qmax - qmin)*1.50
                qbinparams = str(qmin)+","+str(qfullbin)+","+str(qmax)
                SofQW(InputWorkspace=datawsname,OutputWorkspace='SofQWdata',QAxisBinning=qbinparams,Emode="Direct",Efixed=efixed)
                Transpose(InputWorkspace='SofQWdata',OutputWorkspace='SofQWdata')
                wsofe = Rebin2D(InputWorkspace='SofQWdata',Axis1Binning=qbinparams,Axis2Binning=Erange)
                SaveAscii(Filename=friendlynamebase+"_iofe.dat",InputWorkspace='wsofe')
                self.datatext += "Data have been saved as a iofe.dat file, FILENAME="+friendlynamebase+"_iofe.dat.\n" 
                #change permissions of the directory and file
                changepermissions(friendlynamebase+"_iofe.dat")

            #need to add i(q,e) save option
            #need to add .jpg save option.

            #writing the summary file
            if 'summary' in self.save:
                #check if the friendlyname directory exists.
                summaryfilename = friendlynamebase+"_summary.txt"
                parentdir = os.path.dirname(summaryfilename)    
                if os.path.isdir(parentdir) == False:
                    os.mkdir(parentdir,0755)
                sumfile = open(summaryfilename, 'w')
                sumfile.write("-----------VANADIUM CALIBRATION AND MASKING-----------\n")
                sumfile.write(self.calibrationtext)
                sumfile.write("\n-----------DATA REDUCTION-----------------------------\n")
                sumfile.write(self.datatext)
                sumfile.close()
                changepermissions(summaryfilename)


    def SetInstrument(self,instrument):
        if instrument not in ['ARCS','CNCS','HYSPEC','SEQUOIA']:
            raise ValueError("Instrument not defined")
        else:
            self.instrument=instrument
            self.LoadInstrumentSettings(instrument)

    def SetFilterBadPulses(self,value):
        self.fileterbadpulses=value

    def SetMask(self,algorithm=None,**kwargs):
        pass

    def LoadInstrumentSettings(self,instrumentname):
        """
        instrumentname is a string
        """
        modulename = instrumentname.lower()+'default'
        #We have the instrument name.
        #Try to import the default file for that instrument
        #First try to import from the /SNS/"instrumentname"/shared
        try:
            path = "/SNS/"+instrumentname+"/shared/"+modulename+'.py'
            if (os.path.isfile(path)):
                m = imp.load_source(modulename,path)
            else:
                m = imp.load_source(modulename,os.path.abspath(os.curdir)+"/"+modulename+'.py')
        except:
            raise RuntimeError("Could not find instrument definition module"+os.path.abspath(os.curdir)+"/"+modulename+'.py')
        params=m.instrumentparameters()
        self.LoadParameters(params)


    def t0fromei(self,ei,instrumentname):
        modulename = instrumentname.lower()+'default'
        try:
            path = "/SNS/"+instrumentname+"/shared/"+modulename+'.py'
            if (os.path.isfile(path)):
                m = imp.load_source(modulename,path)
            else:
                m = imp.load_source(modulename,os.path.abspath(os.curdir)+"/"+modulename+'.py')
        except:
            raise RuntimeError("Could not find instrument definition module"+os.path.abspath(os.curdir)+"/"+modulename+'.py')
        return m.t0fromei(ei)


    def LoadParameters(self,params):
        """
        Load parameters from a dictionary
        """
        if params.has_key('instrument'):
            if self.instrument==None:
                self.instrument=params['instrument']
            else:
                if self.instrument!=params['instrument']:
                    raise RuntimeError("instruments not compatible")
        if params.has_key('filterbadpulses'):
            self.filterbadpulses=params['filterbadpulses']
        if params.has_key('vanruns'):
            self.vanruns=params['vanruns']
        if params.has_key('units'):
            self.units=params['units']
        if params.has_key('vanmin'):
            self.vanmin=params['vanmin']
        if params.has_key('vanmax'):
            self.vanmax=params['vanmax']
        if params.has_key('processedfilename'):
            self.processedfilename=params['processedfilename']
        if params.has_key('maskfilename'):
            self.maskfilename=params['maskfilename']
        if params.has_key('mask'):
            self.mask=params['mask']
        if params.has_key('normalizedcalibration'):
            self.normalizedcalibration=params['normalizedcalibration']
        if params.has_key('ipts'):
            self.ipts=params['ipts']
        if params.has_key('runs'):
            self.runs=params['runs']
        if params.has_key('efixed'):
            self.efixed=params['efixed']
        if params.has_key('t0'):
            self.t0=params['t0']
        if params.has_key('calce'):
            self.calce=params['calce']
        if params.has_key('ei-mon1-spec'):
            self.eimon1=params['ei-mon1-spec']
        if params.has_key('ei-mon2-spec'):
            self.eimon2=params['ei-mon2-spec']
        if params.has_key('emin'):
            self.emin=params['emin']
        if params.has_key('emax'):
            self.emax=params['emax']
        if params.has_key('ebin'):
            self.ebin=params['ebin']
        if params.has_key('qstep'):
            self.qstep=params['qstep']
        if params.has_key('kiokf'):
            self.kiokf=params['kiokf']
        if params.has_key('tibg'):
            self.tibg=params['tibg']
        if params.has_key('tibgstart'):
            self.tibgstart=params['tibgstart']
        if params.has_key('tibgstop'):
            self.tibgstop=params['tibgstop']
        if params.has_key('grouping'):
            self.grouping=params['grouping']
        if params.has_key('powderanglestep'):
            self.powderanglestep=params['powderanglestep']      
        if params.has_key('goniometermotor'):
            self.goniometermotor=params['goniometermotor']
        if params.has_key('goniometermotoroffset'):
            self.goniometermotoroffset=params['goniometermotoroffset']
        if params.has_key('goniometermotoraxis'):
            self.goniometermotoraxis=params['goniometermotoraxis']
        if params.has_key('goniometermotordirection'):
            self.goniometermotordirection=params['goniometermotordirection']
        if params.has_key('save'):
            self.save=params['save']
        if params.has_key('friendlyname'):
            self.friendlyname=params['friendlyname']
        if params.has_key('scantype'):
            self.scantype=params['scantype']
        if params.has_key('logvalue'):
            self.logvalue=params['logvalue']
        if params.has_key('logvaluestep'):
            self.logvaluestep=params['logvaluestep']
        if params.has_key('logvaluemin'):
            self.logvaluemin=params['logvaluemin']
        if params.has_key('logvaluemax'):
            self.logvaluemax=params['logvaluemax']
        if params.has_key('friendlynamelogs'):
            self.friendlynamelogs=params['friendlynamelogs']
        if params.has_key('vanpath'):
            self.vanpath= params['vanpath']
        if params.has_key('datapath'):
            self.datapath = params['datapath']


    def PerformCalibration(self):
        #first check if there is currently a file ready to use
        #if the processedfilename exists in the current directory,
        #then load it and exit.
        if (self.processedfilename!=None and (os.path.isfile(self.processedfilename))):
            path=os.getcwd()+'/'+self.processedfilename
            Load(Filename = path, OutputWorkspace = 'calibration')
            self.calibrationtext += "Calibration loaded from "+path+"\n"
            print "Calibration file "+path+" loaded."
        #otherwise, this file does not exist, and it must be made.
        else:
            self.calibrationtext+="Performing calibration:\n"
            if self.vanruns == None:
                calibration=LoadEmptyInstrument(Filename=config['instrumentDefinition.directory']+self.instrument+'_Definition.xml',
                                DetectorValue='1')
                #if there is any error in loading the empty instrument, then issue an error.
                if calibration == None:
                    raise RuntimeError("Instrument definition could not be loaded")

                #remove any beam monitors from the empty instrument calibration
                i=0
                while(calibration.getDetector(i).isMonitor()):
                    i += 1
                    #i is the index of the first true detector
                #now, crop the workspace of the monitors
                calibration = CropWorkspace(calibration,StartWorkspaceIndex=i)                
    
                self.calibrationtext += "Loaded empty instrument file.   Calibration set to 1.0 \n"
                print "No calibration file loaded, calibration set to 1.0."
            else:
                path=GetPathFromRunNumber(self.instrument,self.vanruns[0])
                calibration=Load(Filename=path)
                if self.filterbadpulses:
                    FilterBadPulses(InputWorkspace = 'calibration', OutputWorkspace = 'calibration')
                self.calibrationtext += "Loaded vanadium run from "+path +"\n"
                print "Calibration file "+path+" loaded."
                if len(self.vanruns) > 1:
                    for i in range(1,len(self.vanruns)):
                        path = GetPathFromRunNumber(self.instrument,self.vanruns[i])
                        calibrationtemp=Load(Filename=path)
                        if self.filterbadpulses:
                            FilterBadPulses(InputWorkspace = 'calibrationtemp', OutputWorkspace = 'calibrationtemp')
                        calibration=calibration+calibrationtemp
                        self.calibrationtext += "Added vanadium run from "+path +"\n"
                if self.filterbadpulses:
                    self.calibrationtext+= "Bad pulses have been filtered from the vanadium calibration.\n"

                #If the file is not being normalized to one than normalise by current.    
                #Normalize the intensity by charge in milliAmphours (1mA hr = 3.6 coul)
                if self.normalizedcalibration!=True:
                    NormaliseByCurrent(InputWorkspace='calibration', OutputWorkspace='calibration')
                    self.calibrationtext += "Normalized vanadium to proton charge.\n"

                #get the total current and store it in a variable.
                totalmuAhr = mtd['calibration'].run().getProtonCharge()
                totalcoul  = totalmuAhr/1000*3.6
                self.calibrationtext += "Proton charge ("+str(totalmuAhr) + " micro-Ah), ("+str(totalcoul)+" C).\n"


                #Change the units of the calibration workspace
                if self.units != None:
                    ConvertUnits(InputWorkspace = 'calibration', OutputWorkspace = 'calibration', Target = self.units, EMode='Elastic')
                #Need to integrate the data from vanmin to vanmax
                #check the limits.
                if self.vanmin == None:
                    vanmin = mtd['calibration'].getTofMin()
                else:
                    vanmin = self.vanmin
                if self.vanmax == None:
                    vanmax = mtd['calibration'].getTofMax()
                else:
                    vanmax = self.vanmax
                #integrate the data.
                Rebin(InputWorkspace = 'calibration', OutputWorkspace = 'calibration', Params=str(vanmin)+','+str(vanmax-vanmin)+','+str(vanmax),PreserveEvents=False)
                calibration=mtd['calibration']
                self.calibrationtext += "Vanadium integrated between "+ calibration.getAxis(0).getUnit().unitID() + " " +str(vanmin)+", "+str(vanmax)+ " " +calibration.getAxis(0).getUnit().label()+"\n"
            #Do the Masking.
            #loop over the self.mask array
            for elem in self.mask:
                self.calibrationtext += "Masking vanadium with "+str(elem) +"\n"
                if elem['algorithm'].lower()=='banktubepixel':
                    myparser=XMLparser(None)
                    #get the bank tubes and pixels
                    if elem.has_key('bank'):
                        bank = myparser.parseto('bank',elem['bank'])
                    else:
                        bank = None
                    if elem.has_key('tube'):
                        tube = myparser.parseto('tube',elem['tube'])
                    else:
                        tube = None
                    if elem.has_key('pixel'):
                        pixel = myparser.parseto('pixel',elem['pixel'])
                    else:
                        pixel = None
                    #mask the btp
                    MaskBTP(Workspace='calibration',Bank=bank,Tube=tube,Pixel=pixel)
#                    btpdict = dict(elem)
#                    btpdict['Workspace'] = 'calibration'
#                    MaskBTP(**btpdict)

                if elem['algorithm'].lower()=='angle':
                    #get the ttmin and max
                    if elem.has_key('twothetamin'):
                        ttmin = float(elem['twothetamin'])
                    else:
                        ttmin = None
                    if elem.has_key('twothetamax'):
                        ttmax = float(elem['twothetamax'])
                    else:
                        ttmax = None
                    #mask the angle
                    print ttmin,ttmax
                    MaskAngle(Workspace='calibration',twothetamin=ttmin,twothetamax=ttmax)
#                    angledict = dict(elem)
#                    angledict['Workspace'] = 'calibration'
#                    MaskAngle(**angledict)

                if elem['algorithm'].lower() == 'mediandetectortest':
                    detectordiagdict = dict(elem)
                    detectordiagdict['InputWorkspace'] = 'calibration'
                    detectordiagdict['OutputWorkspace']= 'maskdetectordiag'
                    #Need to pop off the 'algorithm' keyvalue from the dict.
                    detectordiagdict.pop('algorithm')
                    MedianDetectorTest(**detectordiagdict)
                    MaskDetectors(Workspace='calibration',MaskedWorkspace='maskdetectordiag')
                    DeleteWorkspace(Workspace='maskdetectordiag')

                if elem['algorithm'].lower() == 'finddetectorsoutsidelimits':
                    detectordiagdict = dict(elem)
                    detectordiagdict['InputWorkspace'] = 'calibration'
                    detectordiagdict['OutputWorkspace']= 'maskdetectordiag'
                    #Need to pop off the 'algorithm' keyvalue from the dict.
                    detectordiagdict.pop('algorithm')
                    FindDetectorsOutsideLimits(**detectordiagdict)
                    MaskDetectors(Workspace='calibration',MaskedWorkspace='maskdetectordiag')
                    DeleteWorkspace(Workspace='maskdetectordiag')


            #Merging a prior mask with the current mask.
            if self.maskfilename != None:
                #check if the file exists
                try:
                    Load(Filename = self.maskfilename, OutputWorkspace = 'loadedmask')
                    MaskDetectors(Workspace='calibration',MaskedWorkspace='loadedmask')
                    DeleteWorkspace(Workspace='loadedmask')
                    self.calibrationtext += "Mask vanadium merged with file " + self.maskfilename+'\n'
                except:
                    raise RuntimeError("Could not load preexisting maskfile " + self.maskfilename)


            #This is the section for normalizing the vanadium calibration to fluctuate about 1.0
            if self.normalizedcalibration==True:
                #normalize the calibration file about 1.0
                #get the mean calibration intensity for non-masked points
                datay = mtd['calibration'].extractY()
                meanval = float(datay[datay>0].mean())
                CreateSingleValuedWorkspace(OutputWorkspace='meanval',DataValue=meanval)
                Divide(LHSWorkspace='calibration',RHSWorkspace='meanval',OutputWorkspace='calibration')
                DeleteWorkspace(Workspace='meanval')
                self.calibrationtext += "Calibration normalized to fluctuate about 1.0\n"      

            #save the calibration file as a loadable workspace.            
            SaveNexus(Filename = os.getcwd()+'/'+self.processedfilename, InputWorkspace = 'calibration')
            changepermissions(os.getcwd()+'/'+self.processedfilename)
            self.calibrationtext += "Saving calibration and mask to "+self.processedfilename
            
            sumfile = open(os.getcwd()+'/'+self.processedfilename+'.calibration_summary.txt', 'w')
            sumfile.write(self.calibrationtext)
            sumfile.close()
            changepermissions(os.getcwd()+'/'+self.processedfilename+'.calibration_summary.txt')


#        print self.calibrationtext


    def CreateFriendlyFilename(self,wsname):

        #access the sample log to generate the friendlyname
        friendlyfilename = os.path.abspath(os.curdir)+'/'+self.friendlyname+'/'+self.friendlyname
        if self.friendlynamelogs != None:
            #get the handle to the run
            run = mtd[wsname].run()
            for part in self.friendlynamelogs if not isinstance(self.friendlynamelogs,basestring) else [self.friendlynamelogs]:
                if run.hasProperty(part):
                    value = run.getProperty(part).value
                    try:
                        friendlyfilename += "_"+part+"_"+value
                    except:
                        #splitting up the value in case of deciaml points.  Typically, decimal points in filenames
                        #causes issues with Horace.
                        value=array(value).mean()

                        roundedvalue = "%.2f" % value
                        valuestringwithoutdot = str(roundedvalue).replace('.', 'p')
                        friendlyfilename += "_"+valuestringwithoutdot

#                        intpart = int(value)
#                        decpart = int(round(abs(value-intpart)*100))
#                        friendlyfilename += "_"+str(intpart)+"p"+str(decpart)
        return friendlyfilename



def definegoniometer(names, offsets, directions, axes, workspace):
    """
    Function for defining the goniometer
    Returns the psi value which is written to the nxspe file.
    workspace is a string
    returns [psi, text string of goniometer settings]
    """


    #get a handle to the workspace
    ws = mtd[workspace]

    outputstring = ""
    psivalue     = 0


    #Transform all inputs to LISTS.
    if str(type(offsets)) != "<type 'list'>" and str(type(offsets)) != "<type 'NoneType'>":
        offsets = [offsets]

    if str(type(names)) != "<type 'list'>" and str(type(names)) != "<type 'NoneType'>":
        names = [names]

    if str(type(directions)) != "<type 'list'>" and str(type(directions)) != "<type 'NoneType'>":
        directions = [directions]

    if str(type(axes)) != "<type 'list'>" and str(type(axes)) != "<type 'NoneType'>":
        axes = [axes]


    #first case, no motor name given (i.e. None)
    if names == None:
        if offsets == None:
            #no motor name, and no offset, just exit
            outputstring = "No goniometer set.\n"
            return [psivalue, outputstring]
        else:
           #check the number of offsets = N axes = Ndirections
            if (len(offsets)==len(directions) and len(directions)==len(axes) and len(offsets)<7):
                #because no motor names given, we need to make a fake log.
                #for loop over the angles
                #list of 6 empty strings that will be filled in
                anglelist = ["","","","","",""]
                try:
                    for i in range(len(offsets)):
                        AddSampleLog(Workspace=workspace,LogName="angle"+str(i),
                                    LogText=str(offsets[i]),LogType="Number Series")
                        anglelist[i] = "angle"+str(i)+","+axes[i]+","+str(directions[i])
                        #print anglelist[i]
                except:
                    raise RuntimeError("Could not find goniometer axis(axes)")
                SetGoniometer(Workspace=workspace,Axis0=anglelist[0],Axis1=anglelist[1],Axis2=anglelist[2],
                              Axis3=anglelist[3],Axis4=anglelist[4],Axis5=anglelist[5])
                outputstring = "The following axes have been set:\n"

                for i in range(len(offsets)):
                    tempstr = "CCW"
                    if directions[i] == -1:
                        tempstr = "CW"
                    outputstring += "   Axis"+str(i)+" along the "+ axes[i] +" direction, "+tempstr+", rotation angle="+str(offsets[i])+"\n"

                psivalue = offsets[0]
                return [psivalue, outputstring]
            else:
                raise ValueError("Number of angle offsets, directions and axes do not match.")

    #other big case, Motor name is given
    else:
        #Are there any offsets listed
        if offsets==None:
            #create offsets = 0 for ALL motor names listed
            offsets = zeros(len(names))
        #check if the noffsets = Naxes = Ndirectiosn = Nnames
        if (len(offsets)==len(directions) and len(directions)==len(axes) and len(axes)==len(names) and len(offsets)<7):
            #everything is ready
            anglelist = ["","","","","",""]
            anglevalues = []
            try:
                for i in range(len(offsets)):
                    #get the correct log from the workspace
                    angle = mean(ws.run().get(names[i]).value) + offsets[i]
                    anglevalues.append(angle)
                    AddSampleLog(Workspace=workspace,LogName="angle"+str(i),
                                LogText=str(angle),LogType="Number Series")
                    anglelist[i] = "angle"+str(i)+","+axes[i]+","+str(directions[i])
            except:
                raise RuntimeError("Could not find goniometer axis "+names[i])
            SetGoniometer(Workspace=workspace,Axis0=anglelist[0],Axis1=anglelist[1],Axis2=anglelist[2],
                              Axis3=anglelist[3],Axis4=anglelist[4],Axis5=anglelist[5])
            outputstring = "The following axes have been set:\n"

            for i in range(len(offsets)):
                tempstr = "CCW"
                if directions[i] == -1:
                    tempstr = "CW"
                outputstring += "   Axis"+str(i)+" along the "+ axes[i] +" direction, "+tempstr+", rotation angle="+names[i]+"+" +str(offsets[i])+"="+str(anglevalues[i])+"\n"

            return [anglevalues[0],outputstring]
        else:
            raise ValueError("Number of angle names, offsets, directions and axes do not match.")


def createanglelist(ws,astep):
    """
    Function to create a map of detectors corresponding to angles in a certain range.
    ws is a string refering to the workspace that contains an instrument
    astep is the step size for the angular binning (degrees).
    """

    amin = 0.0
    amax = 180.0
    bin_angles=arange(amin+astep*0.5,amax+astep*0.5,astep)
    a=[[] for i in range(len(bin_angles))] #list of list with detector IDs
    w=mtd[ws]
    origin = w.getInstrument().getSample().getPos()
    for i in range(w.getNumberHistograms()):
        ang=w.getDetector(i).getTwoTheta(origin,V3D(0,0,1))*180/pi
        index=int((ang-amin)/astep)
        if (index>=0) and (index<len(a)) and ((w.getDetector(i).getID())>0):
            a[index].append(w.getDetector(i).getID())
    #create lists with angles and detector ID only for bins where there are detectors 
    ang_list=[]
    detIDlist=[]
    for elem,ang in zip(a,bin_angles):
        if len(elem)>0:
            detIDlist.append(elem)
            ang_list.append(ang)
	# file with grouping information, saved to current directory
    outputdir = os.path.abspath(os.curdir)
    f=open(outputdir+"/powdergroup.map",'w')
    print >>f,len(ang_list)
    for i in range(len(ang_list)):
        print >>f,i
        print >>f,len(detIDlist[i])
        mystring=str(detIDlist[i]).strip(']').strip('[')
        mystring=mystring.replace(',','')
        print >>f,mystring
    f.close()
    #print ang_list
    # par file
    f=open(outputdir+"/powdergroup.par",'w')
    print >>f,len(ang_list)
    for i in range(len(ang_list)):
        print >>f,3.5,ang_list[i],0.0,1.0,1.0,1
    f.close()
    return [ang_list,detIDlist]


def calqrangefromworkspace(workspace):
    """ Direct geometry calculator for determining the max and min
    q-values that will be measured for a given workspace, ws.  The
    workspace must have a log value for incident energy called Ei.
    ws must be a 2D workspace with the x-axis energy.
    THIS IS NOT CHECKED.
    The function uses the incident energy (Ei),
    corresponding to zero energy transfer at scattering angle1
    and hwmin energy transfer at scattering angle2"""

    ws = mtd[workspace]
    
    #Get the Ei from the log
    ei = float(ws.run().get('Ei').value)

    #get angle1 and angle2 (these are the min and max scattering angles).
    angle1 = 180.0
    angle2 = 0.0

    #loop over all the detectors and find the
    for i in range(ws.getNumberHistograms()):
        #get the detectorIDs
        ids = ws.getSpectrum(i).getDetectorIDs()
        for j in ids:
            det = ws.getInstrument().getDetector(j)
            if not det.isMasked():
                angle =  degrees(det.getTwotheta(V3D(0,0,0),V3D(0,0,1)))
                if angle < angle1:
                    angle1 = angle
                if angle > angle2:
                    angle2 = angle
 
    #need the minimum energy transfer the data have been binned to.
    hwmin = ws.readX(0)[0]

    hbarsqrd_over_2m = 2.072

    Qmin = sqrt((1.0/hbarsqrd_over_2m)*(2*ei-0.0-2.0*sqrt(ei*(ei-0.0))*cos(angle1*pi/180.0)))
    Qmax = sqrt((1.0/hbarsqrd_over_2m)*(2*ei-hwmin-2.0*sqrt(ei*(ei-hwmin))*cos(angle2*pi/180.0)))
        
    return array([Qmin, Qmax])



def changepermissions(filename):
    """
    change permissions of the directory and file of filename
    to read, write for everyone
    directory also allows execute.
    """
    parentdir = os.path.dirname(filename)
    try:
        os.chmod(parentdir, stat.S_IRWXG | stat.S_IRWXO | stat.S_IRWXU)
    except:
        print "Not able to change permissions of " + parentdir

    try:
        os.chmod(filename, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH |stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH )
    except:
       print "Not able to change permissions of " + filename


def parsedbl(val):
	"""
	####PARSEDBL##########################################
	#This function parses the string and converts to a double floating
	#point value
	#
	# example:
	#  print parsedbl("123.7")
	#	  123.7
	#
	######################################################
	"""
	return float(val)

def parsebool(val):
	"""
	####PARSEBOOL##########################################
	#This function parses text to see if it is a true statement.
	#It will return true for "yes", "true", "t", "1", "tru", "tr", "y".
	#can be single or double quotes.
	#
	#Returns True for all cases listed here, independent of case (upper/lower/mixed)
	#Returns False for all other cases.
	#
	# example:
	#  print parsebool("y")
	#	 True
	#  print parsebool("test")
	#	 False
	#
	#######################################################
	"""
	return val.lower() in ("yes", "true", "t", "1", "tru", "tr", "y")


if __name__ == "__main__":
    #check number of arguments
    if (len(sys.argv) != 2): 
        print "reduction code requires a datatext file"
        sys.exit()
    if not(os.path.isfile(sys.argv[1])):
        print "data text file ", sys.argv[1], " not found"
        sys.exit()
    dgsreduction(XMLfile=sys.argv[1])
