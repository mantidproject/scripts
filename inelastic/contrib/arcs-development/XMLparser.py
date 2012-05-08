'''
Created on Feb 20, 2012

@author: andrei
'''
from lxml import etree
#from string import lower
from copy import deepcopy
import os
import imp
from mantid.kernel import logger

class XMLparser(object):
    '''
    classdocs
    '''

    def __init__(self,filename):
        '''
        Constructor
        '''
        if filename!=None:
            doc=etree.parse(filename)
            root=doc.getroot()
            if root.tag not in ['dgsreduction']:
                raise RuntimeError("This is not a dgsreduction xml file")
            self.root=root
            defaults=root.find('defaults')
            if defaults==None:
                raise RuntimeError("The defaults section is missing")
            temp={}
            self.recursiveGetAttributes(defaults,temp)
            #First verify that there is an instrument in the tempdefaultdict.
            if not temp.has_key('instrument'):
                raise RuntimeError("The instrument is not defined")
            if temp['instrument'] not in ['ARCS','SEQUOIA','HYSPEC','CNCS']:
                raise ValueError("Unknown instrument "+str(temp['instrument']))    
            else:
                #if there is an instrument value, then try to load the
                #defaults for that instrument.
                self.loadinstrumentdict(temp['instrument'])

            self.defaultdict.update(temp)
            self.calibdict=deepcopy(self.defaultdict)
            calibration=root.find('calibration')
            self.recursiveGetAttributes(calibration,self.calibdict)
            self.datadicts=[]
            listruns = root.findall('scan')
        
            for i in range(len(listruns)):
                tempi=deepcopy(self.defaultdict)
                self.recursiveGetAttributes(listruns[i],tempi)
                tempi['type']=listruns[i].tag
                if tempi['runs']!=[]:
                    self.datadicts.append(tempi)
                
    def recursiveGetAttributes(self,node,dictionary):
        if node.tag!=etree.Comment:
            if node.tag=='mask':
                if not dictionary.has_key('mask'):
                    dictionary['mask']=[]
                if node.attrib!={}:        
                    dictionary['mask'].append(node.attrib)
            else:
                for k in node.attrib.keys():
                    dictionary[k]=self.parseto(k,node.attrib[k])
                if (node.text!=None) and (len(node.text.strip())!=0):
                    dictionary[node.tag]=self.parseto(node.tag,node.text)
                for c in node.iterchildren():
                    self.recursiveGetAttributes(c,dictionary)

    def parseto(self,tag,value):
        if tag == 'scantype':
            tstring = value[0:2].lower()
            if tstring == 'si':
                return 'single'
            elif tstring == 'st':
                return 'step'
            elif tstring == 'sw':
                return 'sweep'
            else:
                raise ValueError("scantype not understood (single, step or sweep).")
                

        if tag in ['t0','efixed','emin','emax','ebin','qstep','tibgstart','tibgstop','vanmin','vanmax',
                   'logvaluemin','logvaluemax','logvaluestep','powderanglestep']:
            return float(value)

        if tag in ['runs','vanruns','pixel','tube','bank']:
            #create empty list to put the runs in
            runs = []
            #split the commas
            parts = value.split(',')
            #now deal with the hyphens
            for p in parts:
                if len(p) > 0:
                    elem = p.split("-")
                if len(elem) == 1:
                    runs.append(int(elem[0]))
                if len(elem) == 2:
                    startelem = int(elem[0])
                    endelem   = int(elem[1])
                    if endelem < startelem:
                        raise ValueError("The element after the hyphen needs to be greater or equal than the first element")
                    elemlist  = range(startelem,endelem+1)
                    runs.extend(elemlist)
            return runs
        if tag in ['calce','kiokf','tibg','normalizedcalibration','filterbadpulses']:
            return (value.lower() in ("yes", "true", "t", "1", "tru", "tr", "y"))
        if tag=='units':
            #replace all whitespace with nothing
            value.replace(' ', '')
            #this is a list of strings of the supported filetypes
            unittypesexist = ['DeltaE', 'DeltaE_inWavenumber', 'Energy', 'Energy_inWavenumber',
                              'Momentum', 'MomentumTransfer', 'QSquared', 'TOF', 'Wavelength',
                              'dspacing']
            #get a string all lowercase.
            unitlower = []
            for u in unittypesexist:
                unitlower.append(u.lower())
            #Find the matching value of the unit list.
            i = unitlower.index(value.lower())
            return unittypesexist[i]
        if tag=='save':
            #replace all whitespace with commas
            value.replace(' ', ',')
            #split up parts by comma delimited values
            parts = value.split(',')
            #this is a list of strings of the supported filetypes
            filetypesexist = ['phx', 'spe', 'nxspe', 'par', 'jpg', 'nxs', 'iofq','iofe','summary']
            #loop over the elements and fill up a list of the
            #filetypes to return.
            filestoprocess = []
            for p in parts:
                p = p.strip(' "')
                if p.lower() in filetypesexist:
                    filestoprocess.append(p.lower())
                else:
                    if len(p) > 0:
                        logger.warning("*****FILETYPE "+p+" DOES NOT EXIST. This file will be ignored.*****")
            return filestoprocess
        if tag=='friendlynamelogs':
            #convert to an arry of strings, split by the commas
            parts = value.replace(' ', '').split(',')
            return parts
        if tag in ['goniometermotor','goniometermotoraxis','goniometermotoroffset','goniometermotordirection']:
            try:
                return eval(value)
            except:
                return value            
        return value


    def loadinstrumentdict(self,instrumentname):
        """
        instrumentname is a string
        returns a dictionary of instrument defaults
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
            raise RuntimeError("Could not find instrument definition module "+os.path.abspath(os.curdir)+"/"+modulename+'.py')

        self.defaultdict=m.instrumentparameters()




if __name__ == "__main__":
    filename='test1b.xml'
    parsed=XMLparser(filename)
    
    print "-----------default--------------"
    print parsed.defaultdict 
    print "-----------calibration--------------"
    print parsed.calibdict
    print "-----------data--------------"
    print parsed.datadicts
