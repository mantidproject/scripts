'''
This is provided for testing purposes. It allows you to run the single system test
within Mantid environent
It can be useful for debugging because the errors do not alway 'get out' of
the sub-process used for running the tests in the regular way
'''
#import mantidsimple
#print mantidsimple.__file__
from mantid.simpleapi import *
import sys
import os
import inspect


this_dir = sys.path[0]
print "Script resides in : ",this_dir
os.chdir(r'd:\Data\MantidSystemTests\SystemTests')
stressmodule_dir = r'd:\Data\MantidSystemTests\StressTestFramework'
tests_dir = r'd:/Data/MantidSystemTests/SystemTests/AnalysisTests'

#reference_dir = 'C:/Backup/Backup_folder1/work/code/Mantid/git/mantid/Test/systemTests/SystemTests/AnalysisTests/ReferenceResults'
#mtdpy_header_dir = 'C:/Backup/Backup_folder1/work/code/Mantid/mantid/Code/Mantid/debug'

sys.path.append(r'd:/Data/Mantid_GIT/Code/win64/bin/Release')
#sys.path.append(r'C:\Backup\Backup_folder1\work\code\Mantid\builds\all\bin\Debug')
sys.path.append(stressmodule_dir)
sys.path.append(r'd:\Data\MantidSystemTests\SystemTests\Data\SANS2D')
sys.path.append(r'd:\Data\MantidSystemTests\SystemTests\Data\LOQ')
config.setDataSearchDirs('d:/Data/MantidSystemTests/SystemTests/Data')
# Find these first
sys.path.insert(0,tests_dir)
#Failed: 

modlToRun = ['Diffraction_Workflow_Test']
# My: 'SXDAnalysis','WishDiffuseScattering','Diffraction_Workflow_Test',


# not migrated 'EQSANSLive' EQSANSEff EQSANSIQOutput EQSANSTrans

print "-------------------------------------------------------------------------------"


for i in xrange(0,len(modlToRun)):
    print "-------------------------------------------------------------------------------"        
    print "---> Testing modue "+modlToRun[i]
    print "-------------------------------------------------------------------------------"        
    
    module = __import__(modlToRun[i])
    reload(module)
    
    testsToRun = []
    
    for name, obj in inspect.getmembers(module):
        if inspect.isclass(obj) and obj.__module__ == modlToRun[i] and not inspect.isabstract(obj) :
           #print "Testing class: "+modlToRun[i]+" Name: "+name         
            testsToRun.append(name)

            
    if len(testsToRun)==0:
        testsToRun.append(modlToRun[i]);
        #print "--->  Testing class:  ",testsToRun[0]


    for j in xrange(0,len(testsToRun)):
        print "--->  Testing class:  ",testsToRun[j]    
        with open("CompletedTests.txt", "a") as myfile:
                myfile.write("Started: "+modlToRun[i]+"\t\t Test: "+testsToRun[j]+"\t finished:\t ")
        
        
        className = testsToRun[j]
        try:
            testClass = getattr(module, className)()
            testClass.execute()
        except:
            with open("CompletedTests.txt", "a") as myfile:
                myfile.write("Test: -----------------------------\n")
            #continue
            os.chdir(this_dir)
            raise
	    
        outcome = testClass.doValidation()    
        print 'Test result: ' + str(outcome)
        
        with open("CompletedTests.txt", "a") as myfile:
                myfile.write("Test: "+testsToRun[j]+"\n")

        
    print "-------------------------------------------------------------------------------"        
    
os.chdir(this_dir)
