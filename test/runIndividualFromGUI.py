'''
This is provided for testing purposes. It allows you to run the single system test
within Mantid environent
It can be useful for debugging because the errors do not alway 'get out' of
the sub-process used for running the tests in the regular way
'''
from mantidsimple import *
import sys


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
#sys.path.append(r'C:\Backup\Backup_folder1\work\code\Mantid\git\mantid\Test\systemtests\Data\SANS2D')
#sys.path.append(r'C:\Backup\Backup_folder1\work\code\Mantid\git\mantid\Test\systemtests\Data\LOQ')
# Find these first
sys.path.insert(0,tests_dir)
modlToRun = ['ISISDirectInelastic'];
testToRun = ['MARIReductionFromWorkspace'];


module = __import__(modlToRun[0])
reload(module)
#reload(sys.modules['isis_reduction_steps'])

className = testToRun[0]
if ( len(testToRun) > 1 ) : className = testToRun[1]
try:
    testClass = getattr(module, className)()
    testClass.execute()
except:
    os.chdir(this_dir)
    raise

outcome = testClass.doValidation()
print 'Test result: ' + str(outcome)
os.chdir(this_dir)
