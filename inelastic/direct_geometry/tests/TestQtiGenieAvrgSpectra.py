from mantid.simpleapi import *
from qtiGenie import avrg_spectra
import math
import numpy
import unittest

# 
class TestQtiGenieAvrgSpectra(unittest.TestCase):

    def setUp(self):        
        dataX = [-1,0,1,2,3,-1,0,1,2,3,-1,0,1,2,3,-1,0,1,2,3]
        dataY = [2,2,2,2,3,3,3,3,1,1,1,1,10,10,10,10]
        dataE = [None]*(len(dataY))
        for i in range(0,len(dataY)):
            dataE[i] = math.sqrt(dataY[i])
            self.ws="TestWS"
        
        self.dataX = dataX;
        CreateWorkspace(OutputWorkspace=self.ws,DataX=dataX,DataY=dataY,DataE=dataE, NSpec=4,UnitX="DeltaE")

    def tearDown(self):
        DeleteWorkspace(self.ws)
        pass


    def test_avrgN_WSName(self):
        sums,stats=avrg_spectra(self.ws)
        # single spectra with specified number of elements
        self.assertEqual(len(sums),4)
        self.assertEqual(stats,[4,0,0])
        for i in range(0,4):
            self.assertEqual(sums[i],(2+3+1+10)/4)

    def test_avrgN_WSPtr(self):
    # should work for workspace handle too
        ws = mtd[self.ws]
        sums,stats=avrg_spectra(ws)
        # single spectra with specified number of elements
        self.assertEqual(len(sums),4)
        self.assertEqual(stats,[4,0,0])
        for i in range(0,4):
            self.assertEqual(sums[i],(2+3+1+10)/4)
            
    def test_avrg3Sigma(self):
        sums,stats=avrg_spectra(self.ws,0,3,True)
        # single spectra with specified number of elements
        self.assertEqual(len(sums),4)
        self.assertEqual(stats,[4,0,0])     
        for i in range(0,stats[0]):
            self.assertAlmostEqual(sums[i],(4./(1./2+1./3+1.+0.1)))

    def test_singleValSigma(self):
        rwsName="ReducedWS"        
        Integration(InputWorkspace=self.ws,OutputWorkspace=rwsName,RangeLower=str(self.dataX[0]),RangeUpper=str(self.dataX[4]),IncludePartialBins='1')
        sums,stats=avrg_spectra(rwsName,0,3,True)
        # single spectra with slingle element
        #self.assertEqual(len(sums),1)
        self.assertEqual(stats,[4,0,0])     
        self.assertAlmostEqual(sums,stats[0]*(4./(1./2+1./3+1.+0.1)))
        DeleteWorkspace(rwsName)  

    def test_singleValAvrg(self):
        rwsName="ReducedWS"        
        Integration(InputWorkspace=self.ws,OutputWorkspace=rwsName,RangeLower="-1",RangeUpper="3",IncludePartialBins='1')
        sums,stats=avrg_spectra(rwsName)
        # single spectra with slingle element
        #self.assertEqual(len(sums),1)
        self.assertEqual(stats,[4,0,0])  
        # result should be normal average of the worskspace
        self.assertAlmostEqual(sums,(self.dataX[4]-self.dataX[0])*(2+3+1+10)/stats[0])
        DeleteWorkspace(rwsName)  


if __name__ == '__main__':
    unittest.main()