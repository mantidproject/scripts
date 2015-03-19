import unittest
import ncs

import mantid
from mantid.simpleapi import LoadVesuvio,CropWorkspace, mtd

class CsHS204Test(unittest.TestCase):
    
    def setUp(self):
        runs = "15039-15045"
        spectra = "135"
        diff_type="SingleDifference" # Allowed values=Single,Double,Thick
        ip_file = "IP0004_10.par"
        
        raw_ws = LoadVesuvio(Filename=runs, SpectrumList=spectra,
                             Mode=diff_type,InstrumentParFile=ip_file)
        raw_ws = CropWorkspace(raw_ws,XMin=50.0,XMax=562.0)
    
    def tearDown(self):
        mtd.clear()

    def test_k_fixed_no_background(self):
        fit_options = self._setup_fit_options()
        raw_ws = mtd['raw_ws']
        raw_ws = ncs.preprocess(raw_ws, fit_options)
        reduced_chi_square, params_ws = ncs.run_fit(raw_ws, fit_options)
        message = ncs.display_fit_output(reduced_chi_square, params_ws,fit_options)

        expected = [\
            {'Name': 'f0.Width', 'Value': 3.9527823429702624, 'Error': 0.18552724464340611},\
            {'Name': 'f0.FSECoeff', 'Value': 0.465839866544787, 'Error': 0.0},\
            {'Name': 'f0.C_0', 'Value': 21.228162833736132, 'Error': 0.7484061245669142},
            {'Name': 'f1.Width', 'Value': 10.0, 'Error': 0.0},\
            {'Name': 'f1.Intensity', 'Value': 4.31589459135, 'Error': 0.426413174495},\
            {'Name': 'f2.Width', 'Value': 13.0, 'Error': 0.0},\
            {'Name': 'f2.Intensity', 'Value': 2.92303015759, 'Error': 0.463509032992},\
            {'Name': 'f3.Width', 'Value': 30.0, 'Error': 0.0},\
            {'Name': 'f3.Intensity', 'Value': 0.931308549279, 'Error': 0.22624298123},\
            {'Name': 'Cost function value', 'Value': 3.28100193547, 'Error': 0.0}\
        ]
        self._check_params(params_ws,expected)
        
        expected_msg = \
        """

Reduced Chi-Square =3.281002

Fitting in the TOF space
--------------------------------------------------------------------------------
The mass M(1)=1.007900

Parameters values in the Y space:

See results log for resolution parameters w_L1, w_L0 etc
St. dev. of momentum distr. = 3.952782 +/- 0.185527
Scatt. int. (area, normalised) = 0.722086 +/- 0.004042
Hermite polynomial expansion coefficient c0 = 1.000000 +/- 0.000000
Hermite polynomial expansion coefficient c2 = 0.000000 +/- 0.000000
Hermite polynomial expansion coefficient c4 = 0.000000 +/- 0.000000
Hermite polynomial expansion coefficient a0 = c0/(2^0*0!) = 1.000000 +/ 0.000000
Hermite polynomial expansion coefficient a2 = c2/(2^2*1!) = 0.000000 +/ 0.000000
Hermite polynomial expansion coefficient a4 = c4/(2^4*2!) = 0.000000 +/ 0.000000

FSE coefficient k by the k/q He_3(y) expansion member = 0.465840 +/- 0.021865

The coefficient k calculated in a harmonic oscillator model would be k = sigma*sqrt(2)/12 = 0.465840

--------------------------------------------------------------------------------
The mass M(2)=16.000000

Parameters values in the Y space:

See results log for resolution parameters w_L1, w_L0 etc
St. dev. of momentum distr. = 10.000000 +/- 0.000000
Scatt. int. (area, normalised) = 0.146807 +/- 0.010727
--------------------------------------------------------------------------------
The mass M(3)=27.000000

Parameters values in the Y space:

See results log for resolution parameters w_L1, w_L0 etc
St. dev. of momentum distr. = 13.000000 +/- 0.000000
Scatt. int. (area, normalised) = 0.099428 +/- 0.010590
--------------------------------------------------------------------------------
The mass M(4)=133.000000

Parameters values in the Y space:

See results log for resolution parameters w_L1, w_L0 etc
St. dev. of momentum distr. = 30.000000 +/- 0.000000
Scatt. int. (area, normalised) = 0.031679 +/- 0.010457
"""
        self.assertEquals(expected_msg, message)

#=============================================================================================================

    def test_k_fixed_with_background(self):
        fit_options = self._setup_fit_options()
        fit_options.background_order = 2

        raw_ws = mtd['raw_ws']
        raw_ws = ncs.preprocess(raw_ws, fit_options)
        reduced_chi_square, params_ws = ncs.run_fit(raw_ws, fit_options)
        ncs.display_fit_output(reduced_chi_square, params_ws,fit_options)
        
        expected = [\
            {'Name': 'f0.Width', 'Value': 4.11688961577, 'Error': 0.243634042057},\
            {'Name': 'f0.FSECoeff', 'Value': 0.485180094118, 'Error': 0.0},\
            {'Name': 'f0.C_0', 'Value': 20.0966608709, 'Error': 1.21251471312},\
            {'Name': 'f1.Width', 'Value': 10.0, 'Error': 0.0},\
            {'Name': 'f1.Intensity', 'Value': 4.09470837231, 'Error': 0.433950802885},\
            {'Name': 'f2.Width', 'Value': 13.0, 'Error': 0.0},\
            {'Name': 'f2.Intensity', 'Value': 2.93153580371, 'Error': 0.464917511735},\
            {'Name': 'f3.Width', 'Value': 30.0, 'Error': 0.0},\
            {'Name': 'f3.Intensity', 'Value': 0.969159274084, 'Error': 0.23026233512},\
            {'Name': 'f4.A0', 'Value': -0.00656341362348, 'Error': 0.00289823424173},\
            {'Name': 'f4.A1', 'Value': 18.9776869966, 'Error':19.1127629958},\
            {'Name': 'f4.A2', 'Value': 93.5252719638, 'Error': 29084.614777},\
            {'Name': 'Cost function value', 'Value': 3.38368838556, 'Error': 0.0}\
        ]
        self._check_params(params_ws,expected)

    def test_k_free_no_background(self):
        fit_options = self._setup_fit_options()
        #------------------- FIRST MASS ------------------------------------
        fit_options.masses[0]['k_free'] = True
        #-------------------------------------------------------------------
        raw_ws = mtd['raw_ws']
        raw_ws = ncs.preprocess(raw_ws, fit_options)
        reduced_chi_square, params_ws = ncs.run_fit(raw_ws, fit_options)
        message = ncs.display_fit_output(reduced_chi_square, params_ws,fit_options)

        expected = [\
            {'Name': 'f0.Width', 'Value': 4.23594841483436, 'Error': 0.21738339873184126},\
            {'Name': 'f0.FSECoeff', 'Value': -76.53672476106186, 'Error': 20.37825852906559},\
            {'Name': 'f0.C_0', 'Value': 22.602517742831605, 'Error': 0.8424355253357285},\
            {'Name': 'f1.Width', 'Value': 10.0, 'Error': 0.0},\
            {'Name': 'f1.Intensity', 'Value': 4.315249149863499, 'Error': 0.42641324871786956},\
            {'Name': 'f2.Width', 'Value': 13.0, 'Error': 0.0},\
            {'Name': 'f2.Intensity', 'Value': 2.9233511235166656, 'Error': 0.463509049047},\
            {'Name': 'f3.Width', 'Value': 30.0, 'Error': 0.0},\
            {'Name': 'f3.Intensity', 'Value': 0.9310641107710553, 'Error': 0.226242997134},\
            {'Name': 'Cost function value', 'Value': 3.2673708431428428, 'Error': 0.0}\
        ]
        self._check_params(params_ws,expected)
        
        expected_msg = \
"""

Reduced Chi-Square =3.267371

Fitting in the TOF space
--------------------------------------------------------------------------------
The mass M(1)=1.007900

Parameters values in the Y space:

See results log for resolution parameters w_L1, w_L0 etc
St. dev. of momentum distr. = 4.235948 +/- 0.217383
Scatt. int. (area, normalised) = 0.734511 +/- 0.003924
Hermite polynomial expansion coefficient c0 = 1.000000 +/- 0.000000
Hermite polynomial expansion coefficient c2 = 0.000000 +/- 0.000000
Hermite polynomial expansion coefficient c4 = 0.000000 +/- 0.000000
Hermite polynomial expansion coefficient a0 = c0/(2^0*0!) = 1.000000 +/ 0.000000
Hermite polynomial expansion coefficient a2 = c2/(2^2*1!) = 0.000000 +/ 0.000000
Hermite polynomial expansion coefficient a4 = c4/(2^4*2!) = 0.000000 +/ 0.000000

FSE coefficient k by the k/q He_3(y) expansion member = -3.386204 +/- 1.027802

The coefficient k calculated in a harmonic oscillator model would be k = sigma*sqrt(2)/12 = 0.499211

--------------------------------------------------------------------------------
The mass M(2)=16.000000

Parameters values in the Y space:

See results log for resolution parameters w_L1, w_L0 etc
St. dev. of momentum distr. = 10.000000 +/- 0.000000
Scatt. int. (area, normalised) = 0.140232 +/- 0.011058
--------------------------------------------------------------------------------
The mass M(3)=27.000000

Parameters values in the Y space:

See results log for resolution parameters w_L1, w_L0 etc
St. dev. of momentum distr. = 13.000000 +/- 0.000000
Scatt. int. (area, normalised) = 0.095000 +/- 0.010934
--------------------------------------------------------------------------------
The mass M(4)=133.000000

Parameters values in the Y space:

See results log for resolution parameters w_L1, w_L0 etc
St. dev. of momentum distr. = 30.000000 +/- 0.000000
Scatt. int. (area, normalised) = 0.030257 +/- 0.010815
"""
        self.assertEquals(expected_msg,message)
#=============================================================================================================

    def test_k_free_with_background(self):
        fit_options = self._setup_fit_options()
        #------------------- FIRST MASS --------------------------
        fit_options.masses[0]['k_free'] = True
        #-------------------------------------------------------------------
        fit_options.background_order = 2
        raw_ws = mtd['raw_ws']
        raw_ws = ncs.preprocess(raw_ws, fit_options)
        reduced_chi_square, params_ws = ncs.run_fit(raw_ws, fit_options)
        ncs.display_fit_output(reduced_chi_square, params_ws,fit_options)

        expected = [\
            {'Name': 'f0.Width', 'Value': 4.12370499357, 'Error': 0.280633041767},
            {'Name': 'f0.FSECoeff', 'Value': -32.2147778697, 'Error': 21.4446620068},\
            {'Name': 'f0.C_0', 'Value': 20.6729781391, 'Error': 1.46423308815},\
            {'Name': 'f1.Width', 'Value': 10.0, 'Error': 0.0},\
            {'Name': 'f1.Intensity', 'Value': 4.11299150262, 'Error': 0.434739226944},\
            {'Name': 'f2.Width', 'Value': 13.0, 'Error': 0.0},\
            {'Name': 'f2.Intensity', 'Value': 2.93613295977, 'Error': 0.465065906128},\
            {'Name': 'f3.Width', 'Value': 30.0, 'Error': 0.0},\
            {'Name': 'f3.Intensity', 'Value': 0.963632631744, 'Error': 0.230636343676},\
            {'Name': 'f4.A0', 'Value': -0.00518193477871, 'Error': 0.00300080764282},\
            {'Name': 'f4.A1', 'Value': 15.8505498861, 'Error': 19.129998072},\
            {'Name': 'f4.A2', 'Value': -30.6569391098, 'Error': 29179.377345725108},\
            {'Name': 'Cost function value', 'Value': 3.34158168582, 'Error': 0.0}\
        ]
        self._check_params(params_ws,expected)

    def _setup_fit_options(self):
        fit_options = ncs.FitOptions()
        fit_options.workspace_index = 0
        fit_options.smooth_points = None
        fit_options.bad_data_error = 1e6
        ## Background Polynomial (none)
        fit_options.background_order = None

        mass1 = {'value':1.0079, 'widths':[2,5,7], 'function':'GramCharlier', 
                  'hermite_coeffs':[1,0,0],'k_free':False, 'sears_flag':1}
        mass2 = {'value':16.0, 'widths':10, 'function':'Gaussian', }
        mass3 = {'value':27.0, 'widths':13, 'function':'Gaussian', }
        mass4 = {'value':133.0, 'widths':30, 'function':'Gaussian'}
        fit_options.masses = [mass1, mass2, mass3, mass4]
        
        ## Intensity constraints
        fit_options.constraints = ([0,1,0,-4])
        
        return fit_options


    def _check_params(self, params_ws, expected_rows):
         
        for index, expected_row in enumerate(expected_rows):
            cur_param_row = params_ws.row(index)
            self.assertEquals(expected_row['Name'],cur_param_row['Name'])
            self.assertAlmostEqual(expected_row['Value'],cur_param_row['Value'], places=8)
            self.assertAlmostEqual(expected_row['Error'],cur_param_row['Error'], places=8)

##########################################################

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(CsHS204Test))
    runner = unittest.TextTestRunner()
    runner.failfast = False
    #Run using either runner
    res = runner.run(suite)
    if len(res.errors) != 0 or len(res.failures) != 0:
        raise RuntimeError("Tests failed. See results log")
