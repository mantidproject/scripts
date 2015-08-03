"""These are more integration tests as they will require that the test data is available
and that mantid can be imported
"""
import unittest

from mantid.api import (WorkspaceGroup, MatrixWorkspace)
import mantid
from vesuvio.workflow import fit_tof

class FitTofTest(unittest.TestCase):

    def test_fit_single_spectrum_tof_and_no_background(self):
        flags = self._create_test_flags(background=False)
        runs = "15039-15045"

        fit_results = fit_tof(runs, flags)
        self.assertTrue(isinstance(fit_results, tuple))
        self.assertEqual(4, len(fit_results))

        fitted_wsg = fit_results[0]
        self.assertTrue(isinstance(fitted_wsg, WorkspaceGroup))
        self.assertEqual(2, len(fitted_wsg))

        fitted_ws = fitted_wsg[0]
        self.assertTrue(isinstance(fitted_ws, MatrixWorkspace))
        self.assertEqual(7, fitted_ws.getNumberHistograms())
        # self.assertAlmostEqual(50.0, fitted_ws.readX(0)[0])
        # self.assertAlmostEqual(562.0, fitted_ws.readX(0)[-1])

        # self.assertAlmostEqual(0.000928695463881635, fitted_ws.readY(0)[0])
        # self.assertAlmostEqual(0.00722948549525415, fitted_ws.readY(0)[-1])
        # self.assertAlmostEqual(1.45746507977816e-05, fitted_ws.readY(1)[0])
        # self.assertAlmostEqual(7.33791942084561e-05, fitted_ws.readY(1)[-1])

        fitted_params = fit_results[1]
        self.assertTrue(isinstance(fitted_params, MatrixWorkspace))
        self.assertEqual(10, fitted_params.getNumberHistograms())

        chisq_values = fit_results[2]
        self.assertTrue(isinstance(chisq_values, list))
        self.assertEqual(1, len(chisq_values))

        exit_iteration = fit_results[3]
        self.assertTrue(isinstance(exit_iteration, int))


    def test_fit_single_spectrum_tof_including_background(self):
        flags = self._create_test_flags(background=True)
        runs = "15039-15045"

        fit_results = fit_tof(runs, flags)
        self.assertTrue(isinstance(fit_results, tuple))
        self.assertEqual(4, len(fit_results))

        fitted_wsg = fit_results[0]
        self.assertTrue(isinstance(fitted_wsg, WorkspaceGroup))
        self.assertEqual(2, len(fitted_wsg))

        fitted_ws = fitted_wsg[0]
        self.assertTrue(isinstance(fitted_ws, MatrixWorkspace))
        self.assertEqual(8, fitted_ws.getNumberHistograms())
        # self.assertAlmostEqual(50.0, fitted_ws.readX(0)[0])
        # self.assertAlmostEqual(562.0, fitted_ws.readX(0)[-1])

        # self.assertAlmostEqual(0.000928695463881635, fitted_ws.readY(0)[0])
        # self.assertAlmostEqual(0.00722948549525415, fitted_ws.readY(0)[-1])
        # self.assertAlmostEqual(-0.00756178413274695, fitted_ws.readY(1)[0])
        # self.assertAlmostEqual(0.00355843687365601, fitted_ws.readY(1)[-1])

        fitted_params = fit_results[1]
        self.assertTrue(isinstance(fitted_params, MatrixWorkspace))
        self.assertEqual(14, fitted_params.getNumberHistograms())

        chisq_values = fit_results[2]
        self.assertTrue(isinstance(chisq_values, list))
        self.assertEqual(1, len(chisq_values))

        exit_iteration = fit_results[3]
        self.assertTrue(isinstance(exit_iteration, int))


    def test_fit_bank_by_bank_for_forward_spectra_tof_and_no_background(self):
        flags = self._create_test_flags(background=False)
        flags['fit_mode'] = 'bank'
        flags['spectra'] = 'forward'
        runs = "15039-15045"

        fit_results = fit_tof(runs, flags)
        self.assertTrue(isinstance(fit_results, tuple))
        self.assertEquals(4, len(fit_results))

        fitted_banks = fit_results[0]
        self.assertTrue(isinstance(fitted_banks, list))
        self.assertEqual(8, len(fitted_banks))

        bank1 = fitted_banks[0]
        self.assertTrue(isinstance(bank1, WorkspaceGroup))

        bank1_data = bank1[0]
        self.assertTrue(isinstance(bank1_data, MatrixWorkspace))

        # self.assertAlmostEqual(50.0, bank1_data.readX(0)[0])
        # self.assertAlmostEqual(562.0, bank1_data.readX(0)[-1])

        # self.assertAlmostEqual(0.000107272755986595, bank1_data.readY(1)[0])
        # self.assertAlmostEqual(0.000585633970072128, bank1_data.readY(1)[-1])

        bank8 = fitted_banks[-1]
        self.assertTrue(isinstance(bank8, WorkspaceGroup))

        bank8_data = bank8[0]
        self.assertTrue(isinstance(bank8_data, MatrixWorkspace))

        # self.assertAlmostEqual(50.0, bank8_data.readX(0)[0])
        # self.assertAlmostEqual(562.0, bank8_data.readX(0)[-1])

        # self.assertAlmostEqual(0.000596850729120898, bank8_data.readY(1)[0])
        # self.assertAlmostEqual(0.000529343513813141, bank8_data.readY(1)[-1])

        chisq_values = fit_results[2]
        self.assertTrue(isinstance(chisq_values, list))
        self.assertEqual(8, len(chisq_values))

        exit_iteration = fit_results[3]
        self.assertTrue(isinstance(exit_iteration, int))


    def test_fit_spectra_by_spectra_for_forward_spectra_tof_and_no_background(self):
        flags = self._create_test_flags(background=False)
        flags['fit_mode'] = 'spectra'
        flags['spectra'] = '143-144'
        runs = "15039-15045"

        fit_results = fit_tof(runs, flags)
        self.assertTrue(isinstance(fit_results, tuple))
        self.assertEquals(4, len(fit_results))

        fitted_spec = fit_results[0]
        self.assertTrue(isinstance(fitted_spec, list))
        self.assertEqual(2, len(fitted_spec))

        spec143 = fitted_spec[0]
        self.assertTrue(isinstance(spec143, WorkspaceGroup))

        spec143_data = spec143[0]
        self.assertTrue(isinstance(spec143_data, MatrixWorkspace))

        # self.assertAlmostEqual(50.0, spec143_data.readX(0)[0])
        # self.assertAlmostEqual(562.0, spec143_data.readX(0)[-1])

        # self.assertAlmostEqual(2.37897941103748e-06, spec143_data.readY(1)[0])
        # self.assertAlmostEqual(3.58226563303213e-05, spec143_data.readY(1)[-1])

        spec144 = fitted_spec[-1]
        self.assertTrue(isinstance(spec144, WorkspaceGroup))

        spec144_data = spec144[0]
        self.assertTrue(isinstance(spec144_data, MatrixWorkspace))

        # self.assertAlmostEqual(50.0, spec144_data.readX(0)[0])
        # self.assertAlmostEqual(562.0, spec144_data.readX(0)[-1])

        # self.assertAlmostEqual(5.57952304659615e-06, spec144_data.readY(1)[0])
        # self.assertAlmostEqual(6.00056973529846e-05, spec144_data.readY(1)[-1])

        chisq_values = fit_results[2]
        self.assertTrue(isinstance(chisq_values, list))
        self.assertEqual(2, len(chisq_values))

        exit_iteration = fit_results[3]
        self.assertTrue(isinstance(exit_iteration, int))


    def _create_test_flags(self, background):
        runs = "15039-15045"
        flags = dict()
        flags['fit_mode'] = 'spectrum'
        flags['spectra'] = '135'

        mass1 = {'value': 1.0079, 'function': 'GramCharlier', 'width': [2, 5, 7],
                  'hermite_coeffs': [1,0,0], 'k_free': 0, 'sears_flag': 1}
        mass2 = {'value': 16.0, 'function': 'Gaussian', 'width': 10}
        mass3 = {'value': 27.0, 'function': 'Gaussian', 'width': 13}
        mass4 = {'value': 133.0, 'function': 'Gaussian', 'width': 30}
        flags['masses'] = [mass1, mass2, mass3, mass4]
        flags['intensity_constraints'] = [0, 1, 0, -4]
        if background:
            flags['background'] = {'function': 'Polynomial', 'order':3}
        else:
            flags['background'] = None
        flags['ip_file'] = 'IP0004_10.par'
        flags['diff_mode'] = 'single'
        flags['gamma_correct'] = True
        flags['ms_flags'] = dict()
        flags['ms_flags']['SampleWidth'] = 10.0
        flags['ms_flags']['SampleHeight'] = 10.0
        flags['ms_flags']['SampleDepth'] = 0.5
        flags['ms_flags']['SampleDensity'] = 241
        flags['fit_minimizer'] = 'Levenberg-Marquardt,AbsError=1e-08,RelError=1e-08'

        return flags


if __name__ == '__main__':
    unittest.main()
