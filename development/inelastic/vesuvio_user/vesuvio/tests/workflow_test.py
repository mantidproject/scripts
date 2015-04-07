"""These are more integration tests as they will require that the test data is available
and that mantid can be imported
"""
import unittest

import mantid
from vesuvio.workflow import fit_tof

class FitTofTest(unittest.TestCase):

    def test_fit_single_spectrum_tof_and_no_background(self):
        flags = self._create_test_flags(background=False)
        runs = "15039-15045"

        fit_results = fit_tof(runs, flags)
        self.assertTrue(isinstance(fit_results, mantid.api.WorkspaceGroup))
        self.assertEqual(2, fit_results.size())

        fitted_ws = fit_results[0]
        self.assertEqual(7, fitted_ws.getNumberHistograms())
        self.assertAlmostEqual(50.0, fitted_ws.readX(0)[0])
        self.assertAlmostEqual(562.0, fitted_ws.readX(0)[-1])

        self.assertAlmostEqual(0.000928695463881635, fitted_ws.readY(0)[0])
        self.assertAlmostEqual(0.00722948549525415, fitted_ws.readY(0)[-1])
        self.assertAlmostEqual(1.45746507977816e-05, fitted_ws.readY(1)[0])
        self.assertAlmostEqual(7.33791942084561e-05, fitted_ws.readY(1)[-1])

        fitted_params = fit_results[1]
        self.assertEqual(10, fitted_params.rowCount())

    def test_fit_single_spectrum_tof_including_background(self):
        flags = self._create_test_flags(background=True)
        runs = "15039-15045"

        fit_results = fit_tof(runs, flags)
        self.assertTrue(isinstance(fit_results, mantid.api.WorkspaceGroup))

        fitted_ws = fit_results[0]
        self.assertEqual(8, fitted_ws.getNumberHistograms())
        self.assertAlmostEqual(50.0, fitted_ws.readX(0)[0])
        self.assertAlmostEqual(562.0, fitted_ws.readX(0)[-1])

        self.assertAlmostEqual(0.000928695463881635, fitted_ws.readY(0)[0])
        self.assertAlmostEqual(0.00722948549525415, fitted_ws.readY(0)[-1])
        self.assertAlmostEqual(-0.00756178413274695, fitted_ws.readY(1)[0])
        self.assertAlmostEqual(0.00355843687365601, fitted_ws.readY(1)[-1])

        fitted_params = fit_results[1]
        self.assertEqual(14, fitted_params.rowCount())

    def test_fit_bank_by_bank_for_forward_spectra_tof_and_no_background(self):
        flags = self._create_test_flags(background=False)
        flags['fit_mode'] = 'bank'
        flags['spectra'] = 'forward'
        runs = "15039-15045"

        fit_results = fit_tof(runs, flags)
        self.assertTrue(isinstance(fit_results, list))
        self.assertEquals(8, len(fit_results))

        bank1 = fit_results[0]
        bank1_data = bank1[0]
        self.assertAlmostEqual(50.0, bank1_data.readX(0)[0])
        self.assertAlmostEqual(562.0, bank1_data.readX(0)[-1])

        self.assertAlmostEqual(0.000107272755986595, bank1_data.readY(1)[0])
        self.assertAlmostEqual(0.000585633970072128, bank1_data.readY(1)[-1])

        bank8 = fit_results[-1]
        bank8_data = bank8[0]
        self.assertAlmostEqual(50.0, bank8_data.readX(0)[0])
        self.assertAlmostEqual(562.0, bank8_data.readX(0)[-1])

        self.assertAlmostEqual(0.000596850729120898, bank8_data.readY(1)[0])
        self.assertAlmostEqual(0.000529343513813141, bank8_data.readY(1)[-1])

    def test_fit_spectra_by_spectra_for_forward_spectra_tof_and_no_background(self):
        flags = self._create_test_flags(background=False)
        flags['fit_mode'] = 'spectra'
        flags['spectra'] = '143-150'
        runs = "15039-15045"

        fit_results = fit_tof(runs, flags)
        self.assertTrue(isinstance(fit_results, list))
        self.assertEquals(8, len(fit_results))

        spec143 = fit_results[0]
        spec143_data = spec143[0]
        self.assertAlmostEqual(50.0, spec143_data.readX(0)[0])
        self.assertAlmostEqual(562.0, spec143_data.readX(0)[-1])

        self.assertAlmostEqual(2.37897941103748e-06, spec143_data.readY(1)[0])
        self.assertAlmostEqual(3.58226563303213e-05, spec143_data.readY(1)[-1])

        spec150 = fit_results[-1]
        spec150_data = spec150[0]
        self.assertAlmostEqual(50.0, spec150_data.readX(0)[0])
        self.assertAlmostEqual(562.0, spec150_data.readX(0)[-1])

        self.assertAlmostEqual(5.57952304659615e-06, spec150_data.readY(1)[0])
        self.assertAlmostEqual(6.00056973529846e-05, spec150_data.readY(1)[-1])

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

        return flags


if __name__ == '__main__':
    unittest.main()
