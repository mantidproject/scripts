"""These are more integration tests as they will require that the test data is available
and that mantid can be imported
"""
import unittest

from vesuvio.processsing import process_data

class ProcessingTest(unittest.TestCase):

    def test_process_data_with_standard_user_args_and_no_background(self):
        flags = self._create_test_flags(background=False)
        runs = "15039-15045"

        fitted_ws, fitted_params = process_data(runs, flags)

        self.assertEqual(7, fitted_ws.getNumberHistograms())
        self.assertAlmostEqual(50.0, fitted_ws.readX(0)[0])
        self.assertAlmostEqual(562.0, fitted_ws.readX(0)[-1])

        self.assertAlmostEqual(0.000928695463881635, fitted_ws.readY(0)[0])
        self.assertAlmostEqual(0.00722948549525415, fitted_ws.readY(0)[-1])
        self.assertAlmostEqual(1.45746507977816e-05, fitted_ws.readY(1)[0])
        self.assertAlmostEqual(7.33791942084561e-05, fitted_ws.readY(1)[-1])

        self.assertEqual(10, fitted_params.rowCount())

    def test_process_data_with_standard_user_args_and_no_background(self):
        flags = self._create_test_flags(background=True)
        runs = "15039-15045"

        fitted_ws, fitted_params = process_data(runs, flags)

        self.assertEqual(8, fitted_ws.getNumberHistograms())
        self.assertAlmostEqual(50.0, fitted_ws.readX(0)[0])
        self.assertAlmostEqual(562.0, fitted_ws.readX(0)[-1])

        self.assertAlmostEqual(0.000928695463881635, fitted_ws.readY(0)[0])
        self.assertAlmostEqual(0.00722948549525415, fitted_ws.readY(0)[-1])
        self.assertAlmostEqual(-0.00756178413274695, fitted_ws.readY(1)[0])
        self.assertAlmostEqual(0.00355843687365601, fitted_ws.readY(1)[-1])

        self.assertEqual(14, fitted_params.rowCount())

    def _create_test_flags(self, background):
        runs = "15039-15045"
        flags = dict()
        flags['fit_mode'] = 'bank'

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
