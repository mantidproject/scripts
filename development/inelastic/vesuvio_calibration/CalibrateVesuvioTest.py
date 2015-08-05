#pylint: disable=too-many-instance-attributes,too-many-public-methods,invalid-name,too-many-arguments

import unittest
import numpy as np
import scipy.constants
import scipy.stats

from mantid.api import FileFinder, WorkspaceGroup
from mantid.simpleapi import *
from CalibrateVesuvio import (calculate_r_theta, FRONTSCATTERING_RANGE, DETECTOR_RANGE, \
                              BACKSCATTERING_RANGE, ENERGY_ESTIMATE, MEV_CONVERSION, \
                              U_FRONTSCATTERING_SAMPLE, U_FRONTSCATTERING_BACKGROUND, \
                              U_BACKSCATTERING_SAMPLE, U_BACKSCATTERING_BACKGROUND,\
                              U_MASS, U_PEAK_ENERGIES)

class EVSCalibrationTest(unittest.TestCase):

    _run_range = None
    _spec_list = None
    _mode = None
    _function = None
    _background = None
    _energy_estimates = None
    _mass = None
    _d_spacings = None

    def load_ip_file(self):
        param_names = ['spectrum', 'theta', 't0', 'L0', 'L1']
        file_data = np.loadtxt(self._parameter_file, skiprows=3, usecols=[0,2,3,4,5], unpack=True)

        params = {}
        for name, column in zip(param_names, file_data):
            params[name] = column

        return params

    def _setup_copper_test(self):
        self._run_range = "17087-17088"
        self._background = "17086"
        #Mass of copper in amu
        self._mass = 63.546
        #d-spacings of a copper sample
        self._d_spacings = np.array([2.0865, 1.807, 1.278, 1.0897])
        self._d_spacings.sort()
        self._energy_estimates = np.array(ENERGY_ESTIMATE)
        self._mode = 'FoilOut'

    def _setup_lead_test(self):
        self._run_range = "17083-17084"
        self._background = "17086"
        #Mass of a lead in amu
        self._mass = 207.19
        #d-spacings of a lead sample
        self._d_spacings = np.array([1.489, 1.750, 2.475, 2.858])
        self._d_spacings.sort()
        self._energy_estimates = np.array(ENERGY_ESTIMATE)
        self._mode = 'FoilOut'

    def _setup_niobium_test(self):
        self._run_range = "17089-17091"
        self._background = "17086"
        #Mass of a lead in amu
        self._mass = 92.906
        #d-spacings of a lead sample
        self._d_spacings = np.array([2.3356, 1.6515, 1.3484, 1.1678])
        self._d_spacings.sort()
        self._energy_estimates = np.array(ENERGY_ESTIMATE)
        self._mode = 'FoilOut'

    def _setup_E1_fit_test(self):
        self._spec_list = DETECTOR_RANGE
        self._d_spacings = []
        self._background = ''
        self._mode = 'SingleDifference'

    def _setup_uranium_test(self):
        self._function= 'Gaussian'
        self._mass = U_MASS
        self._d_spacings = []
        self._energy_estimates = np.array(U_PEAK_ENERGIES)
        self._energy_estimates.sort()
        self._mode = 'FoilOut'

    def tearDown(self):
        mtd.clear()


class EVSCalibrationAnalysisTests(EVSCalibrationTest):

    _output_workspace = None

    def setUp(self):
        self._calc_L0 = False
        self._parameter_file = FileFinder.getFullPath("IP0005.par")
        self._calibrated_params = self.load_ip_file()
        self._iterations = 2

    def test_copper(self):
        self._setup_copper_test()
        self._output_workspace = "copper_analysis_test"

        params_table = self.run_evs_calibration_analysis()
        self.assert_theta_parameters_match_expected(params_table)

    def test_lead(self):
        self._setup_lead_test()
        self._output_workspace = "lead_analysis_test"

        params_table = self.run_evs_calibration_analysis()
        self.assert_theta_parameters_match_expected(params_table)

    # def test_niobium(self):
    #     self._setup_niobium_test()
    #     self._output_workspace = "niobium_analysis_test"
    #     self._iterations = 1

    #     params_table = self.run_evs_calibration_analysis()
    #     self.assert_theta_parameters_match_expected(params_table)

    def test_copper_with_uranium(self):
        self._setup_copper_test()
        self._calc_L0 = True
        self._output_workspace = "copper_analysis_test"

        params_table = self.run_evs_calibration_analysis()
        self.assert_theta_parameters_match_expected(params_table)

    def test_lead_with_uranium(self):
        self._setup_lead_test()
        self._calc_L0 = True
        self._output_workspace = "lead_analysis_test"

        params_table = self.run_evs_calibration_analysis()
        self.assert_theta_parameters_match_expected(params_table)

    def tearDown(self):
        mtd.clear()

    #------------------------------------------------------------------
    # Misc helper functions
    #------------------------------------------------------------------

    def assert_theta_parameters_match_expected(self, params_table, rel_tolerance=0.4):
        thetas = params_table.column('theta')
        actual_thetas = self._calibrated_params['theta']

        self.assertFalse(np.isnan(thetas).any())
        self.assertFalse(np.isinf(thetas).any())
        np.testing.assert_allclose(actual_thetas, thetas, rtol=rel_tolerance)

    def run_evs_calibration_analysis(self):
        EVSCalibrationAnalysis(self._output_workspace, Samples=self._run_range, Background=self._background,
                               InstrumentParameterFile=self._parameter_file, Mass=self._mass, DSpacings=self._d_spacings,
                               Iterations=self._iterations, CalculateL0=self._calc_L0)

        last_fit_iteration = self._output_workspace + '_Iteration_%d' % (self._iterations-1)
        return mtd[last_fit_iteration]


class EVSCalibrationFitTests(EVSCalibrationTest):

    _output_workspace = None

    def setUp(self):
        self._function = 'Voigt'
        self._parameter_file = FileFinder.getFullPath("IP0005.par")
        self._calibrated_params = self.load_ip_file()
        self._mode = 'FoilOut'
        self._energy_estimates = np.array([ENERGY_ESTIMATE])

    def test_fit_bragg_peaks_copper(self):
        self._setup_copper_test()
        self._spec_list = DETECTOR_RANGE
        self._output_workspace = "copper_bragg_peak_fit"

        params_table = self.run_evs_calibration_fit()
        expected_values = self.calculate_theta_peak_positions()
        self.assert_fitted_positions_match_expected(expected_values, params_table, rel_tolerance=0.7)

    def test_fit_bragg_peaks_lead(self):
        self._setup_lead_test()
        self._spec_list = DETECTOR_RANGE
        self._output_workspace = "lead_bragg_peak_fit"

        expected_values = self.calculate_theta_peak_positions()
        params_table = self.run_evs_calibration_fit()
        self.assert_fitted_positions_match_expected(expected_values, params_table, rel_tolerance=0.7)

    def test_fit_peaks_copper(self):
        self._setup_copper_test()
        self._setup_E1_fit_test()
        self._output_workspace = "copper_peak_fit"

        expected_values = self.calculate_energy_peak_positions()
        params_table = self.run_evs_calibration_fit()
        self.assert_fitted_positions_match_expected(expected_values, params_table, rel_tolerance=0.7)

    def test_fit_peaks_lead(self):
        self._setup_copper_test()
        self._setup_E1_fit_test()
        self._output_workspace = "lead_peak_fit"

        expected_values = self.calculate_energy_peak_positions()
        params_table = self.run_evs_calibration_fit()
        self.assert_fitted_positions_match_expected(expected_values, params_table, rel_tolerance=0.7)

    def test_fit_frontscattering_uranium(self):
        self._setup_uranium_test()
        self._run_range = U_FRONTSCATTERING_SAMPLE
        self._background = U_FRONTSCATTERING_BACKGROUND
        self._spec_list = FRONTSCATTERING_RANGE
        self._output_workspace = 'uranium_peak_fit_front'

        expected_values = self.calculate_energy_peak_positions()
        params_table = self.run_evs_calibration_fit()

        self.assert_fitted_positions_match_expected(expected_values, params_table, rel_tolerance=0.1)

    def test_fit_backscattering_uranium(self):
        self._setup_uranium_test()
        self._run_range = U_BACKSCATTERING_SAMPLE
        self._background = U_BACKSCATTERING_BACKGROUND
        self._spec_list = BACKSCATTERING_RANGE
        self._output_workspace = 'uranium_peak_fit_back'

        expected_values = self.calculate_energy_peak_positions()
        params_table = self.run_evs_calibration_fit()
        self.assert_fitted_positions_match_expected(expected_values, params_table, rel_tolerance=0.01, ignore_zero=True)

    #------------------------------------------------------------------
    # Misc Helpers functions
    #------------------------------------------------------------------

    def assert_fitted_positions_match_expected(self, expected_positions, params_table,
                                               rel_tolerance=1e-7, abs_tolerance=0,
                                               ignore_zero=False):
        """
        Check that each of the fitted peak positions match the expected
        positions in time of flight calculated from the parameter file
        within a small tolerance.
        """
        if isinstance(params_table, WorkspaceGroup):
            params_table = params_table.getNames()[0]
            params_table = mtd[params_table]

        column_names = self.find_all_peak_positions(params_table)

        for name, expected_position in zip(column_names, expected_positions):
            position = np.array(params_table.column(name))

            if ignore_zero:
                expected_position, position = self.mask_bad_detector_readings(expected_position, position)

            self.assertFalse(np.isnan(position).any())
            self.assertFalse(np.isinf(position).any())
            np.set_printoptions(threshold=np.nan)
            np.testing.assert_allclose(expected_position, position, rtol=rel_tolerance, atol=abs_tolerance)

    def mask_bad_detector_readings(self, expected_positions, actual_positions):
        """
        Mask values that are very close to zero.

        This handles the case where some of the uranium runs have a missing entry
        for one of the detectors.
        """
        zero_mask = np.where(actual_positions > 1e-10)
        expected_positions = expected_positions[zero_mask]
        actual_positions = actual_positions[zero_mask]
        return expected_positions, actual_positions

    def run_evs_calibration_fit(self):
        EVSCalibrationFit(Samples=self._run_range, Background=self._background,
                          Function=self._function, Mode=self._mode, SpectrumRange=self._spec_list,
                          InstrumentParameterFile=self._parameter_file, DSpacings=self._d_spacings,
                          Energy=self._energy_estimates,
                          OutputWorkspace=self._output_workspace, CreateOutput=False)

        return mtd[self._output_workspace + '_Peak_Parameters']

    def find_all_peak_positions(self, params_table):
        filter_errors_func = lambda name: ('LorentzPos' in name or 'PeakCentre' in name) and 'Err' not in name
        column_names = params_table.getColumnNames()
        column_names = filter(filter_errors_func, column_names)
        return column_names

    def calculate_energy_peak_positions(self):
        """
        Using the calibrated values to calculate the expected
        position of the L0/L1/E1 peak in time of flight.
        """
        lower, upper = self._spec_list[0] - DETECTOR_RANGE[0], self._spec_list[1] - DETECTOR_RANGE[0] +1

        L0 = self._calibrated_params['L0'][lower:upper]
        L1 = self._calibrated_params['L1'][lower:upper]
        t0 = self._calibrated_params['t0'][lower:upper]
        thetas = self._calibrated_params['theta'][lower:upper]
        r_theta = calculate_r_theta(self._mass, thetas)

        t0 /= 1e+6

        energy_estimates = np.copy(self._energy_estimates)
        energy_estimates = energy_estimates.reshape(1, energy_estimates.size).T
        energy_estimates = energy_estimates * MEV_CONVERSION

        v1 = np.sqrt(2 * energy_estimates / scipy.constants.m_n)
        tof = ((L0 * r_theta + L1) / v1) + t0
        tof *= 1e+6

        return tof

    def calculate_theta_peak_positions(self):
        """
        Using the calibrated values of theta calculate the expected
        peak position in time of flight.
        """
        L0 = self._calibrated_params['L0']
        L1 = self._calibrated_params['L1']
        t0 = self._calibrated_params['t0']
        thetas = self._calibrated_params['theta']

        t0 /= 1e+6

        d_spacings = np.copy(self._d_spacings)
        d_spacings *= 1e-10
        d_spacings = d_spacings.reshape(1, d_spacings.size).T

        lambdas = 2 * d_spacings * np.sin(np.radians(thetas) / 2)
        tof = (lambdas * scipy.constants.m_n * (L0 + L1)) / scipy.constants.h + t0
        tof *= 1e+6

        return tof

if __name__ == '__main__':
    unittest.main()
