#pylint: disable=no-init

from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *

import math
import numpy as np
import scipy.constants as sc
import scipy.stats as stats

#==============================================================================

NEUTRON_MASS_AMU = sc.value('neutron mass in u')

BACKSCATTERING = range(3, 135)
FORWARDSCATTERING = range(135, 199)

MODES = ['SingleDifference', 'DoubleDifference', 'ThickDifference']

BACKSCATTERING = range(3, 135)
FRONTSCATTERING = range(135, 199)

U_FRONTSCATTERING_SAMPLE = ['14025']
U_FRONTSCATTERING_BACKGROUND = ['12570']

U_BACKSCATTERING_SAMPLE = ['12570']
U_BACKSCATTERING_BACKGROUND = ['12571']

# Uranium peak details taken from: 10.1016/j.nima.2010.09.079
U_PEAKS = [
#    Er      dEr    vr          t      dTr
#    meV     meV    m/usec      usec   usec
    [36684., 133.6, 83768.0e-6, 131.  ,0.24],
    [20874., 91.7,  63190.0e-6, 173.8 ,0.38],
    [6672.,  52.4,  35725.0e-6, 307.7 ,1.21]
]

#U mass in amu
U_MASS_IN_AMU = 238.0289

#==============================================================================

class VesuvioGeometryEnergyResolution(PythonAlgorithm):

    _ipf_filename = None
    _evs_ws = None
    _evs = None
    _sample_pos = None
    _l0_dist = None
    _output_table = None

#------------------------------------------------------------------------------

    def summary(self):
        return "Calculates the geometry and energy resolution for VESUVIO."

#------------------------------------------------------------------------------

    def PyInit(self):
        # Sample data
        self.declareProperty(StringArrayProperty("Samples", Direction.Input),
                             doc="Sample run numbers to fit peaks to.")

        self.declareProperty(StringArrayProperty("Background", Direction.Input),
                             doc="Run numbers to use as a background.")

        self.declareProperty(FloatArrayProperty("DSpacings", Direction.Input),
                             doc="List of dSpacings")

        # Optional parameter file
        self.declareProperty(FileProperty("InstrumentParFile", "", action=FileAction.Load,
                                          extensions=["dat", "par"]),
                             doc="A parameter file. If provided the values are used to correct "
                                 "the default instrument values and attach the t0 values to each "
                                 "detector")

        # Output parameter tables
        # self.declareProperty(WorkspaceGroupProperty("ResolutionGroup", "", Direction.Output),
        #                      doc="Resolution values per bank")

        # Output value table
        self.declareProperty(ITableWorkspaceProperty("OutputWorkspace", "", Direction.Output),
                             doc="Mean resolution parameters")

#------------------------------------------------------------------------------

    def PyExec(self):
        self._ipf_filename = self.getPropertyValue("InstrumentParFile")
        self._evs_ws = LoadEmptyVesuvio(InstrumentParFile=self._ipf_filename)

        # Calculate source to sample distance
        self._evs = self._evs_ws.getInstrument()
        source_pos = self._evs.getSource().getPos()
        self._sample_pos = self._evs.getSample().getPos()
        self._l0_dist = source_pos.distance(self._sample_pos)

        # Create mean output table
        self._output_table = CreateEmptyTableWorkspace(OutputWorkspace=self.getPropertyValue("OutputWorkspace"))

        self._output_table.addColumn('str',   'Parameter')
        self._output_table.addColumn('float', 'Mean')
        self._output_table.addColumn('float', 'StdDev')
        self._output_table.addColumn('float', 'StdErr')
        self._output_table.addColumn('float', 'WeightedMean')
        self._output_table.addColumn('float', 'WeightedErr')

        # Calculate resolutions
        data = {}
        data.update(self._process_l_t0())
        data.update(self._process_theta_dw())
        data.update(self._process_de1())

        #TODO

        print data

        self.setProperty("OutputWorkspace", self._output_table)

#------------------------------------------------------------------------------

    def _process_l_t0(self):
        """
        Caluclates delta L and delta t0.

        @return Dictionary of resolution data.
        """

        alg = AlgorithmManager.create("EVSCalibrationFit")
        alg.initialize()
        alg.setRethrows(True)
        alg.setProperty('Samples', U_FRONTSCATTERING_SAMPLE)
        alg.setProperty('Background', U_FRONTSCATTERING_BACKGROUND)
        alg.setProperty('SpectrumRange', [FRONTSCATTERING[0], FRONTSCATTERING[-1]])
        alg.setProperty('Mass', U_MASS_IN_AMU)
        alg.setProperty('Energy', [p[0] for p in U_PEAKS])
        alg.setProperty('InstrumentParameterFile', self._ipf_filename)
        alg.setProperty('OutputWorkspace', '__l0_fit')
        alg.setProperty('CreateOutput', False)
        alg.execute()

        alg = AlgorithmManager.create("EVSCalibrationFit")
        alg.initialize()
        alg.setRethrows(True)
        alg.setProperty('Samples', U_BACKSCATTERING_SAMPLE)
        alg.setProperty('Background', U_BACKSCATTERING_BACKGROUND)
        alg.setProperty('SpectrumRange', [BACKSCATTERING[0], BACKSCATTERING[-1]])
        alg.setProperty('Mass', U_MASS_IN_AMU)
        alg.setProperty('Energy', [p[0] for p in U_PEAKS])
        alg.setProperty('InstrumentParameterFile', self._ipf_filename)
        alg.setProperty('OutputWorkspace', '__t0_fit')
        alg.setProperty('CreateOutput', False)
        alg.execute()

        t0_peak_params = mtd['__t0_fit_Peak_Parameters']
        l0_peak_params = mtd['__l0_fit_Peak_Parameters']

        # For forward scattering only L1 is used as gamma rays are detected
        forward_l = [self._l0_dist] * len(FRONTSCATTERING)
        # For backscattering total flight path is used as neutrons are detected
        back_l0 = [self._l0_dist + self._get_l1(det_no) for det_no in BACKSCATTERING]

        # Calculate resolution
        back_params = self._calculate_l_t0(BACKSCATTERING, t0_peak_params, back_l0)
        forward_params = self._calculate_l_t0(FRONTSCATTERING, l0_peak_params, forward_l)

        # Get data from param tables
        max_float = np.finfo(float).max
        back_l_data = np.clip(back_params.column('A1'), 0, max_float)
        back_t0_data = np.clip(back_params.column('A0'), 0, max_float)
        forward_l_data = np.clip(forward_params.column('A1'), 0, max_float)
        forward_t0_data = np.clip(forward_params.column('A0'), 0, max_float)

        self._add_param_table_row('Backscattering L', back_l_data)
        self._add_param_table_row('Forward scattering L', forward_l_data)
        self._add_param_table_row('Backscattering t0', back_t0_data)
        self._add_param_table_row('Forward scattering t0', forward_t0_data)

        data = {
            'dt0': np.ravel([back_t0_data, forward_t0_data]),
            'dL0': np.ravel([back_l_data, forward_l_data])
        }

        return data

#------------------------------------------------------------------------------

    def _calculate_l_t0(self, detectors, calibration_params, distance):
        """
        Calculates resolution for given detectors.

        @param detectors List of detector numbers
        @param calibration_params Parameters from calibration fits
        @param distance List of distances for detectors (L for backscattering, L1 for forward scattering)
        @return Workspace containing resolution parameters per detector
        """
        # Create a workspace of the three peaks against 1/v1^2
        wks = WorkspaceFactory.Instance().create('Workspace2D', len(detectors), 3, 3)

        for detector in detectors:
            det_index = detector - detectors[0]

            x_data = []
            for peak in range(3):
                peak_position = calibration_params.getItem(peak).column('f1.PeakCentre')[det_index]
                x_data.append(1.0/(distance[det_index]/peak_position)**2)

                params = calibration_params.getItem(peak)
                sigma = params.column('f1.Sigma')[det_index]
                sigma_err = params.column('f1.Sigma_Err')[det_index]

                u_peak = U_PEAKS[peak]
                wks.dataY(det_index)[peak] = (sigma ** 2) - (u_peak[4]**2)
                wks.dataE(det_index)[peak] = 2*sigma*sigma_err

            wks.setX(det_index, np.array(x_data))

        AnalysisDataService.Instance().addOrReplace('__bank_data', wks)

        # Perform a linear fit of each spectra
        fit = AlgorithmManager.Instance().create('PlotPeakByLogValue')
        fit.initialize()
        fit.setChild(True)
        fit.setProperty('Function', 'name=LinearBackground')
        fit.setProperty('Input', '__bank_data,v0:132')
        fit.setProperty('OutputWorkspace', 'backscattering_params')
        fit.execute()

        DeleteWorkspace('__bank_data')
        params = fit.getProperty('OutputWorkspace').value

        # Process fit parameters
        for index, detector in enumerate(detectors):
            params.setCell(index, 0, detector)

            t0_val = params.cell(index, 1)
            l_dist = params.cell(index, 3)

            # Set negative values to zero, otherwise take square root
            if t0_val > 0:
                t0_val = np.sqrt(t0_val)
            else:
                t0_val = 0

            if l_dist > 0:
                l_dist = np.sqrt(l_dist)
            else:
                l_dist = 0

            params.setCell(index, 1, t0_val)
            params.setCell(index, 3, l_dist)

        return params

#------------------------------------------------------------------------------

    def _process_theta_dw(self):
        """
        Caluclates delta theta and effecitve detector width.

        @return Dictionary of resolution data.
        """

        # Run fitting
        alg = AlgorithmManager.create("EVSCalibrationFit")
        alg.initialize()
        alg.setRethrows(True)
        alg.setProperty('Samples', self.getPropertyValue('Samples'))
        alg.setProperty('Background', self.getPropertyValue('Background'))
        alg.setProperty('Function', 'Gaussian')
        alg.setProperty('Mode', 'FoilOut')
        alg.setProperty('InstrumentParameterFile', self.getPropertyValue('InstrumentParFile'))
        alg.setProperty('DSpacings', self.getProperty('DSpacings').value)
        alg.setProperty('OutputWorkspace', 'delta_theta')
        alg.setProperty('CreateOutput', False)
        alg.execute()
        params = mtd['delta_theta_Peak_Parameters']

        back_theta_data = []
        forward_theta_data = []
        back_dw_data = []
        forward_dw_data = []

        for spec_no in BACKSCATTERING:
            delta_theta, effective_width = self._calc_delta_theta(params, spec_no)
            back_theta_data.append(delta_theta)
            back_dw_data.append(effective_width)

        back_theta_data = np.array(back_theta_data)
        back_dw_data = np.array(back_dw_data)

        for spec_no in FORWARDSCATTERING:
            delta_theta, effective_width = self._calc_delta_theta(params, spec_no)
            forward_theta_data.append(delta_theta)
            forward_dw_data.append(effective_width)

        forward_theta_data = np.array(forward_theta_data)
        forward_dw_data = np.array(forward_dw_data)

        self._add_param_table_row('Backscattering dTheta', back_theta_data)
        self._add_param_table_row('Forward scattering dTheta', forward_theta_data)
        self._add_param_table_row('Backscattering dW', back_dw_data)
        self._add_param_table_row('Forward scattering dW', forward_dw_data)

        data = {
            'dTheta': np.ravel([back_theta_data, forward_theta_data]),
            'EffectiveDetWidth': np.ravel([back_dw_data, forward_dw_data])
        }

        return data

#------------------------------------------------------------------------------

    def _calc_delta_theta(self, params, spec_no):
        """
        Calculates delta Theta for a given spectrum number using Gaussian fits
        to (4) Bragg peaks for a Lead sample.

        @param params Fit parameter table
        @param spec_no Spectrum no to calculate for
        @return Average delta Theta over all fitted peaks
        """
        spec_col_idx = params.column('Spectrum').index(spec_no)

        sigma_col_names = [c for c in params.getColumnNames() if 'Sigma' in c and 'Err' not in c]
        peak_centre_col_names = [c for c in params.getColumnNames() if 'PeakCentre' in c and 'Err' not in c]

        #TODO
        delta_L0 = 0.0
        delta_L1 = 0.0
        delta_t0 = 0.0

        l1_dist = self._get_l1(spec_no)

        delta_thetas = []
        delta_w = []
        for sigma_name, peak_centre_name in zip(sigma_col_names, peak_centre_col_names):
            sigma = params.cell(sigma_name, spec_col_idx)
            peak_centre = params.cell(peak_centre_name, spec_col_idx)
            theta = self._get_theta(spec_no)

            if spec_no in BACKSCATTERING:
                delta_theta = 2 * math.tan(theta / 2) * \
                              math.sqrt((sigma / peak_centre)**2 \
                                        - ((delta_L0**2 + delta_L1**2) / (self._l0_dist + l1_dist)**2) \
                                        - (delta_t0**2 / peak_centre**2))
            elif spec_no in FORWARDSCATTERING:
                delta_theta = math.tan(theta) * (sigma / peak_centre)
            else:
                raise RuntimeError()

            delta_thetas.append(np.degrees(abs(delta_theta)))
            delta_w.append(abs(delta_theta) * l1_dist * 100)

        delta_thetas = np.array(delta_thetas)
        delta_thetas_avg = np.mean(delta_thetas)
        logger.debug('Spectrum {0}: dThetas={1}, mean={2}'.format(spec_no, delta_thetas, delta_thetas_avg))

        delta_w = np.array(delta_w)
        delta_w_avg = np.mean(delta_w)

        return (delta_thetas_avg, delta_w_avg)

#------------------------------------------------------------------------------

    def _process_de1(self):
        """
        Caluclates delta E1 for Lorentzian and Gaussian components.

        @return Dictionary of resolution data.
        """
        data = {}

        for mode in MODES:
            self._energy_res_fit(mode)

            de1_gauss, de1_lorentz = self._calculate_de1(mtd['{0}_Peaks_Peak_0_Parameters'.format(mode)])

            short_mode = ''.join([c for c in mode if c.isupper()])
            data['dE1_%s_Gauss' % short_mode] = de1_gauss
            data['dE1_%s_Lorentz' % short_mode] = de1_lorentz

            de1_gauss = de1_gauss[np.isfinite(de1_gauss)]
            de1_lorentz = de1_lorentz[np.isfinite(de1_lorentz)]

            self._add_param_table_row('%s_Back_dE1_Gauss' & short_mode,
                                      de1_gauss[:FORWARDSCATTERING[0]])
            self._add_param_table_row('%s_Back_dE1_Lorentz' % short_mode,
                                      de1_lorentz[:FORWARDSCATTERING[0]])

            if mode == 'SingleDifference':
                self._add_param_table_row('%s_Forward_dE1_Gauss' % short_mode,
                                          de1_gauss[FORWARDSCATTERING[0]:])
                self._add_param_table_row('%s_Forward_dE1_Lorentz' % short_mode,
                                          de1_lorentz[FORWARDSCATTERING[0]:])

        return data

#------------------------------------------------------------------------------

    def _energy_res_fit(self, mode):
        """
        Performs calibration fitting.

        @param mode Mode in which to fit data
        """
        ws_name = '{0}_Peaks'.format(mode)

        alg = AlgorithmManager.create('EVSCalibrationFit')
        alg.initialize()
        alg.setRethrows(True)
        alg.setProperty('Samples', self.getPropertyValue('Samples'))
        alg.setProperty('Background', self.getPropertyValue('Background'))
        alg.setProperty('Mode', mode)
        alg.setProperty('Function', 'Voigt')
        alg.setProperty('InstrumentParameterFile', self.getPropertyValue('InstrumentParFile'))
        alg.setProperty('CreateOutput', False)
        alg.setProperty('OutputWorkspace', ws_name)
        alg.execute()

#------------------------------------------------------------------------------

    def _calculate_de1(self, parameters):
        """
        Calculates delta-E1 for Lorentzian and Gaussian components.

        @param parameters Fit parameters
        @return Tuple of (dE1_Gauss, dE1_Lorentz)
        """
        spectrum_nos = parameters.column('Spectrum')
        positions = parameters.column('f1.LorentzPos')
        widths_lorentz = parameters.column('f1.LorentzFWHM')
        widths_gauss = parameters.column('f1.GaussianFWHM')

        #TODO
        delta_theta = 0.0
        delta_L0 = 0.0
        delta_L1 = 0.0
        delta_t0 = 0.0

        de1_gauss = []
        de1_lorentz = []

        for spec, pos, fwhm_lorentz, fwhm_gauss in zip(spectrum_nos, positions, widths_lorentz, widths_gauss):
            l1_dist = self._get_l1(spec-1)

            # Lorentzian component
            if pos == 0.0 or fwhm_lorentz == 0.0:
                delta_e1_lorentz = np.nan
            else:
                delta_e1_lorentz = self._convert_to_energy(l1_dist, pos, fwhm_lorentz)

            # Gaussian component
            if pos == 0.0 or fwhm_gauss == 0.0:
                delta_e1_gauss = np.nan
            else:
                delta_e1_gauss = self._convert_to_energy(l1_dist, pos, fwhm_gauss)
                vel = delta_e1_gauss**2 - delta_theta**2 - delta_L0**2 - delta_L1**2 - delta_t0**2
                if vel < 0:
                    delta_e1_gauss = np.nan
                else:
                    delta_e1_gauss = math.sqrt(vel) / 1.775

            de1_gauss.append(delta_e1_gauss)
            de1_lorentz.append(delta_e1_lorentz)

        return (np.array(de1_gauss), np.array(de1_lorentz))

#------------------------------------------------------------------------------

    def _convert_to_energy(self, l1_dist, position, fwhm):
        """
        Converts a peak width to energy.

        @param l1_dist Distance from sample to detector (m)
        @param position Peak position in time of flight (us)
        @param fwhm FWHM in time of flight (us)
        """
        # Calculate energy at the peak position
        d_time = position * sc.micro # us to s
        neutron_v1 = (self._l0_dist + l1_dist) / d_time # ms-1
        peak_e = 5.2276e-6 * neutron_v1**2

        # Calculate HWHM of peak in meV
        energy = peak_e * fwhm / position

        return energy

#------------------------------------------------------------------------------

    def _get_l1(self, spec_no):
        """
        Gets the L1 distance for a given spectrum number.

        @param det_no The detector number
        @return The L1 distance
        """
        ws_idx = self._evs_ws.getIndexFromSpectrumNumber(spec_no)
        det_pos = self._evs_ws.getDetector(ws_idx).getPos()
        dist = self._sample_pos.distance(det_pos)
        return dist

#------------------------------------------------------------------------------

    def _get_theta(self, spec_no):
        """
        Gets Theta (in degrees) for a given spectrum number.

        @param spec_no Spectrum number containing detector
        @return Theta in degrees for detector
        """
        ws_idx = self._evs_ws.getIndexFromSpectrumNumber(spec_no)
        detector = self._evs_ws.getDetector(ws_idx)
        theta = self._evs_ws.detectorTwoTheta(detector)
        return math.degrees(theta)

#------------------------------------------------------------------------------

    def _add_param_table_row(self, name, data):
        """
        Adds a parameter row to the output table.

        @param name Name of the parameter
        @param data Data values
        """
        mean = np.mean(data)
        weights = 1.0 / np.abs(np.repeat(mean, data.size) - data)
        weighted_data = data * weights / np.sum(weights)

        self._output_table.addRow([name,
                                   mean,
                                   np.std(data),
                                   stats.sem(data),
                                   np.average(data, weights=weights),
                                   stats.sem(weighted_data)])

#==============================================================================

AlgorithmFactory.subscribe(VesuvioGeometryEnergyResolution)
