#pylint: disable=no-init

from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *

import math
import numpy as np
import scipy.stats as stats

BACKSCATTERING = range(3, 135)
FORWARDSCATTERING = range(135, 199)

#----------------------------------------------------------------------------------------

class VesuvioThetaResolution(PythonAlgorithm):

    _evs_ws = None
    _evs = None
    _sample_pos = None
    _l0_dist = None
    _fit_params = None
    _output_table = None

    def summary(self):
        return "Calculates the resolution of the scattering angle for VESUVIO."

    def PyInit(self):
        # Sample data
        self.declareProperty(StringArrayProperty("Samples", Direction.Input),
                             doc="Sample run numbers to fit peaks to.")

        self.declareProperty(StringArrayProperty("Background", Direction.Input),
                             doc="Run numbers to use as a background.")

        self.declareProperty(FloatArrayProperty("DSpacings", Direction.Input),
                             doc="List of dSpacings")

        # Optional parameter file
        self.declareProperty(FileProperty("InstrumentParFile", "", action=FileAction.OptionalLoad,
                                          extensions=["dat", "par"]),
                             doc="An optional IP file. If provided the values are used to correct "
                                 "the default instrument values and attach the t0 values to each "
                                 "detector")

        # Resolution params
        self.declareProperty("DeltaL0", 0.0,
                             doc="Resolution in the incident flight path length")

        self.declareProperty("DeltaL1", 0.0,
                             doc="Resolution in the final flight path length")

        self.declareProperty("DeltaT0", 0.0,
                             doc="Resolution in the pulse offset time")

        # Output parameter tables
        self.declareProperty(ITableWorkspaceProperty("ForwardParameters", "", Direction.Output),
                             doc="Resolution values for forward scattering")
        self.declareProperty(ITableWorkspaceProperty("BackParameters", "", Direction.Output),
                             doc="Resolution parameters for backscattering")

        # Output value table
        self.declareProperty(ITableWorkspaceProperty("OutputWorkspace", "", Direction.Output),
                             doc="Mean resolution parameters")

#----------------------------------------------------------------------------------------

    def PyExec(self):
        self._evs_ws = LoadEmptyVesuvio(InstrumentParFile=self.getPropertyValue("InstrumentParFile"))

        # Calculate source to sample distance
        self._evs = self._evs_ws.getInstrument()
        source_pos = self._evs.getSource().getPos()
        self._sample_pos = self._evs.getSample().getPos()
        self._l0_dist = source_pos.distance(self._sample_pos)

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
        self._fit_params = mtd['delta_theta_Peak_Parameters']

        back_params = self._create_output_table(self.getPropertyValue("BackParameters"))
        forward_params = self._create_output_table(self.getPropertyValue("ForwardParameters"))

        for spec_no in BACKSCATTERING:
            row = [spec_no]
            row.extend(self._calc_delta_theta(spec_no))
            back_params.addRow(row)
        for spec_no in FORWARDSCATTERING:
            row = [spec_no]
            row.extend(self._calc_delta_theta(spec_no))
            forward_params.addRow(row)

        DeleteWorkspace(self._evs_ws)

        self.setProperty("BackParameters", back_params)
        self.setProperty("ForwardParameters", forward_params)

        # Get data from param tables
        back_theta_data = np.array(back_params.column('dTheta'))
        forward_theta_data = np.array(forward_params.column('dTheta'))
        back_dw_data = np.array(back_params.column('dW'))
        forward_dw_data = np.array(forward_params.column('dW'))

        self._output_table = CreateEmptyTableWorkspace(OutputWorkspace=self.getPropertyValue("OutputWorkspace"))

        self._output_table.addColumn('str', 'ResolutionParameter')
        self._output_table.addColumn('float', 'Mean')
        self._output_table.addColumn('float', 'StdDev')
        self._output_table.addColumn('float', 'StdErr')
        self._output_table.addColumn('float', 'WeightedMean')
        self._output_table.addColumn('float', 'WeightedErr')

        self._add_param_table_row('Backscattering dTheta', back_theta_data)
        self._add_param_table_row('Forward scattering dTheta', forward_theta_data)
        self._add_param_table_row('Backscattering dW', back_dw_data)
        self._add_param_table_row('Forward scattering dW', forward_dw_data)

        self.setProperty("OutputWorkspace", self._output_table)

#----------------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------------

    def _calc_delta_theta(self, spec_no):
        """
        Calculates delta Theta for a given spectrum number using Gaussian fits
        to (4) Bragg peaks for a Lead sample.

        @param spec_no Spectrum no to calculate for
        @return Average delta Theta over all fitted peaks
        """
        spec_col_idx = self._fit_params.column('Spectrum').index(spec_no)

        sigma_col_names = [c for c in self._fit_params.getColumnNames() if 'Sigma' in c and 'Err' not in c]
        peak_centre_col_names = [c for c in self._fit_params.getColumnNames() if 'PeakCentre' in c and 'Err' not in c]

        delta_L0 = self.getProperty('DeltaL0').value
        delta_L1 = self.getProperty('DeltaL1').value
        delta_t0 = self.getProperty('DeltaT0').value
        l1_dist = self._get_l1(spec_no)

        delta_thetas = []
        delta_w = []
        for sigma_name, peak_centre_name in zip(sigma_col_names, peak_centre_col_names):
            sigma = self._fit_params.cell(sigma_name, spec_col_idx)
            peak_centre = self._fit_params.cell(peak_centre_name, spec_col_idx)
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

        return [delta_thetas_avg, delta_w_avg]

#----------------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------------

    def _create_output_table(self, name):
        """
        Creates an output table workspace.

        @param name Nmae of the output workspace
        @return Created table workspace
        """
        table = CreateEmptyTableWorkspace(OutputWorkspace=name)
        table.addColumn('int', 'Spectrum')
        table.addColumn('float', 'dTheta')
        table.addColumn('float', 'dW')
        return table

#----------------------------------------------------------------------------------------

AlgorithmFactory.subscribe(VesuvioThetaResolution)
