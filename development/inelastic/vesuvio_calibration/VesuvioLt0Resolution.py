#pylint: disable=no-init

from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *

import numpy as np
import scipy.stats as stats

BACKSCATTERING = range(3, 135)
FORWARDSCATTERING = range(135, 199)

# Uranium peak details taken from: 10.1016/j.nima.2010.09.079
U_PEAKS = [
#    Er      dEr    vr          t      dTr
#    meV     meV    m/usec      usec   usec
    [36684., 133.6, 83768.0e-6, 131.  ,0.24],
    [20874., 91.7,  63190.0e-6, 173.8 ,0.38],
    [6672.,  52.4,  35725.0e-6, 307.7 ,1.21]
]

#----------------------------------------------------------------------------------------

class VesuvioLt0Resolution(PythonAlgorithm):

    _output_table = None

    def summary(self):
        return "Calculates the resolution of the total flight path and t0 for VESUVIO."

    def PyInit(self):
        # Input parameter tables
        self.declareProperty(WorkspaceGroupProperty("ParametersT0", "", Direction.Input),
                             doc="Peak parameters fitted for t0 in EVSCalibrationAnalysis")
        self.declareProperty(WorkspaceGroupProperty("ParametersL0", "", Direction.Input),
                             doc="Peak parameters fitted for L0 in EVSCalibrationAnalysis")

        # Optional parameter file
        self.declareProperty(FileProperty("InstrumentParFile", "", action=FileAction.OptionalLoad,
                                          extensions=["dat", "par"]),
                             doc="An optional IP file. If provided the values are used to correct "
                                 "the default instrument values and attach the t0 values to each "
                                 "detector")

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
        t0_peak_params = mtd[self.getPropertyValue("ParametersT0")]
        l0_peak_params = mtd[self.getPropertyValue("ParametersL0")]

        evs_ws = LoadEmptyVesuvio(InstrumentParFile=self.getPropertyValue("InstrumentParFile"))

        # Calculate source to sample distance
        evs = evs_ws.getInstrument()
        source_pos = evs.getSource().getPos()
        sample_pos = evs.getSample().getPos()
        l0_dist = source_pos.distance(sample_pos)

        # For forward scattering only L1 is used as gamma rays are detected
        forward_l = [l0_dist] * len(FORWARDSCATTERING)
        # For backscattering total flight path is used as neutrons are detected
        back_l0 = [l0_dist + self._get_l1(evs_ws, det_no) for det_no in BACKSCATTERING]

        DeleteWorkspace(evs_ws)

        # Calculate resolution
        back_fit_params = self._calculate(BACKSCATTERING, t0_peak_params, back_l0)
        forward_fit_params = self._calculate(FORWARDSCATTERING, l0_peak_params, forward_l)

        # Rename the headers and discard the cost function values
        back_params = CreateEmptyTableWorkspace(OutputWorkspace=self.getPropertyValue("BackParameters"))
        forward_params = CreateEmptyTableWorkspace(OutputWorkspace=self.getPropertyValue("ForwardParameters"))
        self._add_output_columns(back_params, back_fit_params)
        self._add_output_columns(forward_params, forward_fit_params)

        self.setProperty("BackParameters", back_params)
        self.setProperty("ForwardParameters", forward_params)

        # Get data from param tables
        max_float = np.finfo(float).max
        back_l_data = np.clip(back_params.column('L'), 0, max_float)
        back_t0_data = np.clip(back_params.column('t0'), 0, max_float)
        forward_l_data = np.clip(forward_params.column('L'), 0, max_float)
        forward_t0_data = np.clip(forward_params.column('t0'), 0, max_float)

        self._output_table = CreateEmptyTableWorkspace(OutputWorkspace=self.getPropertyValue("OutputWorkspace"))

        self._output_table.addColumn('str', 'ResolutionParameter')
        self._output_table.addColumn('float', 'Mean')
        self._output_table.addColumn('float', 'StdDev')
        self._output_table.addColumn('float', 'StdErr')
        self._output_table.addColumn('float', 'WeightedMean')
        self._output_table.addColumn('float', 'WeightedErr')

        self._add_param_table_row('Backscattering L', back_l_data)
        self._add_param_table_row('Backscattering t0', back_t0_data)
        self._add_param_table_row('Forward scattering L', forward_l_data)
        self._add_param_table_row('Forward scattering t0', forward_t0_data)

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

    def _get_l1(self, evs_ws, det_no):
        """
        Gets the L1 distance for a given detector.

        @param evs_ws Workspace containg the VESUVIO instrument
        @param det_no The detector number
        @return The L1 distance
        """
        evs = evs_ws.getInstrument()
        sample_pos = evs.getSample().getPos()
        det_pos = evs_ws.getDetector(det_no).getPos()
        dist = sample_pos.distance(det_pos)
        return dist

#----------------------------------------------------------------------------------------

    def _calculate(self, detectors, calibration_params, distance):
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

#----------------------------------------------------------------------------------------

    def _add_output_columns(self, table, fit_table):
        """
        Move the columns from the fit parameter workspace to the output workspace.

        @param table Output table workspace
        @param fit_table Fit parameter table workspace
        """
        table.addColumn('float', 'Detector')
        table.addColumn('float', 't0')
        table.addColumn('float', 't0_Err')
        table.addColumn('float', 'L')
        table.addColumn('float', 'L_Err')

        for row_idx in range(fit_table.rowCount()):
            old_row = fit_table.row(row_idx)
            table.addRow([old_row['axis-1'],
                          old_row['A0'],
                          old_row['A0_Err'],
                          old_row['A1'],
                          old_row['A1_Err']])

#----------------------------------------------------------------------------------------

AlgorithmFactory.subscribe(VesuvioLt0Resolution)
