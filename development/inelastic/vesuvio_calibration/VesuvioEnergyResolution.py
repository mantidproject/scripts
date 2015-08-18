#pylint: disable=no-init

from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *

import math
import numpy as np
import scipy.constants as sc
import scipy.stats as stats

NEUTRON_MASS_AMU = sc.value('neutron mass in u')

BACKSCATTERING = range(3, 135)
FORWARDSCATTERING = range(135, 199)

MODES = ['SingleDifference', 'DoubleDifference', 'ThickDifference']

#----------------------------------------------------------------------------------------

class VesuvioEnergyResolution(PythonAlgorithm):

    _evs_ws = None
    _evs = None
    _sample_pos = None
    _l0_dist = None
    _output_table = None

    def summary(self):
        return "Calculates the resolution of the total flight path and t0 for VESUVIO."

    def PyInit(self):
        # Sample data
        self.declareProperty(StringArrayProperty("Samples", Direction.Input),
                             doc="Sample run numbers to fit peaks to.")

        self.declareProperty(StringArrayProperty("Background", Direction.Input),
                             doc="Run numbers to use as a background.")

        # Resolution params
        self.declareProperty("DeltaTheta", 0.0,
                             doc="Resolution in the scattering angle")

        self.declareProperty("DeltaL0", 0.0,
                             doc="Resolution in the incident flight path length")

        self.declareProperty("DeltaL1", 0.0,
                             doc="Resolution in the final flight path length")

        self.declareProperty("DeltaT0", 0.0,
                             doc="Resolution in the pulse offset time")

        # Optional parameter file
        self.declareProperty(FileProperty("InstrumentParFile", "", action=FileAction.OptionalLoad,
                                          extensions=["dat", "par"]),
                             doc="An optional IP file. If provided the values are used to correct "
                                 "the default instrument values and attach the t0 values to each "
                                 "detector")

        # Output parameter tables
        self.declareProperty(WorkspaceGroupProperty("ResolutionGroup", "", Direction.Output),
                             doc="Resolution values for forward scattering")

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

        # Create mean output table
        self._output_table = CreateEmptyTableWorkspace(OutputWorkspace=self.getPropertyValue("OutputWorkspace"))

        self._output_table.addColumn('str', 'ResolutionParameter')
        self._output_table.addColumn('float', 'Mean')
        self._output_table.addColumn('float', 'StdDev')
        self._output_table.addColumn('float', 'StdErr')
        self._output_table.addColumn('float', 'WeightedMean')
        self._output_table.addColumn('float', 'WeightedErr')

        # Calculate resolution
        res_workspaces = []
        for mode in MODES:
            self._do_fit(mode)

            out_ws = CreateEmptyTableWorkspace(OutputWorkspace='{0}_EnergyResolution'.format(mode))
            res_workspaces.append(out_ws)

            out_ws.addColumn('int', 'Spectrum')
            out_ws.addColumn('float', 'DeltaE1_Lorentz')
            out_ws.addColumn('float', 'DeltaE1_Gauss')

            self._calculate_de1(mtd['{0}_Peaks_Peak_0_Parameters'.format(mode)], out_ws)

            de1_gauss = np.array(out_ws.column('DeltaE1_Gauss'))
            de1_lorentz = np.array(out_ws.column('DeltaE1_Lorentz'))

            de1_gauss = de1_gauss[np.isfinite(de1_gauss)]
            de1_lorentz = de1_lorentz[np.isfinite(de1_lorentz)]

            self._add_param_table_row('{0} Backscattering dE1_Gauss'.format(mode),
                                      de1_gauss[:FORWARDSCATTERING[0]])
            self._add_param_table_row('{0} Backscattering dE1_Lorentz'.format(mode),
                                      de1_lorentz[:FORWARDSCATTERING[0]])
            self._add_param_table_row('{0} Forward scattering dE1_Gauss'.format(mode),
                                      de1_gauss[FORWARDSCATTERING[0]:])
            self._add_param_table_row('{0} Forward scattering dE1_Lorentz'.format(mode),
                                      de1_lorentz[FORWARDSCATTERING[0]:])

        res_grp = GroupWorkspaces(InputWorkspaces=res_workspaces,
                                  OutputWorkspace=self.getPropertyValue('ResolutionGroup'))

        DeleteWorkspace(self._evs_ws)

        self.setProperty("ResolutionGroup", res_grp)
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

    def _do_fit(self, mode):
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

    def _calculate_de1(self, parameters, out_ws):
        """
        Calculates delta-E1 for Lorentzian and Gaussian components.

        @param parameters Fit parameters
        @param out_ws Table workspace to add results to
        """
        spectrum_nos = parameters.column('Spectrum')
        positions = parameters.column('f1.LorentzPos')
        widths_lorentz = parameters.column('f1.LorentzFWHM')
        widths_gauss = parameters.column('f1.GaussianFWHM')

        delta_theta = self.getProperty("DeltaTheta").value
        delta_L0 = self.getProperty("DeltaL0").value
        delta_L1 = self.getProperty("DeltaL1").value
        delta_t0 = self.getProperty("DeltaT0").value

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
                v = delta_e1_gauss**2 - delta_theta**2 - delta_L0**2 - delta_L1**2 - delta_t0**2
                if v < 0:
                    delta_e1_gauss = np.nan
                else:
                    delta_e1_gauss = math.sqrt(v)

            out_ws.addRow([spec, delta_e1_lorentz, delta_e1_gauss])

#----------------------------------------------------------------------------------------

    def _convert_to_energy(self, l1_dist, position, fwhm):
        """
        Converts a peak width to energy.

        @param l1_dist Distance from sample to detector (m)
        @param position Peak position in time of flight (us)
        @param fwhm FWHM in time of flight (us)
        """
        # Calculate energy at the peak position
        d_time = position / sc.micro # us to s
        neutron_v1 = (self._l0_dist + l1_dist) / d_time # ms-1
        peak_e = 0.5 * sc.m_n * neutron_v1**2 # Joule
        peak_e /= sc.value('electron volt') / sc.milli # Joule to meV

        # Calculate HWHM of peak in meV
        energy = peak_e * fwhm / position

        return energy

#----------------------------------------------------------------------------------------

AlgorithmFactory.subscribe(VesuvioEnergyResolution)
