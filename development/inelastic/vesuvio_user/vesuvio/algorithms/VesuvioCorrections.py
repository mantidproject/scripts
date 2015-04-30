# pylint: disable=no-init
from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import ExtractSingleSpectrum, CalculateGammaBackground, DeleteWorkspace
import numpy as np
from vesuvio.algorithms.base import VesuvioBase, TableWorkspaceDictionaryFacade
from vesuvio.fitting import parse_fit_options


class VesuvioCorrections(VesuvioBase):

    def summary(self):
        return "Apply post fitting steps to vesuvio data"

    def PyInit(self):
        # Inputs
        self.declareProperty(MatrixWorkspaceProperty("InputWorkspace", "", direction=Direction.Input),
                             doc="Input TOF workspace")

        self.declareProperty("WorkspaceIndex", 0,
                             doc="Index of spectrum to calculate corrections for")

        # Gamma background
        self.declareProperty("GammaBackground", True, direction=Direction.Input,
                             doc="If true, correct for the gamma background")

        self.declareProperty(ITableWorkspaceProperty("FitParameters", "", direction=Direction.Input,
                                                     optional=PropertyMode.Optional),
                             doc="Table containing the calculated fit parameters for the data in the workspace")

        float_length_validator = FloatArrayLengthValidator()
        float_length_validator.setLengthMin(1)
        self.declareProperty(FloatArrayProperty("Masses", float_length_validator),
                             doc="Mass values for fitting")

        self.declareProperty("MassProfiles", "", StringMandatoryValidator(),
                             doc="Functions used to approximate mass profile. "
                                 "The format is function=Function1Name,param1=val1,param2=val2;function=Function2Name,param3=val3,param4=val4")

        self.declareProperty("IntensityConstraints", "",
                             doc="A semi-colon separated list of intensity constraints defined as lists e.g "
                                 "[0,1,0,-4];[1,0,-2,0]")

        # Multiple scattering
        self.declareProperty("MultipleScattering", True, direction=Direction.Input,
                             doc="If true, correct for the effects of multiple scattering")

        # Outputs
        self.declareProperty(MatrixWorkspaceProperty("OutputWorkspace", "", direction=Direction.Output),
                             doc="The name of the output workspace")


    def validateInputs(self):
        errors = dict()

        # The gamma correction requires the table of fitted parameters
        if self.getProperty("GammaBackground").value and self.getProperty("FitParameters").value is None:
            errors["FitParameters"] = "Gamma background correction requires a set of parameters from a fit of the data"

        return errors


    def PyExec(self):
        input_ws = self.getPropertyValue("InputWorkspace")
        spec_idx = self.getProperty("WorkspaceIndex").value
        self._output_ws = self.getPropertyValue("OutputWorkspace")

        ExtractSingleSpectrum(InputWorkspace=input_ws,
                              OutputWorkspace=self._output_ws,
                              WorkspaceIndex=spec_idx)

        if self.getProperty("GammaBackground").value:
            self._gamma_correction()

        if self.getProperty("MultipleScattering").value:
            self._ms_correction()

        self.setProperty("OutputWorkspace", self._output_ws)


    def _gamma_correction(self):
        fit_opts = parse_fit_options(mass_values=self.getProperty("Masses").value,
                                     profile_strs=self.getProperty("MassProfiles").value,
                                     constraints_str=self.getProperty("IntensityConstraints").value)

        params_ws_name = self.getPropertyValue("FitParameters")
        params_dict = TableWorkspaceDictionaryFacade(mtd[params_ws_name])
        func_str = fit_opts.create_function_str(params_dict)

        num_hist = mtd[self._output_ws].getNumberHistograms()

        CalculateGammaBackground(InputWorkspace=self._output_ws,
                                 ComptonFunction=func_str,
                                 BackgroundWorkspace='__background',
                                 CorrectedWorkspace=self._output_ws)
        DeleteWorkspace('__background')


    def _ms_correction(self):
        pass


# -----------------------------------------------------------------------------------------
AlgorithmFactory.subscribe(VesuvioCorrections)
