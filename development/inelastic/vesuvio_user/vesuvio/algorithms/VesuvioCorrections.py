# pylint: disable=no-init
from mantid.kernel import *
from mantid.api import *
import numpy as np
from vesuvio.algorithms.base import VesuvioBase


class VesuvioCorrections(VesuvioBase):

    def summary(self):
        return "Apply post fitting steps to vesuvio data"

    def PyInit(self):
        # Inputs
        self.declareProperty(MatrixWorkspaceProperty("InputWorkspace", "", direction=Direction.Input),
                             doc="Input TOF workspace")

        self.declareProperty("GammaBackground", True, direction=Direction.Input,
                             doc="If true, correct for the gamma background")
        self.declareProperty(ITableWorkspaceProperty("FitParameters", "", direction=Direction.Input,
                                                     optional=PropertyMode.Optional),
                             doc="Table containing the calculated fit parameters for the data in the workspace")

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
        pass

# -----------------------------------------------------------------------------------------
AlgorithmFactory.subscribe(VesuvioCorrections)
