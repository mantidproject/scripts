# pylint: disable=no-init
from mantid.kernel import *
from mantid.api import *

from vesuvio.algorithms.base import VesuvioBase


class VesuvioPreprocess(VesuvioBase):

    def summary(self):
        return "Apply preprocessing steps to loaded vesuvio data"

    def PyInit(self):
        # Inputs
        self.declareProperty(MatrixWorkspaceProperty("InputWorkspace", "", Direction.Input),
                             doc="Input TOF workspace from LoadVesuvio")

        smooth_opts = ["Neighbour", "None"]
        self.declareProperty("Smoothing", smooth_opts[0], StringListValidator(smooth_opts),
                             doc="Defines the smoothing method.")
        self.declareProperty("SmoothingOptions", "NPoints=3",
                             doc="Override the default smoothing options")

        # Outputs
        self.declareProperty(MatrixWorkspaceProperty("OutputWorkspace", "", Direction.Output),
                             doc="The name of the output workspace")

    def validateInputs(self):
        errors = dict()
        smoothing = self.getProperty("Smoothing").value
        if smoothing == "Neighbour":
            options = self.getProperty("SmoothingOptions").value
            if not options.startswith("NPoints="):
                errors["SmoothingOptions"] = "Invalid value for smoothing option. It must begin the format NPoints=3"

        return errors

    def PyExec(self):
        data = self._apply_smoothing(self.getProperty("InputWorkspace").value)

        self.setProperty("OutputWorkspace", data)


    def _apply_smoothing(self, data):
        smoothing = self.getProperty("Smoothing").value
        options = self.getProperty("SmoothingOptions").value
        if smoothing == "None":
            return data
        elif smoothing == "Neighbour":
            npts = int(options[-1])
            return self._execute_child_alg("SmoothData", InputWorkspace=data, OutputWorkspace=data,
                                           NPoints=npts)


# -----------------------------------------------------------------------------------------
AlgorithmFactory.subscribe(VesuvioPreprocess)
