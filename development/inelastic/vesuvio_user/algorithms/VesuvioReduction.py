# pylint: disable=no-init

from mantid.kernel import *
from mantid.api import *

# Loading difference modes
_DIFF_MODES = ["double", "single"]
# Fitting modes
_FIT_MODES = ["bank", "spectrum"]


class VesuvioReduction(DataProcessorAlgorithm):


    def summary(self):
        return "Processes runs for Vesuvio at ISIS"

    def PyInit(self):
        # Inputs
        runs_length_validator = IntArrayLengthValidator(lenmin=1, lenmax=2)
        self.declareProperty(IntArrayProperty("Runs", validator=runs_length_validator),
                     doc="The range of run numbers that should be loaded.")

        float_length_validator = FloatArrayLengthValidator()
        float_length_validator.setLengthMin(1)
        self.declareProperty(FloatArrayProperty("Masses", float_length_validator),
                             doc="Mass values for fitting")
        self.declareProperty(FloatArrayProperty("Widths", float_length_validator),
                             doc="Default widths & limits for masses")
        str_length_validator = StringArrayLengthValidator()
        str_length_validator.setLengthMin(1)
        self.declareProperty(StringArrayProperty("MassProfiles", str_length_validator),
                             doc="Functions used to approximate mass profile. "
                                 "The length should match the number of masses")

        # ----- Optional ------
        self.declareProperty("Spectra", "forward", StringMandatoryValidator(),
                             doc="The spectrum numbers to load. "
                                 "A dash will load a range and a semicolon delimits spectra to sum. "
                                 "The keyword forward and backward indicate either all of forward/backward")

        self.declareProperty("FitMode", "bank", StringListValidator(_FIT_MODES),
                             doc="Fit either bank-by-bank or detector-by-detector")

        self.declareProperty(FileProperty("IPFilename","", action=FileAction.OptionalLoad,
                                          extensions=["dat"]),
                             doc="An optional IP file. If provided the values are used to correct "
                                 "the default instrument values and attach the t0 values to each "
                                 "detector")

        self.declareProperty("SumSpectra", False,
                             doc="If true then the final output is a single spectrum containing "
                                 "the sum of all of the requested spectra. All detector angles/"
                                 "parameters are averaged over the individual inputs")

        self.declareProperty("DifferenceMode", "double", StringListValidator(_DIFF_MODES),
                             doc="The difference option. Valid values: %s" % str(_DIFF_MODES))

        # Outputs
        # self.declareProperty(WorkspaceProperty("OutputWorkspace", "", Direction.Output),
        #                      doc="The name of the output workspace.")


    def validateInputs(self):
        messages = dict()

        # Number of masses & functions must equal
        masses = self.getProperty("Masses").value
        profiles = self.getProperty("MassProfiles").value
        if len(masses) != len(profiles):
            messages["MassProfiles"] = "Number of functions must match the number of masses"

        return messages

    def PyExec(self):
        pass

# -----------------------------------------------------------------------------------------
AlgorithmFactory.subscribe(VesuvioReduction)