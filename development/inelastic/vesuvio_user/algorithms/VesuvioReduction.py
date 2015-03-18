# pylint: disable=no-init

from mantid.kernel import *
from mantid.api import *

# Loading difference modes
_DIFF_MODES = ("double", "single")
# Fitting modes
_FIT_MODES = ("bank", "spectrum")
# Spectra ranges
_BACKWARD_SPECTRA = (3, 134)
_BACKWARD_BANKS = ((3, 46), (47, 90), (91, 134))
_FORWARD_SPECTRA = (135, 198)
_FORWARD_BANKS = ((135, 142), (143, 150), (151, 158), (159, 166),
                  (167, 174), (175, 182), (183, 190), (191, 198))
# Crop range (VMS defaults)
_TOF_RANGE = [50, 562]

class VesuvioReduction(DataProcessorAlgorithm):


    def summary(self):
        return "Processes runs for Vesuvio at ISIS"

    def PyInit(self):
        # Inputs
        self.declareProperty("Runs", "", StringMandatoryValidator(),
                     doc="The range of run numbers that should be loaded.")

        float_length_validator = FloatArrayLengthValidator()
        float_length_validator.setLengthMin(1)
        self.declareProperty(FloatArrayProperty("Masses", float_length_validator),
                             doc="Mass values for fitting")
        self.declareProperty(FloatArrayProperty("FixedWidths", float_length_validator),
                             doc="Fixed widths for each mass. If the width varies then specify 0 in that "
                                 "position and use the WidthConstraints property.")
        str_length_validator = StringArrayLengthValidator()
        str_length_validator.setLengthMin(1)
        self.declareProperty(StringArrayProperty("MassProfiles", str_length_validator),
                             doc="Functions used to approximate mass profile. "
                                 "The length should match the number of masses")

        # ----- Optional ------
        self.declareProperty(FloatArrayProperty("WidthConstraints"),
                             doc="Range constraints for widths during the fit in the format min,default,max. "
                                 "Provide these 3 numbers per masses that has a 0 entry for FixedWidths")

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

        self.declareProperty("DifferenceMode", "single", StringListValidator(_DIFF_MODES),
                             doc="The difference option. Valid values: %s" % str(_DIFF_MODES))

        # Outputs
        self.declareProperty(WorkspaceProperty("OutputWorkspace", "", Direction.Output),
                             doc="The name of the output workspace.")


    def validateInputs(self):
        """Cross-check the inputs"""
        messages = dict()

        # Number of masses & functions must equal
        masses = self.getProperty("Masses").value
        profiles = self.getProperty("MassProfiles").value
        if len(profiles) != len(masses):
            messages["MassProfiles"] = "Number of functions must match the number of masses"

        # Number of fixed width values & masses must equal
        fixed_widths = self.getProperty("FixedWidths").value
        if len(fixed_widths) != len(masses):
            messages["FixedWidths"] = "Number of fixed widths must match the number of masses"

        # Width constraints. 3 values per unfixed width
        widths_ranges = self.getProperty("WidthConstraints").value
        nfixed = len(filter(lambda x: x > 0, fixed_widths))
        nwidth_vals = len(widths_ranges)
        expected_nwidth_vals = 3*(len(fixed_widths) - nfixed)
        if nwidth_vals != expected_nwidth_vals:
            messages["WidthConstraints"] = "Constraints should be min,default,max for each unfixed width"

        return messages

    def PyExec(self):
        """Run the processing"""
        loaded_data = self._load_and_crop_data()

        self.setProperty("OutputWorkspace", loaded_data)

    def _load_and_crop_data(self):
        """
           Load the data and return the workspace
        """
        spectra = self.getProperty("Spectra").value
        if spectra == "forward":
            spectra = "{0}-{1}".format(*_FORWARD_SPECTRA)
        elif spectra == "backward":
            spectra = "{0}-{1}".format(*_BACKWARD_SPECTRA)

        if self.getProperty("DifferenceMode").value == "double":
            diff_mode = "DoubleDifference"
        else:
            diff_mode = "SingleDifference"

        kwargs = {"Filename": self.getProperty("Runs").value,
                  "Mode": diff_mode, "InstrumentParFile": self.getProperty("IPFilename").value,
                  "SpectrumList": spectra}
        loader = self._execute_child_alg("LoadVesuvio", **kwargs)
        full_range = loader.getProperty("OutputWorkspace").value
        cropper = self._execute_child_alg("CropWorkspace", InputWorkspace=full_range,
                                          XMin=_TOF_RANGE[0], XMax=_TOF_RANGE[1])
        return cropper.getProperty("OutputWorkspace").value

    def _execute_child_alg(self, name, **kwargs):
        alg = self.createChildAlgorithm(name)
        for name, value in kwargs.iteritems():
            alg.setProperty(name, value)
        alg.execute()
        return alg

# -----------------------------------------------------------------------------------------
AlgorithmFactory.subscribe(VesuvioReduction)