# pylint: disable=no-init
from mantid.kernel import *
from mantid.api import *

from vesuvio.fitting import parse_fit_options

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

        self.declareProperty(FileProperty("IPFilename","", action=FileAction.Load,
                                          extensions=["dat"]),
                             doc="An optional IP file. If provided the values are used to correct "
                                 "the default instrument values and attach the t0 values to each "
                                 "detector")

        float_length_validator = FloatArrayLengthValidator()
        float_length_validator.setLengthMin(1)
        self.declareProperty(FloatArrayProperty("Masses", float_length_validator),
                             doc="Mass values for fitting")

        self.declareProperty("MassProfiles", "", StringMandatoryValidator(),
                             doc="Functions used to approximate mass profile. "
                                 "The format is Function1,param1=val1,param2=val2;Function2,param3=val3,param4=val4")

        # ----- Optional ------
        self.declareProperty("Spectra", "forward", StringMandatoryValidator(),
                             doc="The spectrum numbers to load. "
                                 "A dash will load a range and a semicolon delimits spectra to sum. "
                                 "The keyword forward and backward indicate either all of forward/backward")

        self.declareProperty("FitMode", "bank", StringListValidator(list(_FIT_MODES)),
                             doc="Fit either bank-by-bank or detector-by-detector")

        self.declareProperty("DifferenceMode", "single", StringListValidator(list(_DIFF_MODES)),
                             doc="The difference option. Valid values: %s" % str(_DIFF_MODES))

        # Outputs
        self.declareProperty(WorkspaceProperty("FittedWorkspace", "", Direction.Output),
                             doc="The name of the fitted workspaces.")
        self.declareProperty(WorkspaceProperty("FittedParameters", "", Direction.Output),
                             doc="The name of the fitted parameter workspaces.")

    def PyExec(self):
        """Run the processing"""
        loaded_data = self._load_and_crop_data()
        results = self._fit_tof(loaded_data)

        self.setProperty("FittedWorkspace", results[2])
        self.setProperty("FittedParameters", results[1])

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
        full_range = self._execute_child_alg("LoadVesuvio", **kwargs)
        return self._execute_child_alg("CropWorkspace", InputWorkspace=full_range,
                                       XMin=_TOF_RANGE[0], XMax=_TOF_RANGE[1])

    def _fit_tof(self, tof_data):
        """
        Runs a fit against the loaded data
        """
        fit_opts = parse_fit_options(self.getProperty("Masses").value,
                                     self.getProperty("MassProfiles").value)

        #!!!! TEMPORARY!!!
        workspace_index = 0
        return self._run_fit(tof_data, workspace_index, fit_opts)

    # -----------------------------------------------------------------------------------------
    def _run_fit(self, tof_ws, workspace_index, fit_options):
        """
            Run the Fit algorithm with the given options on the input data
            @param data_ws :: The workspace containing the data to fit too
            @param fit_options :: An object of type FitOptions containing
                                  the parameters
        """
        if fit_options.global_fit:
            raise RuntimeError("Global fit not implemented yet")
        else:
            return self._run_fit_impl(tof_ws, workspace_index, fit_options, simulation=False)

    def _run_fit_impl(self, data_ws, workspace_index, fit_options, simulation=False):
        """
            Run the Fit algorithm with the given options on the input data
            @param data_ws :: The workspace containing the data to fit too
            @param workspace_index :: The spectra to fit against
            @param fit_options :: An object of type FittingOptions containing
                                  the parameters
            @param simulation :: If true then only a single iteration is run
        """
        if simulation:
            raise RuntimeError("Simulation not implemented yet")
        else:
            # Run fitting first time using constraints matrix to reduce active parameter set
            function_str = fit_options.create_function_str()
            constraints = fit_options.create_constraints_str()
            ties = fit_options.create_ties_str()
            results = self._do_fit(function_str, data_ws, workspace_index, constraints, ties,
                                   max_iter=5000)

            # Run second time using standard CompositeFunction & no constraints matrix to
            # calculate correct reduced chisq
            param_values = TableWorkspaceDictionaryFacade(results[1])

        function_str = fit_options.create_function_str(param_values)
        max_iter = 0 if simulation else 1
        return self._do_fit(function_str, data_ws, workspace_index, constraints,
                            ties, max_iter=max_iter)

    # ----------------------------------------------------------------------------------------

    def _do_fit(self, function_str, data_ws, index, constraints, ties, max_iter):

        # The tof data is required to be in seconds for the fitting
        # in order to re-use the standard Mantid Polynomial function. This polynomial simply
        # accepts the data "as is" in the workspace so if it is in microseconds then the
        # we would have to either implement a another wrapper to translate or write another
        # Polynomial.
        # The simplest option is to put the data in seconds here and then put it back afterward
        # TODO: This is quite wasteful. Come up with something better!
        data_ws = self._execute_child_alg("ScaleX", InputWorkspace=data_ws, OutputWorkspace=data_ws,
                                          Operation='Multiply',Factor=1e-06)

        outputs = self._execute_child_alg("Fit", Function=function_str, InputWorkspace=data_ws,
                                          WorkspaceIndex=index, Ties=ties,
                                          Constraints=constraints,
                                          CreateOutput=True, OutputCompositeMembers=True, MaxIterations=max_iter,
                                          Minimizer="Levenberg-Marquardt,AbsError=1e-08,RelError=1e-08")

        reduced_chi_square, params, fitted_data = outputs[1], outputs[3], outputs[4]
        fitted_data = self._execute_child_alg("ScaleX", InputWorkspace=fitted_data, OutputWorkspace=fitted_data,
                                              Operation='Multiply', Factor=1e06)
        data_ws = self._execute_child_alg("ScaleX", InputWorkspace=data_ws, OutputWorkspace=data_ws,
                                          Operation='Multiply',Factor=1e06)

        return reduced_chi_square, params, fitted_data

    # ----------------------------------------------------------------------------------------

    def _execute_child_alg(self, name, **kwargs):
        alg = self.createChildAlgorithm(name)
        for name, value in kwargs.iteritems():
            alg.setProperty(name, value)
        alg.execute()
        outputs = list()
        for name in alg.outputProperties():
            outputs.append(alg.getProperty(name).value)
        if len(outputs) == 1:
            return outputs[0]
        else:
            return tuple(outputs)

# -----------------------------------------------------------------------------------------
# Helper to translate from an table workspace to a dictionary. Should be on the workspace
# really ...
# -----------------------------------------------------------------------------------------
class TableWorkspaceDictionaryFacade(object):
    """
    Allows an underlying table workspace to be treated like a read-only dictionary
    """

    def __init__(self, held_object):
        self._table_ws = held_object

    def __getitem__(self, item):
        for row in self._table_ws:
            if row['Name'] == item:
                return row['Value']
        #endfor
        raise KeyError(str(item))

# -----------------------------------------------------------------------------------------
AlgorithmFactory.subscribe(VesuvioReduction)