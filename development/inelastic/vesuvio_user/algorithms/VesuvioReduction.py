# pylint: disable=no-init
import types

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

        self.declareProperty(FileProperty("IPFilename","", action=FileAction.Load,
                                          extensions=["dat"]),
                             doc="An optional IP file. If provided the values are used to correct "
                                 "the default instrument values and attach the t0 values to each "
                                 "detector")

        float_length_validator = FloatArrayLengthValidator()
        float_length_validator.setLengthMin(1)
        self.declareProperty(FloatArrayProperty("Masses", float_length_validator),
                             doc="Mass values for fitting")
        self.declareProperty(FloatArrayProperty("FixedWidths", float_length_validator),
                             doc="Fixed widths for each mass. If the width varies then specify -1 in that "
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

        self.declareProperty("FitMode", "bank", StringListValidator(list(_FIT_MODES)),
                             doc="Fit either bank-by-bank or detector-by-detector")

        self.declareProperty("DifferenceMode", "single", StringListValidator(list(_DIFF_MODES)),
                             doc="The difference option. Valid values: %s" % str(_DIFF_MODES))

        # Outputs
        self.declareProperty(WorkspaceProperty("FittedWorkspace", "", Direction.Output),
                             doc="The name of the fitted workspaces.")
        self.declareProperty(WorkspaceProperty("FittedParameters", "", Direction.Output),
                             doc="The name of the fitted parameter workspaces.")


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
        mass_values = self.getProperty("Masses").value
        profiles = self.getProperty("MassProfiles").value
        fixed_widths = self.getProperty("FixedWidths").value
        width_ranges = self.getProperty("WidthConstraints").value

        fit_opts = _FitOptions()
        fit_opts.masses = []
        for index, (mass, profile, fixed_width) in enumerate(zip(mass_values, profiles, fixed_widths)):
            fit_prop = dict()
            fit_prop["value"] = mass
            fit_prop["function"] = profile
            if profile == "GramCharlier":
                fit_prop["hermite_coeffs"] = [1, 0, 0]
                fit_prop["k_free"] = True
                fit_prop["sears_flag"] = 1
            if fixed_width > 0.0:
                fit_prop["widths"] = fixed_width
            else:
                fit_prop["widths"] = width_ranges[3*index:3*index + 3]
            # save property
            fit_opts.masses.append(fit_prop)

        # TEMPORARY!!!!!!!!!!!
        fit_opts.constraints = ([0, 1, 0, -4],)

        fit_opts.workspace_index = 0
        return self._run_fit(tof_data, fit_opts)

    # -----------------------------------------------------------------------------------------
    def _run_fit(self, tof_ws, fit_options):
        """
            Run the Fit algorithm with the given options on the input data
            @param data_ws :: The workspace containing the data to fit too
            @param fit_options :: An object of type FitOptions containing
                                  the parameters
        """
        if fit_options.global_fit:
            raise RuntimeError("Global fit not implemented yet")
        else:
            return self._run_fit_impl(tof_ws, fit_options, simulation=False)

    def _run_fit_impl(self, data_ws, fit_options, simulation=False):
        """
            Run the Fit algorithm with the given options on the input data
            @param data_ws :: The workspace containing the data to fit too
            @param fit_options :: An object of type FitOptions containing
                                  the parameters
            @param simulation :: If true then only a single iteration is run
        """
        if simulation:
            # Just set what we have been given
            param_values = fit_options.create_param_table()
            constraints = None
            ties = None
        else:
            # Run fitting first time using constraints matrix to reduce active parameter set
            function_str = fit_options.create_function_str()
            constraints = fit_options.create_constraints_str()
            ties = fit_options.create_ties_str()
            results = self._do_fit(function_str, data_ws, fit_options.workspace_index, constraints, ties,
                                   max_iter=5000)

            # Run second time using standard CompositeFunction & no constraints matrix to
            # calculate correct reduced chisq
            param_values = results[1]

        function_str = fit_options.create_function_str(param_values)
        max_iter = 0 if simulation else 1
        return self._do_fit(function_str, data_ws, fit_options.workspace_index, constraints,
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
            return (outputs)

# -----------------------------------------------------------------------------------------
# Helper functions that don't need access to class attributes
# -----------------------------------------------------------------------------------------
class _FitOptions(object):
    """Holds all of the parameters for the reduction & fit"""

    def __init__(self):
        self.smooth_points = None
        self.bad_data_error = None
        self.background_function = None
        self.background_order = None
        self.masses = None
        self.constraints = None
        self.workspace_index = None
        self.output_prefix = None
        self.global_fit = False

    def has_been_set(self, name):
        """Returns true if the given option has been set by the user
        """
        return getattr(self, name) is not None

    def create_param_table(self):
        """Returns a dictionary of parameters that can be used
        to simulate a fit
        @todo cut down copy-and-paste here!
        """
        param_values = {}
        for index, mass_info in enumerate(self.masses):
            par_value_prefix = "f%d." % (index)

            widths = mass_info['widths']
            if hasattr(widths, "__len__"):
                width = widths[1]
            else:
                width = widths
            param_values[par_value_prefix + "Width"] = width

            func_type = mass_info['function']
            if func_type == "GramCharlier":
                if 'FSECoeff' in mass_info:
                    fse = mass_info['FSECoeff']
                else:
                    fse = 0
                param_values[par_value_prefix + "FSECoeff"] = fse
                hermite = mass_info['hermite_coeffs']
                for i,c in enumerate(hermite):
                    if c != 0:
                        param_values["%sC_%d" % (par_value_prefix,2*i)] = c
            elif func_type == "Gaussian":
                param_values[par_value_prefix + "Intensity"] = mass_info["Intensity"]

        return param_values

    # -------------------------------------------------------------------------------------------------------------

    def create_function_str(self, param_values=None):
        """
            Creates the function string to pass to fit

            @param param_values :: A dict/tableworkspace of key/values specifying the parameter values
            that have already been calculated. If None then the ComptonScatteringCountRate
            function, along with the constraints matrix, is used rather than the standard CompositeFunction.
            It is assumed that the standard CompositeFunction is used when running the fit for a
            second time to compute the errors with everything free
        """
        import re
        all_free = (param_values is not None)
        if isinstance(param_values, ITableWorkspace):
            params_ws = param_values
            param_values = {}
            for row in params_ws:
                values = row.values()
                param_name = values[0]
                if self.global_fit:
                    param_name = re.sub('^f\d+\.','',param_name)
                param_values[param_name] = values[1]

        if all_free:
            function_str = "composite=CompositeFunction,NumDeriv=1;"
        else:
            function_str = "composite=ComptonScatteringCountRate,NumDeriv=1%s;"
            matrix_str = self.create_matrix_string(self.constraints)
            if matrix_str == "":
                function_str = function_str % ""
            else:
                function_str = function_str % (",IntensityConstraints=" + matrix_str)

        for index, mass_info in enumerate(self.masses):
            par_value_prefix = "f%d." % (index)
            try:
                func_type = mass_info['function']
            except KeyError:
                raise RuntimeError("Mass details does not contain a function definition")

            if func_type == "GramCharlier":
                local_function_str = \
                    self.create_gram_charlier_function(mass_info, all_free, param_values, par_value_prefix)
            elif func_type == "Gaussian":
                func_name = "GaussianComptonProfile"
                mass = mass_info['value']
                if all_free:
                    width = param_values[par_value_prefix + "Width"]
                else:
                    widths = mass_info['widths']
                    if hasattr(widths, "__len__"):
                        width = widths[1]
                    else:
                        width = widths

                params = "Mass=%f,Width=%f" % (mass,width)
                if all_free:
                    param_name = "Intensity"
                    params += ",%s=%f" % (param_name,param_values[par_value_prefix + param_name])
                local_function_str = "name=%s,%s;" % (func_name,params)
            else:
                raise ValueError("Unknown function type: '%s'. Allowed values are: GramCharlier, Gaussian")
            function_str += local_function_str

        # Add on a background polynomial if requested
        if self.has_been_set("background_order"):
            if not isinstance(self.background_order, types.IntType):
                raise RuntimeError("background_order parameter should be an integer, found '%s'" % type(self.background_order))
            if self.has_been_set("background_function"):
                if not isinstance(self.background_function, types.StringType):
                    raise RuntimeError("background_function parameter should be a string, found '%s'" % type(self.background_function))
                background_func = self.background_function
            else:
                background_func = self._defaults["background_function"]
            background_str = "name=%s,n=%d" % (background_func,self.background_order)
            if all_free:
                func_index = len(self.masses)
                for power in range (0,self.background_order+1):
                    param_name = 'A%d' % (power)
                    comp_par_name = 'f%d.%s' % (func_index,param_name)
                    background_str += ",%s=%f" % (param_name,param_values[comp_par_name])
            function_str += "%s" % background_str.rstrip(",")

        return function_str.rstrip(";")

    def create_gram_charlier_function(self, mass_info, all_free, param_values, par_value_prefix):
        """
        Adds a GramCharlierComptonProfile to the function str and returns
        the updated version
        @param mass_info Dictionary of mass details
        @param all_free True if all should be left free and initial parameters set (see param_values)
        @param values Dictionary of named parameter values from a previous fit
        @param par_value_prefix The prefix within the composite for this function, i.e f0.
        @returns the function string
        """
        func_name = "GramCharlierComptonProfile"
        mass = mass_info['value']
        if all_free:
            width = param_values[par_value_prefix + "Width"]
        else:
            widths = mass_info['widths']
            if hasattr(widths, "__len__"):
                width = widths[1]
            else:
                width = widths

        def to_space_sep_str(collection):
            _str = ""
            for item in collection:
                _str += " " + str(item)
            return _str.lstrip()
        hermite_coeffs = mass_info['hermite_coeffs']
        hermite_str = to_space_sep_str(hermite_coeffs)

        params = "Mass=%f,HermiteCoeffs=%s,Width=%f" % (mass,hermite_str,width)
        if all_free:
            par_names = ["FSECoeff"]
            for i,c in enumerate(hermite_coeffs):
                if c > 0:
                    par_names.append("C_%d" % (2*i))
            for par_name in par_names:
                params += ",%s=%f" % (par_name,param_values[par_value_prefix + par_name])

        return "name=%s,%s;" % (func_name, params)

    def create_matrix_string(self, constraints_tuple):
        """Returns a string for the value of the Matrix of intensity
        constraint values
        """
        if self.constraints is None or len(self.constraints) == 0:
            return ""

        if hasattr(self.constraints[0], "__len__"):
            nrows = len(self.constraints)
            ncols = len(self.constraints[0])
        else:
            nrows = 1
            # without trailing comma a single-element tuple is automatically
            # converted to just be the element
            ncols = len(self.constraints)
            # put back in sequence

        matrix_str = "\"Matrix(%d|%d)%s\""
        values = ""
        for row in self.constraints:
            for val in row:
                values += "%f|" % val
        values = values.rstrip("|")
        matrix_str = matrix_str % (nrows, ncols, values)
        return matrix_str

    def create_constraints_str(self):
        """Returns the string of constraints for this Fit
        """
        constraints = ""
        for index, mass_info in enumerate(self.masses):
            # Constraints
            func_index = index
            par_name = "f%d.Width" % func_index
            widths = mass_info['widths']
            if hasattr(widths, "__len__"):
                constraints += "%f < %s < %f," % (widths[0], par_name, widths[2])

        return constraints.rstrip(",")

    def create_ties_str(self):
        """Returns the string of ties for this Fit
        """

        ties = ""
        # Widths
        for index, mass_info in enumerate(self.masses):
            func_index = index
            par_value_prefix = "f%d." % (func_index)
            par_name = "%sWidth" % par_value_prefix
            widths = mass_info['widths']
            if not hasattr(widths, "__len__"):
                # Fixed width
                ties += "%s=%f," % (par_name,widths)

            func_type = mass_info['function']
            if func_type == "GramCharlier":
                if 'k_free' not in mass_info:
                    raise RuntimeError("GramCharlier requested for mass %d but no k_free argument was found" % (index+1))
                k_free = mass_info['k_free']
                ## FSE constraint
                if k_free:
                    continue

                if 'sears_flag' not in mass_info:
                    raise RuntimeError("Fixed k requested for mass %d but no sears_flag argument was found" % (index+1))
                sears_flag = mass_info['sears_flag']
                par_name = "%sFSECoeff" % par_value_prefix
                if sears_flag == 1:
                    value = "%sWidth*%s" % (par_value_prefix,"sqrt(2)/12")#math.sqrt(2.)/12.0)
                else:
                    value = "0"
                ties += "%s=%s," % (par_name, value)

        return ties.rstrip(",")

    def create_global_function_str(self, n, param_values=None):
        """
            Creates the function string to pass to fit for a multi-dataset (global) fitting

            @param n :: A number of datasets (spectra) to be fitted simultaneously.

            @param param_values :: A dict/tableworkspace of key/values specifying the parameter values
            that have already been calculated. If None then the ComptonScatteringCountRate
            function, along with the constraints matrix, is used rather than the standard CompositeFunction.
            It is assumed that the standard CompositeFunction is used when running the fit for a
            second time to compute the errors with everything free
        """

        # create a local function to fit a single spectrum
        f = self.create_function_str(param_values)
        # insert an attribute telling the function which spectrum it should be applied to
        i = f.index(';')
        # $domains=i means "function index == workspace index"
        fun_str = f[:i] + ',$domains=i' + f[i:]

        # append the constrints and ties within the local function
        fun_str += ';constraints=(' + self.create_constraints_str() + ')'
        ties = self.create_ties_str()
        if len(ties) > 0:
            fun_str += ';ties=(' + ties + ')'

        # initialize a string for composing the global ties
        global_ties = 'f0.f0.Width'
        # build the multi-dataset function by joining local functions of the same type
        global_fun_str = 'composite=MultiDomainFunction'
        for i in range(n):
            global_fun_str += ';(' + fun_str + ')'
            if i > 0:
                global_ties = 'f' + str(i) + '.f0.Width=' + global_ties
        # add the global ties
        global_fun_str += ';ties=(' + global_ties + ')'

        return global_fun_str

    def __str__(self):
        """Returns a string representation of the object
        """
        self.generate_function_str()



# -----------------------------------------------------------------------------------------
AlgorithmFactory.subscribe(VesuvioReduction)