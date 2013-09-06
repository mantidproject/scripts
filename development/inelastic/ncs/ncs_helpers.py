##
## Just helpers for now to keep alot of the unimportant details out of the way
## This will all be refactored when it gets incorporated properly
##
import numpy as np

#----------------------------------------------------------------------------------------
def preprocess(data_ws, options):
    """Runs a set of preprocessing steps on the workspace
    before fitting based on the given options
    @param data_ws :: The workspace containing the data to process
    @param fit_options :: An object of type ReductionOptions containing
                          the required parameters
    @returns The result workspace
    """
    if options.has_been_set("smooth_points"):
        pts = options.smooth_points
        if pts > 0:
            print "Smoothing data using %d neighbours" % pts
            from mantid.simpleapi import SmoothData
            data_ws = SmoothData(InputWorkspace=data_ws, OutputWorkspace=data_ws, NPoints=[pts])
        else:
            raise ValueError("Invalid number of points for smoothing: '%d'."
                             " Value must be greater than zero. Set to None to turn off smoothing" % pts)

    if options.has_been_set("bad_data_error"):
        mask_data(data_ws, options.bad_data_error)
        
    return data_ws

#----------------------------------------------------------------------------------------

def mask_data(workspace, error_threshold, verbose=False):
    """
        Masks data points on the workspace that have an error above
        the given error threshold
        @param workspace :: The workspace containing the data. The mask will occur place
        @param error_threshold :: Above this value, the data point will be masked
                                  using MaskBins
        @param verbose :: If true then each bin that is masked is printed to the screen
    """
    from mantid.simpleapi import MaskBins
    print "Masking data with errors > %f" % error_threshold

    # The data shouldn't be too big, clone it to numpy and use its search capabilities
    errors = workspace.extractE()
    if len(errors.shape) != 2:
        raise RuntimeError("Expected 2D array of errors, found %dD" % len(errors.shape))

    indices = np.where(errors > error_threshold) # Indices in 2D matrix where errors above threshold
    # The output is a tuple of 2 arrays where the indices from each array are paired to 
    # give the correct inndex in the original 2D array.
    x_data = workspace.readX(0)

    for ws_index, bin_index in zip(indices[0],indices[1]):
        xmin, xmax = x_data[bin_index],x_data[bin_index+1]
        if verbose:
            print "Masking index %d between [%f,%f]" % (ws_index,xmin,xmax)
            MaskBins(InputWorkspace=workspace, SpectraList=[int(ws_index)],
                     XMin=xmin,XMax=xmax)

#----------------------------------------------------------------------------------------

def run_fitting(data_ws, fit_options):
    """
        Run the Fit algorithm with the given options on the input data
        @param data_ws :: The workspace containing the data to fit too
        @param fit_options :: An object of type FitOptions containing
                              the parameters 
    """
    from mantid.simpleapi import RenameWorkspace, mtd
    
    # Move to validate member
    # If FSE free, then sears flags should be off
    if fit_options.k_free:
        fit_options.sears_flag = 0

    c_free = fit_options.hermite_coeffs
    terms = ""
    for index, value in enumerate(c_free):
        if value == 1:
            terms += "H_%d " % (2*index)
    print "Hermite expansion will contain terms=%s" % terms
    print "Sears flag=%d" % (fit_options.sears_flag)
    if fit_options.has_been_set("background_order"):
        print "Including background using Chebyshev polynomial of order=%d" % fit_options.background_order

    #### Run fitting first time using constraints matrix to reduce active parameter set ######
    function_str = fit_options.create_function_str()
    constraints = fit_options.create_constraints_str()
    ties = fit_options.create_ties_str()
    _do_fit(function_str, data_ws, constraints, ties, max_iter=5000)
    
    #### Run fitting first time using constraints matrix to reduce active parameter set ######
    params_ws = mtd["__fit_Parameters"]
    param_values = {}
    for row in params_ws:
        values = row.values()
        param_values[values[0]] = values[1]

    function_str = fit_options.create_function_str(param_values)
    reduced_chi_square = _do_fit(function_str, data_ws, constraints, ties, max_iter=1)

    ws_prefix = "__fit"
    fit_suffixes = ('_Parameters','_NormalisedCovarianceMatrix','_Workspace')
    for suffix in fit_suffixes:
        RenameWorkspace(InputWorkspace=ws_prefix + suffix,OutputWorkspace="raw" + suffix)
    
    return reduced_chi_square, mtd["raw_Parameters"]

def _do_fit(function_str, data_ws, constraints, ties, max_iter):
    from mantid.simpleapi import Fit

    results = Fit(function_str,data_ws,Ties=ties,Constraints=constraints,Output="__fit",
                  CreateOutput=True,OutputCompositeMembers=True,MaxIterations=max_iter, 
                  Minimizer="Levenberg-Marquardt,AbsError=1e-08,RelError=1e-08")
    return results[1] # reduced chi-squared

#----------------------------------------------------------------------------------------
class ReductionOptions(object):
    """Holds all of the parameters for the reduction & fit"""

    # Main options
    _options = {}
    
    def __init__(self):
        options = {
          "smooth_points":None,
          "bad_data_error":None,
          "background_order":None,
          "hermite_coeffs":None,
          "k_free":None,
          "sears_flag":None,
          "masses" : None,
          "widths":None,
          "widths_lower":None,
          "widths_upper":None,
          "constraints":None,
          "workspace_index":None,
          "output_prefix":None
         }
        object.__setattr__(self, "_options", options)

    def has_been_set(self, name):
        """Returns true if the given option has been set by the user
        """
        return self.__getattr__(name) is not None

    def validate_options(self):
        """Checks that the current options
        are consistent"""
        pass
        
    def create_function_str(self, param_values=None):
        """
            Creates the function string to pass to fit
        
            @param param_values :: A dict of key/values specifying the parameter values
            that have already been calculated. If None then the ComptonScatteringCountRate
            function, along with the constraints matrix, is used rather than the standard CompositeFunction.
            It is assumed that the standard CompositeFunction is used when running the fit for a
            second time to compute the errors with everything free
        """
        def to_space_sep_str(collection):
            _str = ""
            for item in collection:
                _str += " " + str(item)
            return _str.lstrip()
        hermite_str = to_space_sep_str(self.hermite_coeffs)
        ws_index = self.workspace_index

        all_free = (param_values is not None)
        
        if all_free:
            function_str = "composite=CompositeFunction,NumDeriv=1;"
        else:
            function_str = "composite=ComptonScatteringCountRate,NumDeriv=1%s;"
            matrix_str = self.create_matrix_string(self.constraints)
            if matrix_str == "":
                function_str = function_str % ""
            else:
                function_str = function_str % (",IntensityConstraints=" + matrix_str)

        # Currently, it is assumed that the first mass is always  proton/deuterium
        for index, mass in enumerate(self.masses):
            function_str += "name=%s,%s;"
            par_value_prefix = "f%d." % (index)
            if all_free:
                width = param_values[par_value_prefix + "Width"]
            else:
                widths = self.widths[index]
                if hasattr(widths, "__len__"):
                    width = widths[1]
                else:
                    width = widths

            if index == 0:
                func_name = "GramCharlierComptonProfile"
                params = "WorkspaceIndex=%d,Mass=%f,HermiteCoeffs=%s,Width=%f" \
                  % (ws_index,mass,hermite_str,width)
                if all_free:
                    par_names = ["FSECoeff"]
                    for i,c in enumerate(self.hermite_coeffs):
                        if c > 0:
                            par_names.append("C_%d" % (2*i))
                    for par_name in par_names:
                        params += ",%s=%f" % (par_name,param_values[par_value_prefix + par_name])
            else:
                func_name = "GaussianComptonProfile"
                params = "WorkspaceIndex=%d,Mass=%f,Width=%f" % (ws_index,mass,width)
                if all_free:
                    par_name = "Intensity"
                    params += ",%s=%f" % (par_name,param_values[par_value_prefix + par_name])
            function_str = function_str % (func_name, params)
        
        return function_str.rstrip(";")

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
        for index, mass in enumerate(self.masses):
            # Constraints
            func_index = index
            par_name = "f%d.Width" % func_index
            widths = self.widths[index]
            if hasattr(widths, "__len__"):
                constraints += "%f < %s < %f," % (widths[0], par_name, widths[2])

        return constraints.rstrip(",")
   
    def create_ties_str(self):
        """Returns the string of ties for this Fit
        """
        
        ties = ""
        # Widths
        for index, mass in enumerate(self.masses):
            func_index = index
            par_name = "f%d.Width" % func_index
            widths = self.widths[index]
            if not hasattr(widths, "__len__"):
                # Fixed width
                ties += "%s=%f," % (par_name,widths)
            
        ## FSE constraint
        if not self.k_free:
            ties += "%s=%s"
            par_name = "f0.FSECoeff"
            if self.sears_flag == 1:
                value = "f0.Width*sqrt(2)/12"
            else:
                value = "0"
            ties = ties % (par_name, value)

        return ties.rstrip(",")

    def __setattr__(self, name, value):
        """Only allows those attributes that we have defined. 
        Raises a AttributeError if the name is unknown
        """
        if name in self._options:
            # Make sure the constraints tuple is still a tuple
            # If first element looks like sequence we have more than 1 constraint
            # or the tuple was created with a trailing comma and is 1 in size.
            # If the first element is not a sequence a single-element tuple was created without a 
            # trailing comma so the tuple was stripped
            if name == "constraints" and (value is not None) \
                    and (len(value) > 0 and not hasattr(value[0], "__len")):
                    value = (value,) # note trailing comma
            self._options[name] = value
        else:
            raise AttributeError("Unknown attribute %s. "
                            "Allowed names (%s)" % (name, str(self._options.keys())))


            
    def __getattr__(self, name):
        """Return the named attribute's value.
        Raises a NameError if the name is unknown
        """
        try:
            return self._options[name]
        except KeyError:
            raise AttributeError("Unknown attribute %s. "
                            "Allowed names (%s)" % (name, str(self._options.keys())))
            
    def __str__(self):
        """Returns a string representation of the object
        """
        self.generate_function_str()
    
#----------------------------------------------------------------------------------------
        
    
