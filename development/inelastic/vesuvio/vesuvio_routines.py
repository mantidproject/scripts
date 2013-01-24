##
## Functions to aid VESUVIO reduction
##
import numpy as np

def preprocess(data_ws, options):
    """Runs a set of preprocessing steps on the workspace
    before fitting based on the given options
    @param data_ws :: The workspace containing the data to process
    @param fit_options :: An object of type ReducitonOptions containing
                          the required parameters
    @returns The result workspace
    """
    if options.has_been_set("smooth_points"):
        pts = options.smooth_points
        print "Smoothing data using %d neighbours" % pts
        from mantid.simpleapi import SmoothData
        data = SmoothData(InputWorkspace=data_ws, OutputWorkspace=data_ws, NPoints=[pts])

    if options.has_been_set("bad_data_error"):
        mask_data(data, options.bad_data_error)
        
    return data

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
    from mantid.simpleapi import Fit, RenameWorkspace
    
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

    
    function_str = fit_options.create_function_str()
    constraints = fit_options.create_constraints_str()
    ties = fit_options.create_ties_str()

    ws_prefix = "__fit"
    Fit(function_str,data_ws,Ties=ties,Constraints=constraints,Output=ws_prefix,
        CreateOutput=True,OutputCompositeMembers=True,MaxIterations=5000)
    
    fit_suffixes = ('_Parameters','_NormalisedCovarianceMatrix','_Workspace')
    for suffix in fit_suffixes:
        RenameWorkspace(InputWorkspace=ws_prefix + suffix,OutputWorkspace=fit_options.output_prefix + suffix)

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
        
    def create_function_str(self):
        """
        Creates the function string to pass to fit
        """
        def to_space_sep_str(collection):
            _str = ""
            for item in collection:
                _str += " " + str(item)
            return _str.lstrip()
        hermite_str = to_space_sep_str(self.hermite_coeffs)
        ws_index = self.workspace_index

        function_str = "composite=CompositeFunction,NumDeriv=1;"
        # Currently, it is assumed that the first mass is always  proton/deuterium
        for index, mass in enumerate(self.masses):
            # Main function
            function_str += "name=%s,%s;"
            if index == 0:
                func_name = "GramCharlierComptonProfile"
                params = "WorkspaceIndex=%d,Mass=%f,HermiteCoeffs=%s" \
                  % (ws_index,mass,hermite_str)
            else:
                func_name = "GaussianComptonProfile"
                params = "WorkspaceIndex=%d,Mass=%f" % (ws_index,mass)
            function_str = function_str % (func_name, params)
        
        return function_str.rstrip(";")

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

        ### Stoichiometry/Intensity constraints
        nmasses = len(self.masses)
        c_free = self.hermite_coeffs
        n_c = len(c_free)
        constraints = self.constraints
        # Make sure it is a list of lists
        if type(constraints[0]) != list:
            constraints = [constraints]
        
        for aeq_i in constraints:
            # The first value in the aeq list relates to the first mass. We need to expand the
            # matrix to incorporate the hermite polynomial coeffs that are used instead of the Gaussian
            intensity_pars = ["f%d.Intensity" % index for index in range(1,nmasses)]
            firstmass_aeq = aeq_i[0]
            aeq_i = aeq_i[1:]
            for back_index, coeff in enumerate(reversed(c_free)):
                if coeff == 1:
                    aeq_i.insert(0, firstmass_aeq) #prepend
                    intensity_pars.insert(0,"f0.C_%d" % (2*(n_c - 1 - back_index)) ) #Even coefficents
        
            # Make list of tuples (coeff,par_name) that contribute
            lhs_terms = []
            for coeff, par_name in zip(aeq_i, intensity_pars):
                if abs(coeff) > 0.0:
                    lhs_terms.append((coeff,par_name))
            
            # Pick first par name to go on lhs of tie expression.
            # and rearrange so everything else is on the the RHS of aeq_i*x_i = 0
            lhs = lhs_terms[0]
            rhs = lhs_terms[1:]
            denominator = -1.0*lhs[0]
            rhs = map(lambda term: (term[0]/denominator,term[1]), rhs)
            # Now form p=a*b+c*d...
            rhs_str = ""
            for coeff, param in rhs:
                rhs_str += "%+f*%s" % (coeff,param) # +f prints +/- with number
            ties += lhs[1] + "=" + rhs_str.lstrip("+") + ","
            
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

