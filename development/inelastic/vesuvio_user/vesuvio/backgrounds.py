"""Holds classes that define the backgrounds for fitting
"""

# --------------------------------------------------------------------------------
# Background
# --------------------------------------------------------------------------------

class Background(object):
    """Base class"""
    pass

# --------------------------------------------------------------------------------
# Polynomial
# --------------------------------------------------------------------------------

class PolynomialBackground(object):

    cfunction = "Polynomial"

    def __init__(self, order):
        self.order = order

    def create_fit_function_str(self, param_vals=None, param_prefix=""):
        """Creates a string used by the Fit algorithm for this function

        :param param_vals: A table of values for the parameters that override those set
        on the object already. Default=None
        :param param_prefix: A string prefix for the parameter name in the params_vals list
        """
        vals_provided = (param_vals is not None)
        func_str = "name={0},n={1}".format(self.cfunction, str(self.order))

        if vals_provided:
            for power in range (0,self.order+1):
                param_name = 'A{0}'.format(power)
                func_str += ",{0}={1}".format(param_name, param_vals[param_prefix + param_name])

        return func_str
