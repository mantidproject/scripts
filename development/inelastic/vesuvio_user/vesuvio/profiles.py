"""Holds classes that define the mass profiles
"""
import ast
import collections
import re

# --------------------------------------------------------------------------------
# Mass profile base class
# --------------------------------------------------------------------------------

class MassProfile(object):

    cfunction = None

    def __init__(self, width, mass):
        self.width = width
        self.mass = mass

    def create_fit_function_str(self, param_vals=None, param_prefix=""):
        raise NotImplementedError("MassProfile: Subclasses should overrode create_fitting_str")

    def create_constraint_str(self, param_prefix=""):
        """Returns a constraints string for the Fit algorithm

        :param param_prefix: An optional prefix for the parameter name
        """
        try:
            return "{0:f} < {1} < {2:f}".format(self.width[0], param_prefix + "Width", self.width[2])
        except TypeError:
            return ""

    def create_ties_str(self, param_prefix=""):
        """Return a ties string for the Fit algorithm

        :param param_prefix: An optional prefix for the parameter name
        """
        if not isinstance(self.width, collections.Iterable):
            return "{0}={1:f}".format(param_prefix + "Width", self.width)
        else:
            return ""

# --------------------------------------------------------------------------------
# Gaussian profile
# --------------------------------------------------------------------------------

class GaussianMassProfile(MassProfile):

    cfunction = "GaussianComptonProfile"

    @classmethod
    def from_str(cls, func_str, mass):
        """Attempt to create an object of this type from a string.
        Raises a TypeError if parsing fails
        """
        profile_prefix = "function=Gaussian,"
        if not func_str.startswith(profile_prefix):
            raise TypeError("Gaussian function string should start with 'function=Gaussian,'")

        func_str = func_str[len(profile_prefix):]
        # The only remaining property should be the width
        if func_str.startswith("width="):
            # Trim off width= prefix
            width = ast.literal_eval(func_str[6:])
        else:
            raise TypeError("Unexpected value in function string. Expected width= following function")

        return GaussianMassProfile(width, mass)

    def create_fit_function_str(self, param_vals=None, param_prefix=""):
        """Creates a string used by the Fit algorithm for this profile

        :param param_vals: A table of values for the parameters that override those set
        on the object already. Default=None
        :param param_prefix: A string prefix for the parameter as seen by the Mantid Fit algorithm
        """
        vals_provided = (param_vals is not None)

        if vals_provided:
            def_width = param_vals[param_prefix + "Width"]
        else:
            def_width = self.width
            if isinstance(def_width, list):
                def_width = def_width[1]

        fitting_str = "name={0},Mass={1:f},Width={2:f}".format(self.cfunction, self.mass, def_width)
        if vals_provided:
            param_name = "Intensity"
            intensity_str = "{0}={1:f}".format(param_name, param_vals[param_prefix + param_name])
            fitting_str += "," + intensity_str

        return fitting_str + ";"

# --------------------------------------------------------------------------------
# GramCharlier profile
# --------------------------------------------------------------------------------

class GramCharlierMassProfile(MassProfile):

    cfunction = "GramCharlierComptonProfile"

    def __init__(self, width, mass, hermite_coeffs, k_free, sears_flag):
        super(GramCharlierMassProfile, self).__init__(width, mass)

        self.hermite_co = hermite_coeffs
        self.k_free = k_free
        self.sears_flag = sears_flag

    @classmethod
    def from_str(cls, func_str, mass):
        """Attempt to create an object of this type from a string.
        Raises a TypeError if parsing fails
        """
        profile_prefix = "function=GramCharlier,"
        if not func_str.startswith(profile_prefix):
            raise TypeError("GramCharlier function string should start with 'function=GramCharlier,'")

        key_names = [("width", cls._parse_list),
                     ("hermite_coeffs", cls._parse_list),
                     ("k_free", cls._parse_bool_flag),
                     ("sears_flag", cls._parse_bool_flag)]
        # Possible key names:
        parsed_values = []
        for key, parser in key_names:
            try:
                parsed_values.append(parser(func_str, key))
            except ValueError, exc:
                raise TypeError(str(exc))

        return GramCharlierMassProfile(parsed_values[0], mass, parsed_values[1],
                                       parsed_values[2], parsed_values[3])

    def create_fit_function_str(self, param_vals=None, param_prefix=""):
        """Creates a string used by the Fit algorithm for this profile

        :param param_vals: A table of values for the parameters that override those set
        on the object already. Default=None
        :param param_prefix: A string prefix for the parameter as seen by the Mantid Fit algorithm
        """
        vals_provided = (param_vals is not None)

        if vals_provided:
            def_width = param_vals[param_prefix + "Width"]
        else:
            def_width = self.width
            if isinstance(def_width, list):
                def_width = def_width[1]

        def to_space_sep_str(collection):
            _str = ""
            for item in collection:
                _str += " " + str(item)
            return _str.lstrip()
        hermite_str = to_space_sep_str(self.hermite_co)

        fitting_str = "name={0},Mass={1:f},HermiteCoeffs={2},Width={3:f}".format(self.cfunction,
                                                                                 self.mass, hermite_str,
                                                                                 def_width)
        if vals_provided:
            par_names = ["FSECoeff"]
            for i,c in enumerate(self.hermite_co):
                if c > 0:
                    par_names.append("C_{0}".format(2*i))
            for par_name in par_names:
                fitting_str += ",{0}={1:f}".format(par_name, param_vals[param_prefix + par_name])

        return fitting_str + ";"

    def create_ties_str(self, param_prefix=""):
        """Return a ties string for the Fit algorithm

        :param param_prefix: An optional prefix for the parameter name
        """
        ties = super(GramCharlierMassProfile, self).create_ties_str(param_prefix)
        if not self.k_free:
            # Sears flag controls value of FSECoeff
            param_name = param_prefix + "FSECoeff"
            if self.sears_flag == 1:
                # tie to multiple of the width
                tied_value = param_prefix + "Width*sqrt(2)/12"
            else:
                tied_value = "0"
            if ties != "":
                ties += ","
            ties += "{0}={1}".format(param_name, tied_value)

        return ties

    @classmethod
    def _parse_list(cls, func_str, prop_name):
        """
        Parse a list from a string containing 'prop_name=[]'

        :param prop_name: The string on the lhs of the equality
        :return: The parsed list
        """
        prop_re = re.compile(prop_name + r"=(\[(?:\w(?:,)?(?:\s)?)+\])")
        match = prop_re.search(func_str)
        if match:
            value = ast.literal_eval(match.group(1))
            if not isinstance(value, list):
                raise ValueError("Unexpected format for {0} value. Expected e.g. {0}=[1,0,1]".format(prop_name))
        else:
            raise ValueError("Cannot find {0}= in function_str".format(prop_name))

        return value

    @classmethod
    def _parse_bool_flag(cls, func_str, prop_name):
        """
        Parse an integer from a string containing 'prop_name=1'

        :param prop_name: The string on the lhs of the equality
        :return: The parsed value
        """
        prop_re = re.compile(prop_name + r"=([0,1])")
        match = prop_re.search(func_str)
        if match:
            value = ast.literal_eval(match.group(1))
            if not isinstance(value, int):
                raise ValueError("Unexpected format for {0} value. Expected e.g. {0}=1".format(prop_name))
        else:
            raise ValueError("Cannot find {0}= in function_str".format(prop_name))

        return value


# --------------------------------------------------------------------------------
# Factory function
# --------------------------------------------------------------------------------

def create_from_str(func_str, mass):
    """Try and parse the function string to give the required profile function

        :param func_str: A string of the form 'function=Name,attr1=val1,attr2=val2'
        :param mass: The value of the mass for the profile
    """
    known_types = [GaussianMassProfile, GramCharlierMassProfile]

    errors = dict()
    for cls in known_types:
        try:
            return cls.from_str(func_str, mass)
        except TypeError, exc:
            errors[str(cls)] = str(exc)

    # if we get here we were unable to parse anything acceptable
    msgs = ["{0}: {1}".format(name, error) for name, error in errors.iteritems()]
    raise ValueError("\n".join(msgs))