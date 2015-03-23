"""Holds classes that define the mass profiles
"""
import ast
import re

# --------------------------------------------------------------------------------
# Mass profile interface
# --------------------------------------------------------------------------------

class MassProfile(object):

    cfunction = None

    def __init__(self, width, mass):
        self.width = width
        self.mass = mass

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

        func_str = func_str.lstrip(profile_prefix)
        # The only remaining property should be the width
        if func_str.startswith("width="):
            # Trim off width= prefix
            width = ast.literal_eval(func_str[6:])
        else:
            raise TypeError("Unexpected value in function string. Expected width= following function")

        return GaussianMassProfile(width, mass)

    def create_fitting_str(self, param_vals=None, param_prefix=""):
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

        fitting_str = "name={0},Mass={1},Width={2}".format(self.cfunction, self.mass, def_width)
        if vals_provided:
            param_name = "Intensity"
            intensity_str = "{0}={1}".format(param_name, param_vals[param_prefix + param_name])
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

        func_str = func_str.lstrip(profile_prefix)
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

    def create_fitting_str(self, param_vals=None, param_prefix=""):
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

        fitting_str = "name={0},Mass={1},HermiteCoeffs={2},Width={3}".format(self.cfunction,
                                                                             self.mass, hermite_str,
                                                                             def_width)
        if vals_provided:
            par_names = ["FSECoeff"]
            for i,c in enumerate(self.hermite_co):
                if c > 0:
                    par_names.append("C_%d" % (2*i))
            for par_name in par_names:
                fitting_str += ",%s=%f" % (par_name, param_vals[param_prefix + par_name])

        return fitting_str + ";"


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

def create_profile_from_str(func_str, mass):
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