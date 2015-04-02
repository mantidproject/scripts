"""
Defines functions and classes to start the processing of Vesuvio data. The main entry point that most users should care
about is fit_tof().
"""
import mantid
from mantid.simpleapi import VesuvioTOFFit

# --------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------

def fit_tof(runs, flags):
    """
    The main entry point for user scripts fitting in TOF.

    :param runs: A string specifying the runs to process
    :param flags: A dictionary of flags to control the processing
    :return: Tuple of (fitted workspace, fitted_params)
    """
    # Transform inputs into something the algorithm can understand
    mass_values, profiles_strs = _create_profile_strs_and_mass_list(flags['masses'])
    background_str = _create_background_str(flags.get('background', None))
    intensity_constraints = _create_intensity_constraint_str(flags['intensity_constraints'])

    fitted, params = VesuvioTOFFit(Runs=runs, IPFilename=flags['ip_file'],
                                   Masses=mass_values,
                                   MassProfiles=profiles_strs,
                                   Background=background_str,
                                   IntensityConstraints=intensity_constraints,
                                   DifferenceMode=flags['diff_mode'])
    return fitted, params

# --------------------------------------------------------------------------------
# Private Functions
# --------------------------------------------------------------------------------

def _create_profile_strs_and_mass_list(profile_flags):
    """
    Create a string suitable for the algorithms out of the mass profile flags
    and a list of mass values
    :param profile_flags: A list of dict objects for the mass profile flags
    :return: A string to pass to the algorithm & a list of masses
    """
    mass_values, profiles = [], []
    for mass_prop in profile_flags:
        function_props = ["function={0}".format(mass_prop["function"])]
        del mass_prop["function"]
        for key, value in mass_prop.iteritems():
            if key == 'value':
                mass_values.append(value)
            else:
                function_props.append("{0}={1}".format(key,value))
        profiles.append(",".join(function_props))
    profiles = ";".join(profiles)

    return mass_values, profiles

def _create_background_str(background_flags):
    """
    Create a string suitable for the algorithms out of the background flags
    :param background_flags: A dict for the background (can be None)
    :return: A string to pass to the algorithm
    """
    if background_flags:
        background_props = ["function={0}".format(background_flags["function"])]
        del background_flags["function"]
        for key, value in background_flags.iteritems():
            background_props.append("{0}={1}".format(key,value))
        background_str = ",".join(background_props)
    else:
        background_str = ""

    return background_str

def _create_intensity_constraint_str(intensity_constraints):
    """
    Create a string suitable for the algorithms out of the intensity constraint flags
    :param inten_constr_flags: A list of lists for the constraints (can be None)
    :return: A string to pass to the algorithm
    """
    if intensity_constraints:
        if not isinstance(intensity_constraints[0], list):
            intensity_constraints = [intensity_constraints,]
        # Make each element a string and then join them together
        intensity_constraints = [str(c) for c in intensity_constraints]
        intensity_constraints_str = ";".join(intensity_constraints)
    else:
        intensity_constraints_str = ""

    return intensity_constraints_str