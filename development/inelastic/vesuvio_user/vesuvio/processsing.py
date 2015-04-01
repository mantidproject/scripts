"""
Defines functions and classes to start the processing of Vesuvio data. The main entry point that users should care
about is process_data().
"""
import mantid
from mantid.simpleapi import VesuvioTOFFit

def process_data(runs, flags):
    """
    The main entry point for user scripts. It simple takes the user input options
    and massages them into inputs the the VesuvioReduction algorithm can understand.

    :param runs: A string specifying the runs to process
    :param flags: A dictionary of flags to control the processing
    :return: None
    """
    # Transform inputs into something the algorithm can understand
    profile_descr = flags['masses']
    mass_values, profiles = [], []
    for mass_prop in profile_descr:
        function_props = []
        for key, value in mass_prop.iteritems():
            if key == 'value':
                mass_values.append(value)
            else:
                function_props.append("{0}={1}".format(key,value))
        profiles.append(",".join(function_props))
    profiles = ";".join(profiles)

    # Put the constraint arguments into a semi-colon separated string
    constraints = flags['intensity_constraints']
    if not isinstance(constraints[0], list):
        constraints = list(constraints,)
    constraints = [str(c) for c in constraints]
    intensity_constraints = ";".join(constraints)

    fitted, params = VesuvioTOFFit(Runs=runs, IPFilename=flags['ip_file'],
                                      Masses=mass_values,
                                      MassProfiles=profiles,
                                      IntensityConstraints=intensity_constraints,
                                      DifferenceMode=flags['diff_mode'])
    return fitted, params
