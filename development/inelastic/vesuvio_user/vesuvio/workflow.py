"""
Defines functions and classes to start the processing of Vesuvio data. The main entry point that most users should care
about is fit_tof().
"""
from vesuvio.instrument import VESUVIO

from mantid import mtd
from mantid.simpleapi import (_create_algorithm_function, AlgorithmManager, CropWorkspace,
                              GroupWorkspaces, UnGroupWorkspace, LoadVesuvio, DeleteWorkspace)


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

    # Load
    spectra = flags['spectra']
    fit_mode = flags['fit_mode']
    tof_data = load_and_crop_data(runs, spectra, flags['ip_file'],
                                  flags['diff_mode'], fit_mode)

    # The simpleapi function won't have been created so do it by hand
    VesuvioTOFFit = _create_algorithm_function("VesuvioTOFFit", 1,
                                               AlgorithmManager.createUnmanaged("VesuvioTOFFit"))
    VesuvioCorrections = _create_algorithm_function("VesuvioCorrections", 1,
                                                    AlgorithmManager.createUnmanaged("VesuvioCorrections"))

    output_groups = []
    for index in range(tof_data.getNumberHistograms()):
        suffix = _create_fit_workspace_suffix(index, tof_data, fit_mode, spectra)

        # Corrections
        corrections_args = dict()

        # Need to do a fit first to obtain the parameter table
        pre_correction_pars_name = runs + "_params_pre_correction" + suffix
        corrections_fit_name = "__vesuvio_corrections_fit"
        VesuvioTOFFit(InputWorkspace=tof_data,
                      WorkspaceIndex=index,
                      Masses=mass_values,
                      MassProfiles=profiles_strs,
                      Background=background_str,
                      IntensityConstraints=intensity_constraints,
                      OutputWorkspace=corrections_fit_name,
                      FitParameters=pre_correction_pars_name)
        DeleteWorkspace(corrections_fit_name)
        corrections_args['FitParameters'] = pre_correction_pars_name

        if flags['ms_correct']:
            # Add the mutiple scattering arguments
            corrections_args.update(flags['ms_flags'])

        corrected_data_name = runs + "_tof_corrected" + suffix

        if flags['output_verbose_corrections']:
            corrections_args["CorrectionWorkspaces"] = runs + "_correction" + suffix
            corrections_args["CorrectedWorkspaces"] = runs + "_corrected" + suffix

        if flags['calculate_correction_proportion']:
            corrections_args["LinearFitResult"] = runs + "_correction_fit" + suffix

        VesuvioCorrections(InputWorkspace=tof_data,
                           OutputWorkspace=corrected_data_name,
                           WorkspaceIndex=index,
                           GammaBackground=flags['gamma_correct'],
                           Masses=mass_values,
                           MassProfiles=profiles_strs,
                           IntensityConstraints=intensity_constraints,
                           MultipleScattering=flags['ms_correct'],
                           **corrections_args)

        # Fit
        ws_name = runs + "_data" + suffix
        pars_name = runs + "_params" + suffix
        VesuvioTOFFit(InputWorkspace=corrected_data_name,
                      WorkspaceIndex=0, # Corrected data always has a single histogram
                      Masses=mass_values,
                      MassProfiles=profiles_strs,
                      Background=background_str,
                      IntensityConstraints=intensity_constraints,
                      OutputWorkspace=ws_name,
                      FitParameters=pars_name)
        DeleteWorkspace(corrected_data_name)

        group_name = runs + suffix
        output_workspaces = [ws_name, pars_name, pre_correction_pars_name]

        if flags['output_verbose_corrections']:
            output_workspaces += mtd[corrections_args["CorrectionWorkspaces"]].getNames()
            output_workspaces += mtd[corrections_args["CorrectedWorkspaces"]].getNames()
            UnGroupWorkspace(corrections_args["CorrectionWorkspaces"])
            UnGroupWorkspace(corrections_args["CorrectedWorkspaces"])

        if flags['calculate_correction_proportion']:
            output_workspaces.append(runs + "_correction_fit" + suffix)

        output_groups.append(GroupWorkspaces(InputWorkspaces=output_workspaces, OutputWorkspace=group_name))

    if len(output_groups) > 1:
        return output_groups
    else:
        return output_groups[0]


def load_and_crop_data(runs, spectra, ip_file, diff_mode='single',
                       fit_mode='spectra'):
    """
    @param runs The string giving the runs to load
    @param spectra A list of spectra to load
    @param ip_file A string denoting the IP file
    @param diff_mode Either 'double' or 'single'
    @param fit_mode If bank then the loading is changed to summing each bank to a separate spectrum
    """
    instrument = VESUVIO()
    load_banks = (fit_mode == 'bank')
    output_name = _create_tof_workspace_suffix(runs, spectra)

    if load_banks:
        sum_spectra = True
        if spectra == "forward":
            bank_ranges = instrument.forward_banks
        elif spectra == "backward":
            bank_ranges = instrument.backward_banks
        else:
            raise ValueError("Fitting by bank requires selecting either 'forward' or 'backward' "
                             "for the spectra to load")
        bank_ranges = ["{0}-{1}".format(x, y) for x, y in bank_ranges]
        spectra = ";".join(bank_ranges)
    else:
        sum_spectra = False
        if spectra == "forward":
            spectra = "{0}-{1}".format(*instrument.forward_spectra)
        elif spectra == "backward":
            spectra = "{0}-{1}".format(*instrument.backward_spectra)

    if diff_mode == "double":
        diff_mode = "DoubleDifference"
    else:
        diff_mode = "SingleDifference"

    kwargs = {"Filename": runs,
              "Mode": diff_mode, "InstrumentParFile": ip_file,
              "SpectrumList": spectra, "SumSpectra": sum_spectra,
              "OutputWorkspace": output_name}
    full_range = LoadVesuvio(**kwargs)
    return CropWorkspace(InputWorkspace=full_range, XMin=instrument.tof_range[0],
                         XMax=instrument.tof_range[1], OutputWorkspace=output_name)

# --------------------------------------------------------------------------------
# Private Functions
# --------------------------------------------------------------------------------

def _create_tof_workspace_suffix(runs, spectra):
    return runs + "_" + spectra + "_tof"

def _create_fit_workspace_suffix(index, tof_data, fit_mode, spectra):
    if fit_mode == "bank":
        return "_" + spectra + "_bank_" + str(index+1)
    else:
        spectrum = tof_data.getSpectrum(index)
        return "_spectrum_" + str(spectrum.getSpectrumNo())

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
