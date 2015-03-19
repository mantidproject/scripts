import mantid

import math
import types

RAD2DEG = 180.0/math.pi

EVS_FWD_BANK_START = 135
EVS_SPEC_PER_BANK = 8
SUM_PLOT_SUFFIX = "_plotsum"

def plot(workspaces, **kwargs):
    """
        General plotting command.
          @param workspaces A workspace or list of workspaces references holding the data
                            (string names will work too)

        Accepted Keywords:
        ==================
        The following are exclusive options and only a single one is allowed to be specified at one time:
          spectra: An int or list of ints specifying the spectra numbers to plot
          angles:  A list or tuple of length 2 defining the min/max angles in degrees to plot (inclusive)
          bank:    A bank number between 1-8 (inclusive) defining the spectra to plot
        If multiple workspaces are given the ranges are applied separately to each.

        The following are general and apply to all plots types:
          errors:  If True then error bars are drawn on the plot
          sum:     If supplied and True the given input range on each workspace is summed and 
                   a single curve for each workspace is generated
          fig:     If supplied this window will be used to do the plot rather than creating a new one
          clrfig:  If True then clear the given window before plotting the user input
    """
    raise_error_if_not_in_gui("plot")

    # Work with a list of workspace references only
    workspaces = to_workspace_list(workspaces)

    # parse input
    range_kw, input_range = check_input(workspaces, **kwargs)
    do_sum = kwargs.get("sum", False)
    with_errors = kwargs.get("errors", False)
    plot_win = kwargs.get("fig", None)
    clrfig = kwargs.get("clrfig", True)

    # Get the function that will do the transform from user input to workspace indices
    range_transform_func, range_src = get_transform_func(range_kw, input_range)

    if do_sum:
        process_func = sum_range
    else:
        def pass_through(ws, rng, *args): return ws, rng
        process_func = pass_through

    plot_src = []
    for user_ws in workspaces:
        # Translate range inputs on each workspace to indices
        indices = range_transform_func(user_ws, range_src)
        plot_ws, indices = process_func(user_ws, indices, range_kw, input_range)
        plot_src.append((plot_ws, indices)) 

    # now plot
    from mantidplot import plotSpectrum
    for wksp, indices in plot_src:
        plot_win = plotSpectrum(wksp, indices, error_bars = with_errors,
                                window = plot_win, clearWindow = clrfig)
        if clrfig:
            # If clrfig was true then we only want it to clear for the first plot
            clrfig = False

    return plot_win

#----------------------------------------------------------------------------------------

def raise_error_if_not_in_gui(cmd):
    """
        Checks if we are running in MantidPlot and raises an
        error if not.
        @param cmd Name of the function that is trying to be called
    """
    if not mantid.__gui__:
        raise RuntimeError("Unable to use %s outside of MantidPlot" % cmd)

#----------------------------------------------------------------------------------------

def check_input(workspaces, **kwargs):
    """
        Raises an error if inputs are not a valid combination
    """
    # check for single spectrum workspaces that don't require a keyword
    single_spectrum = True
    for wksp in workspaces:
        if wksp.getNumberHistograms() != 1:
            single_spectrum = False

    if single_spectrum:
        return "single", 0

    range_inputs = ("spectra", "angles", "bank")
    count = 0
    for item in range_inputs:
        if item in kwargs:
            range_kw = item
            count += 1

    if count == 0:
        raise ValueError("Invalid input. Provide one of the following keywords: '%s'" % str(range_inputs))
    elif count > 1:
        raise ValueError("Invalid input. Only one of the following keywords can be provided: '%s'" % str(range_inputs))
    else:
        pass

    return range_kw, kwargs[range_kw]

#----------------------------------------------------------------------------------------

def to_workspace_list(src):
    """
        Returns a list of workspace references from the input
        @param src A string, a list of strings, a single workspace reference 
                   or a list of workspace references
    """
    src_type = type(src)
    if src_type == types.ListType:
        # ensure they are all workspace references
        workspaces = [to_workspace(item) for item in src]
    else:
        workspaces = [to_workspace(src)]

    return workspaces

#----------------------------------------------------------------------------------------

def to_workspace(src):
    """
        Converts a string to a workspace reference
    """
    if type(src) == types.StringType:
        return mantid.AnalysisDataService[src]
    else:
        # assume it is a workspace reference
        return src

#----------------------------------------------------------------------------------------

def get_transform_func(src_type, src_data):
    """
        Returns the function responsible for mapping
        the user input to workspace indices along with the
        correct data for feeding this function
        @param src_type Keyword indicating plotting type
        @param src_data User input that in the case of bank plotting is transformed
                        to spectrum numbers
    """
    # Check what type of plot we want
    if src_type == "single":
        def return_zero(workspace, spectra):
            return 0
        range_transform_func = return_zero
        range_src = src_data
    elif src_type == "spectra":
        range_transform_func = indices_from_spectrum
        range_src = src_data
    elif src_type == "angles":
        range_transform_func = indices_from_angles
        range_src = src_data
    elif src_type == "bank":
        range_transform_func = indices_from_spectrum
        bank = src_data
        range_src = get_spectra_from_bank(bank)
    else:
        # check input should mean we don't get here
        raise ValueError("Unknown plot type requested: '%s'" % (str(src_type)))
    # end if block

    return range_transform_func, range_src

#----------------------------------------------------------------------------------------

def indices_from_spectrum(workspace, spectra):
    """
        Returns the workspace indices corresponding to the given
        spectrum numbers on the workspace
        @param workspace A workspace reference 
        @param spectra An int or list of ints specifying spectrum numbers
    """
    if type(spectra) != types.ListType:
        spectra = [spectra]

    indices = []
    nhist = workspace.getNumberHistograms()
    for i in range(nhist):
        spec = workspace.getSpectrum(i)
        if spec.getSpectrumNo() in spectra:
            indices.append(i)
        if len(indices) == len(spectra):
            break
    # end for loop

    nfound = len(indices)
    if nfound == 0:
        raise ValueError("Spectrum numbers '%s' not found in workspace '%s'" % \
                         (str(spectra),str(workspace)))
    if nfound != len(spectra):
        raise ValueError("Unable to find spectra '%s' in workspace '%s'" % \
                         (str(spectra[nfound:]),str(workspace)))

    return indices

#----------------------------------------------------------------------------------------

def indices_from_angles(workspace, angles):
    """
        Returns the workspace indices corresponding to the given
        angular range on the workspace
        @param workspace A workspace reference 
        @param angles A two-element list or tuple containing the min/max 
                    angular range in degrees (inclusive)
    """
    if len(angles) != 2:
        raise ValueError("Invalid number of angles given. Require min/max range, "
                         "found: %s" % str(angles))

    min, max = angles[0], angles[1]
    indices = []
    nhist = workspace.getNumberHistograms()
    for i in range(nhist):
        try:
            det = workspace.getDetector(i)
        except RuntimeError:
            continue
        theta = workspace.detectorTwoTheta(det)*RAD2DEG
        if theta >= min and theta <= max:
            indices.append(i)
    # end for loop

    if len(indices) == 0:
        raise ValueError("Unable to find any detectors in angular range '%s' in "
                        "workspace '%s'" % (str(angles),str(workspace)))

    return indices

#----------------------------------------------------------------------------------------

def get_spectra_from_bank(bank_no):
    """
        Given a bank number, returns the list of spectra within that bank
        @param Bank number between 1 & 8 (inclusive)
    """
    if bank_no < 1 or bank_no > 8:
        raise ValueError("Invalid bank number '%d'. Valid range: 1-8" % (bank_no))

    start_spec = EVS_FWD_BANK_START + (bank_no - 1)*EVS_SPEC_PER_BANK
    return range(start_spec, start_spec + EVS_SPEC_PER_BANK)

#----------------------------------------------------------------------------------------

def sum_range(user_ws, indices, *args):
    """
        Sums the given spectra on the workspace to give a new workspace with a single 
        spectrum
        @param user_ws Multi-spectrum workspace containing the spectra to sum
        @param indices The list of indices that define the range to sum
    """
    output_name = generate_unique_workspace_name(str(user_ws) + SUM_PLOT_SUFFIX,*args)
    summed =  mantid.simpleapi.SumSpectra(InputWorkspace=user_ws,
                                          OutputWorkspace=output_name,
                                          ListOfWorkspaceIndices=indices,
                                          EnableLogging=False)
    # the final range now only has a single spectrum and it will always be zero
    return summed, 0

#----------------------------------------------------------------------------------------

def generate_unique_workspace_name(prefix, plt_type, input_range):
    """
        Given a prefix generate a unique workspace name by appending an
        incrementing number to the end of the prefix
        @param prefix A string prefix for the name
        @param plt_type User-selected plot type
        @param start start of range 
        @param end end of range
    """
    if hasattr(input_range, "__len__"):
        start, end = input_range[0], input_range[-1]
        name = "%s_%s_%s_%s" % (prefix, str(plt_type), str(start), str(end))
    else:
        start = input_range
        name = "%s_%s_%s" % (prefix, str(plt_type), str(start))
    return name