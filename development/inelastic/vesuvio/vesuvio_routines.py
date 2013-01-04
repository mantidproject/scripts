##
## Functions to aid VESUVIO reduction
##
import numpy as np

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
        workspace = MaskBins(InputWorkspace=workspace, SpectraList=[int(ws_index)],
                             XMin=xmin,XMax=xmax)
