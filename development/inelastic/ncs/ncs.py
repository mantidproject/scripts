##
## In active development to iron out the reduction
##
from mantid.simpleapi import *
import sys
if "ncs_helpers" in sys.modules:
    reload(sys.modules["ncs_helpers"])
from ncs_helpers import ReductionOptions, mask_data, preprocess, run_fitting
import numpy
import math

################################################################################
## Load data (need to correct for positions in IP
## file)
################################################################################

runs = "15039-15045"
spectra = "135"
diff_type="Single" # Allowed values=Single,Double,Thick
ip_file = "IP0004_10_no_header.par"

print "Loading spectra %s from runs %s" % (spectra,runs)
raw_ws = LoadVesuvio(RunNumbers=runs, SpectrumList=spectra,
                     DifferenceType=diff_type,InstrumentParFile=ip_file)

raw_ws = CropWorkspace(raw_ws,XMin=50)
tof = raw_ws.readX(0)
raw_ws = Rebin(raw_ws,Params=[50.0,1,tof[-1]])

#################################################################################
## Preprocessing options
#################################################################################
reducer_options = ReductionOptions()

## Smoothing
reducer_options.smooth_points = 3

## Cropping
reducer_options.bad_data_error = None

#################################################################################
## Fitting options
#################################################################################

## Background (Chebyshev polynomial)
reducer_options.background_order = None

#------------------- FIRST MASS ------------------------------------
## Hermite polynomial 
## List of 1/0 indicating if this term in the polynomial should be included. Length
## should equal the number of terms to include
reducer_options.hermite_coeffs = [1,0,1]

# If true then FSE coefficient remains free during the fitting, otherwise
# it is contrained to either 0 (sears_flag=0) or \sigma_H*sqrt(2)/12 (sears_flag=1) where \sigma_H is the
# standard deviation of the momentum distribution
reducer_options.k_free = False
reducer_options.sears_flag = 1
#-------------------------------------------------------------------

## Mass distributions
reducer_options.masses = [1.00794, 16.0, 27.0, 133.0]
## Widths - A list same length as masses list
##    Constrained Width: The entry should contain 3 numbers (low,value,high)
##    Fixed Width: The entry should be a single value, the width
reducer_options.widths = [(2,5,7), 10, 13, 30]

## Intensity constraints
##   Aeq*X = 0
## where Aeq is a matrix of size (nconstraints,nhermite(active) + nkfree + nmasses)
reducer_options.constraints = ([0,-1,0,4])

#################################################################################
## Run
#################################################################################
reducer_options.workspace_index = 0
reducer_options.output_prefix = runs

raw_ws = preprocess(raw_ws, reducer_options)
run_fitting(raw_ws, reducer_options)
