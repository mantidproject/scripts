##
## Basic Reduction script for the VESUVIO
## instrument.
##
## In active development to iron out the reduction
##
from mantid.simpleapi import *
from vesuvio_routines import mask_data
import numpy
import math
from screen import constrain

################################################################################
## Load data (need to correct for positions in IP
## file)
################################################################################

runs = "15039-15045"
spectra = "135"
diff_type="Single" # Allowed values=Single,Double,Thick
ip_file = "/home/dmn58364/Scripts/inelastic/vesuvio/selena/IP0004_10_no_header.par"

print "Loading spectra %s from runs %s" % (spectra,runs)
raw_ws = LoadVesuvio(RunNumbers=runs, SpectrumList=spectra,
                     DifferenceType=diff_type,InstrumentParFile=ip_file)

#################################################################################
## Smoothing
#################################################################################
smooth_in = False
smooth_pts = 3
if smooth_in:
    print "Smoothing data using %d neighbours" % smooth_pts
    raw_ws = SmoothData(InputWorkspace=raw_ws, NPoints=[smooth_pts])

#################################################################################
## Cropping
#################################################################################
crop_in = True
bad_data_error = 1e6
if crop_in:
    mask_data(raw_ws, bad_data_error)

#################################################################################
## Background (Chebyshev polynomial)
#################################################################################
include_background = False
back_poly_order = None
if include_background:
    back_poly_order = 3
    print "Including background using Chebyshev polynomial order=%d" % back_poly_order
else:
    back_poly_order = 0

#################################################################################
## Hermite polynomial (for the FIRST mass only)
#################################################################################
# List of 1/0 indicating if this term in the polynomial should be included. Length
# should equal the number of terms to include 
c_free = [1, 0, 1]
n_c = len(c_free)
n_cfree = numpy.sum(numpy.array(c_free))

# If true then FSE coefficient remains free during the fitting, otherwise
# it is contrained to either 0 (sears_flag=0) or \sigma_H*sqrt(2)/12 (sears_flag=1) where \sigma_H is the
# standard deviation of the momentum distribution
k_free = False
sears_flag = 1

# If FSE not free, then sears flags should be off
if k_free is False:
    sears_flag = 0
    n_kfree = 0
else:
    n_kfree = 1

print "Number of terms in Hermite expansion=%d." % n_c,
terms = ""
for index, value in enumerate(c_free):
    if value == 1:
	    terms += "H_%d " % (2*index)
print "Expansion will contain terms=%s" % terms
print "Sears flag=%d" % (sears_flag)

#################################################################################
## Mass distribution inputs
#################################################################################
masses = [1.00794, 16.0, 270.0, 133.0]
gauss_width = [5, 10, 13, 30] # The only free width is the first mass
gauss_lower = [2, 10, 13, 30]
gauss_upper = [7, 10, 13, 30]
npeaks = len(masses)

#################################################################################
## Intensity constraints
##   Aeq*X = 0
## where Aeq is a matrix of size (nconstraints,nhermite(active) + nkfree + nmasses)
#################################################################################
aeq = [ [0,-1,0,4] ] # list of lists. First=mass1,second=mass2 etc

#################################################################################
## Run the fit
#################################################################################

ws_index = 0
def to_space_sep_str(collection):
    _str = ""
    for item in collection:
        _str += " " + str(item)
    return _str.lstrip()
mass_str = to_space_sep_str(masses)
hermite_str = to_space_sep_str(c_free)

function_str = "name=NCSCountRate,WorkspaceIndex=%d,Masses=%s,HermiteCoeffs=%s" \
     % (ws_index, mass_str, hermite_str)
ties = ""
constraints = ""

## Width contraints
for index, (wg_low, wg_width, wg_hi) in enumerate(zip(gauss_lower,gauss_width,gauss_upper)):
    par_name = "Sigma_%d" % index
    if wg_low == wg_hi:
        ties += "%s=%f," % (par_name,wg_low)
    else:
        constraints += "%f < %s < %f," % (wg_low, par_name, wg_hi)

## Intensity constraints
for index, value in enumerate(c_free): # These all relate to the first mass
    lhs_par = "C_%d" % index
    if value == 1:
        ties += "%s*%s=0," % (lhs_par,"Intens_0")
    
## FSE constraint
if not k_free:
    ties += "%s=%s,"
    par_name = "FSECoeff"
    if sears_flag == 1:
        value = "Sigma_0*sqrt(2)/12"
    else:
        value = "0"
    ties = ties % (par_name, value)
    
## Stoichiometry
for aeq_i in aeq:
    for index, value in enumerate(aeq_i):
        if abs(value) != 0.0:
            rhs_par = "Intens_%d" % (index)
            ties += "%f*%s=0," % (value, rhs_par)
        
ties = ties.rstrip(",")
constraints = constraints.rstrip(",")

print "-"*15,"Ties","-"*15
print ties
print "-"*15,"Constraints","-"*15
print constraints
