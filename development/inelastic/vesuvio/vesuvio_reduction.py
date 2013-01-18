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

raw_ws /= 1e6
#raw_ws = CropWorkspace(raw_ws,XMin=50)
tof = raw_ws.readX(0)
raw_ws = Rebin(raw_ws,Params=[50.0,1,tof[-1]])

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
if k_free:
    sears_flag = 0
    n_kfree = 1
else:
    n_kfree = 0
    
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
masses = [1.00794, 16.0, 27.0, 133.0]
gauss_width = [5, 10, 13, 30] # The only free width is the first mass
gauss_lower = [2, 10, 13, 30]
gauss_upper = [7, 10, 13, 30]
nmasses = len(masses)

#################################################################################
## Intensity constraints
##   Aeq*X = 0
## where Aeq is a matrix of size (nconstraints,nhermite(active) + nkfree + nmasses)
#################################################################################
aeq = [ [0,-1,0,4] ] # list of lists. First list=mass1,second list=mass2 etc

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

function_str = ""
ties = ""
constraints = ""

for index, mass in enumerate(masses):
    # Main function
    function_str += "name=%s,%s;"
    if index == 0:
        func_name = "GramCharlierComptonProfile"
        params = "WorkspaceIndex=%d,Mass=%f,HermiteCoeffs=%s" % (ws_index,mass,hermite_str)
    else:
        func_name = "GaussianComptonProfile"
        params = "WorkspaceIndex=%d,Mass=%f" % (ws_index,mass)
    function_str = function_str % (func_name, params)
    
    # Constraints/Ties
    func_index = index
    par_name = "f%d.Width" % func_index
    wg_low,wg_width,wg_hi = gauss_lower[index],gauss_width[index],gauss_upper[index]
    if wg_low == wg_hi:
        ties += "%s=%f," % (par_name,wg_low)
    else:
        constraints += "%f < %s < %f," % (wg_low, par_name, wg_hi)
        
### Stoichiometry/Intensity constraints
for aeq_i in aeq:
    # The first value in the aeq list relates to the first mass. We need to expand the
    # matrix to incorporate the hermite polynomial coeffs that are used instead of the Gaussian
    intensity_pars = ["f%d.Intensity" % index for index in range(1,nmasses)]
    firstmass_aeq = aeq_i[0]
    aeq_i = aeq_i[1:]
    for back_index, coeff in enumerate(reversed(c_free)):
        if coeff == 1:
            aeq_i.insert(0, firstmass_aeq) #prepend
            intensity_pars.insert(0,"f0.C_%d" % (2*(n_c - 1 - back_index)) ) #Even coefficents

    # Make list of tuples (coeff,par_name) that contribute
    lhs_terms = []
    for coeff, par_name in zip(aeq_i, intensity_pars):
        if abs(coeff) > 0.0:
            lhs_terms.append((coeff,par_name))
    
    # Pick first par name to go on lhs of tie expression.
    # and rearrange so everything else is on the the RHS of aeq_i*x_i = 0
    lhs = lhs_terms[0]
    rhs = lhs_terms[1:]
    denominator = -1.0*lhs[0]
    rhs = map(lambda term: (term[0]/denominator,term[1]), rhs)
    # Now form p=a*b+c*d...
    tie_i = lhs[1] + "="
    rhs_str = ""
    for coeff, param in rhs:
        rhs_str += "%+f*%s" % (coeff,param) # +f prints +/- with number
    ties += lhs[1] + "=" + rhs_str.lstrip("+") + ","

## FSE constraint
if not k_free:
    ties += "%s=%s"
    par_name = "f0.FSECoeff"
    if sears_flag == 1:
        value = "f0.Width*sqrt(2)/12"
    else:
        value = "0"
    ties = ties % (par_name, value)

function_str = function_str.rstrip(";")
ties = ties.rstrip(",")
constraints = constraints.rstrip(",")

print function_str
print "-"*15,"Ties","-"*15
print ties
print "-"*15,"Constraints","-"*15
print constraints

Fit(function_str,"raw_ws",Ties=ties,Constraints=constraints,Output="fit",CreateOutput=True,OutputCompositeMembers=True,MaxIterations=500)
