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
diff_type="SingleDifference" # Allowed values=Single,Double,Thick
ip_file = "IP0004_10.par"

print "Loading spectra %s from runs %s" % (spectra,runs)
raw_ws = LoadVesuvio(Filename=runs, SpectrumList=spectra,
                                     Mode=diff_type,InstrumentParFile=ip_file)

raw_ws = CropWorkspace(raw_ws,XMin=50.0,XMax=562.0)

#################################################################################
## Preprocessing options
#################################################################################
reducer_options = ReductionOptions()

## Smoothing (can be turned off by setting to None)
reducer_options.smooth_points = None

## Cropping (can be turned off by setting to None)
reducer_options.bad_data_error = 1e6

#################################################################################
## Fitting options
#################################################################################

## Background (Chebyshev polynomial)
reducer_options.background_order = None

#------------------- FIRST MASS ------------------------------------
## Hermite polynomial 
## List of 1/0 indicating if this term in the polynomial should be included. Length
## should equal the number of terms to include
reducer_options.hermite_coeffs = [1,0,0]

# If true then FSE coefficient remains free during the fitting, otherwise
# it is contrained to either 0 (sears_flag=0) or \sigma_H*sqrt(2)/12 (sears_flag=1) where \sigma_H is the
# standard deviation of the momentum distribution
reducer_options.k_free = False
reducer_options.sears_flag = 1
#-------------------------------------------------------------------

## Mass distributions
reducer_options.masses = [1.0079, 16.0, 27.0, 133.0]
## Widths - A list same length as masses list
##    Constrained Width: The entry should contain 3 numbers (low,value,high)
##    Fixed Width: The entry should be a single value, the width
reducer_options.widths = [(2,5,7), 10, 13, 30]

## Intensity constraints
##   Aeq*X = 0
## where Aeq is a matrix of size (nconstraints,nmasses)
reducer_options.constraints = ([0,1,0,-4])

#################################################################################
## Run
#################################################################################
reducer_options.workspace_index = 0
reducer_options.output_prefix = runs

raw_ws = preprocess(raw_ws, reducer_options)
reduced_chi_square, params_ws = run_fitting(raw_ws, reducer_options)

# Uncomment to display raw parameters
#for row in params_ws:
#    print row.values()

###############################################################################
# Process results
################################################################################

nmasses = len(reducer_options.masses)
nc = len(reducer_options.hermite_coeffs)

# Widths
wg_best = []
del_wg_best = []
# Peak areas
Av_best = []
del_Av_best = []
# Expansion coefficents for protons
c_best = []
del_c_best = []
# FSE term
k_best = 0.0
del_k_best = 0.0

m_i, c_i = 0,0
for row in params_ws:
    name,value,error = row.values()
    if "C_" in name: 
        term_index = (int(name[-1])/2)
        c_best.append(value)
        del_c_best.append(error)
    if "C_0" in name or "Intensity" in name:
        Av_best.append(value)
        del_Av_best.append(error)
    elif 'Width' in name:
        wg_best.append(value)
        del_wg_best.append(error)
    elif 'FSE' in name:
        k_best = value
        del_k_best = error

# Fill in zeroes for non-active hermite coeffs
n_active_c = len(c_best)
if nc > n_active_c:
    for i in range(n_active_c,nc):
        c_best.append(0.0)
        del_c_best.append(0.0)

print c_best
#######################################
# Peak areas are normalized to the sum
#######################################
Av_best = numpy.array(Av_best)
Av_best_sum = numpy.sum(Av_best)
Av_best_norm = Av_best/Av_best_sum
del_Av_best_norm = []

# Errors by full differential
jacobian = numpy.empty(shape=(nmasses,nmasses))
for j in range(nmasses):
    for k in range(nmasses):
        epsilon_k = del_Av_best[k]
        if k == j:
            df_j_k_r = (Av_best[k]+0.01*Av_best[k])/(Av_best_sum+0.01*Av_best[k])
            df_j_k_l = (Av_best[k]-0.01*Av_best[k])/(Av_best_sum-0.01*Av_best[k])
        else:
            df_j_k_r = Av_best[k]/(Av_best_sum+0.01*Av_best[k])
            df_j_k_l = Av_best[k]/(Av_best_sum-0.01*Av_best[k])
            
        jacobian[j][k] = (df_j_k_r - df_j_k_l)/epsilon_k

for j in range(nmasses):
    del_Av_best_norm_j = 0.0 
    for k in range(nmasses):
        df_j_k = jacobian[j][k]
        Av_best_k = Av_best[k]
        del_Av_best_k = del_Av_best[k]
        del_Av_best_norm_j += (df_j_k*del_Av_best_k)**2

    del_Av_best_norm.append(math.sqrt(del_Av_best_norm_j))

##########################################
# Expansion coefficients normalized to C_0
##########################################
# Insert zeroes for those that are off
# and compute reduced coefficients
a_best = []
del_a_best = []

c0 = c_best[0]
dc0 = del_c_best[0]
if c0 == 0.0: 
    c0 = 1

if reducer_options.k_free:
    del_k_best = math.fabs(del_k_best*c0 - dc0*k_best)/c0/c0
    k_best = k_best/c0
else:
    del_k_best = del_wg_best[0]*math.sqrt(2)/12

# Normalise expansion coefficients
for index, flag in enumerate(reducer_options.hermite_coeffs):
    c_i  = c_best[index]
    dc_i = del_c_best[index]
    if flag == 0:
        c_i = 0.0
        dc_i = 0.0
        c_best.insert(index, c_i)
        del_c_best.insert(index, dc_i)

    del_c_best[index] = math.fabs(dc_i*c0 - dc0*c_i)/c0/c0
    c_best[index] = c_best[index]/c0
    npoly = 2*index
    denom = (2**npoly)*math.factorial(index)
    a_best.append(c_best[index]/denom)
    del_a_best.append(del_c_best[index]/denom)

######################################
# Display
######################################

print
print 'Reduced Chi-Square =',reduced_chi_square
print    
print 'Fitting in the TOF space'

for i in range(nmasses):
    print'--------------------------------------------------------------------------------------------------------------------------------------'

    print 'The mass M(%d)=%f' % (i+1,reducer_options.masses[i])
    print
    print 'Parameters values in the Y space:'
    print

    # --- Currently on displayed in the results log as we can't yet pull them from Fit ---
    print 'See results log for resolution parameters w_L1, w_L0 etc'
    # disp(['w_L_1 (FWHM)= ',num2str(w_L1(:,i))])
    # disp(['w_L_0 (FWHM)= ',num2str(w_L0(:,i))])
    # disp(['w_T_h_e_t_a (FWHM)= ',num2str(w_theta(:,i))])
    # disp(['w_F_o_i_l_L_o_r_e_n_t_z (FWHM)= ',num2str(wl(:,i))])
    # disp(['w_F_o_i_l_G_A_U_S_S (FWHM) = ',num2str(wE_gauss(:,i))])
    
    print 'St. dev. of momentum distr. = %f +/- %f' % (wg_best[i],del_wg_best[i])
    print 'Scatt. int. (area, normalised) = %f +/- %f' % (Av_best_norm[i],del_Av_best_norm[i])
    
    if i == 0: # First mass
        for u in range(nc):
            print 'Hermite polynomial expansion coefficient c%d = %f +/- %f' % (2*u,c_best[u],del_c_best[u])
  
        for u in range(nc):
            print 'Hermite polynomial expansion coefficient a%d = c%d/(2^%d*%d!) = %f +/ %f' % (2*u,2*u,2*u,u,a_best[u],del_a_best[u])

        print
        print 'FSE coefficient k by the k/q He_3(y) expansion member = %f +/- %f' % (k_best, del_k_best)
        print
        print 'The coefficient k calculated in a harmonic oscillator model would be k = sigma*sqrt(2)/12 = %f' % (wg_best[i]*math.sqrt(2)/12)
    print


