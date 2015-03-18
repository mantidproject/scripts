"""
Main driver script for VESUVIO users

It requires the VesuvioReduction algorithm
"""
import mantid
from mantid.simpleapi import VesuvioReduction

# --------------------------------------------------------------------------------
# Standard flags to modify processing
# --------------------------------------------------------------------------------

# Specify the run(s) to process. Can be either a single number or a list of
# two numbers defining ranges
# Example: 14188 will process the single run
# Example 2: [14188, 14195] will process the range 14188-14195
runs = [15039, 15045]

# Fitting mode. Options are:
#    bank: spectra from each bank are summed and each bank is then fitted separately
#    spectrum: spectra are fitted individually
fit_mode = 'bank'

# Spectra selection. Can be a single number, a list of two numbers defining a range
# or one of the keyword strings "forward" or "backward"
# Example: 135 will only process spectrum 135 (only applies if fit_mode="spectrum")
# Example 2: 135-185 will process all spectra in the range (only applies if fit_mode="spectrum")
# Example 3: "forward" will process all spectra in the forward scattering banks
spectra = 'forward'

# Masses. Each mass is specified by a Python dictionary and then combined into a
# list of masses. The dictionary should have the following keys:
#   'value': the actual mass in amu
#   'widths': a single number for a fixed with or a list of 3 numbers,
#             [min,start,max] for a fitted width
#   'function': Name of the fit function. Options: 'Gaussian', 'GramCharlier'
#
# 'GramCharlier' should generally be used for the first mass only. There are extra
# keys required if it used:
#   'hermite_coeffs': A list of 1/0 to indicate whether a given coefficient is included
#   'k_free': True/False to indicate whether k is either fixed or remains free
#   'sears_flag': If k is fixed then sears_flag=1 fixes k=sqrt(2)/12 whereas sears_flag=0
#                 fixes k=0
mass1 = {'value': 1.0079, 'widths':[2, 5, 7], 'function': 'GramCharlier',
         'hermite_coeffs': [1, 0, 0], 'k_free': True, 'sears_flag': 1}
mass2 = {'value': 16.0, 'widths': 10, 'function': 'Gaussian'}
mass3 = {'value': 27.0, 'widths': 13, 'function': 'Gaussian'}
mass4 = {'value': 133.0, 'widths': 30, 'function': 'Gaussian'}
masses = [mass1, mass2, mass3, mass4]

# Intensity constraints. Can be None or a tuple of lists defining the required
# constraints to be imposed between the intensity values for each mass.
# Example 1: ([0, 1, 0, -4]) defines a single constraint that the intensity of mass 2
#            should be 4 times the intensity of mass 4
constraints = ([0, 1, 0, -4])

# --------------------------------------------------------------------------------
# Advanced flags
# --------------------------------------------------------------------------------

# Calibration file specifying the detector positions and parameter values
ip_file = 'IP0004_10.par'

# Spectra summing. Only applicable if the fit_mode='spectra'. If True then all input
# spectra are summed together and processed
sum_spectra = True

# Differencing mode. You should rarely need to modify this. Options are:
#    double
#    single
diff_mode = 'double'

# --------------------------------------------------------------------------------
# Processing
# --------------------------------------------------------------------------------
#VesuvioReduction(Runs=runs
# , Spectra=spectra)
mantid.api.AlgorithmManager.create("VesuvioReduction").initialize()