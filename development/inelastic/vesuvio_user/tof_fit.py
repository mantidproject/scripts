"""
Main driver script for VESUVIO users to do TOF fitting

It requires the vesuvio script
"""
import vesuvio.workflow
reload(vesuvio.workflow)
from vesuvio.workflow import fit_tof

# --------------------------------------------------------------------------------
# Standard flags to modify processing
# --------------------------------------------------------------------------------

# Specify the run(s) to process. Can be either a single number or a list of
# two numbers defining ranges
# Example: "15039" will process the single run
# Example 2: "15039-15045" will process the sum of range 15039-15045
# Example 3: "15039,15045" will process the sum of 15039,15045
runs = "15039-15045"

# Holds flags to alter how processing occurs
flags = dict()

# Fitting mode. Options are:
#    bank: spectra from each bank are summed and each bank is then fitted separately
#    spectrum: spectra are fitted individually
flags['fit_mode'] = 'spectra'

# Spectra selection. Can be a single number, a list of two numbers defining a range
# or one of the keyword strings "forward", "backward", "all"
# Example: 135 will only process spectrum 135 (only applies if fit_mode="spectra")
# Example 2: 135-185 will process all spectra in the range (only applies if fit_mode="spectra")
# Example 3: "forward" will process all spectra in the forward scattering banks
flags['spectra'] = '143-150'

# Binning options passed to the Rebin algorithm run on raw data after loading
# If set to None no rebinning is performed
flags['bin_parameters'] = None
# flags['bin_parameters'] = '50,2,562'

# Masses and properties for fit. Each mass is a dictionary of properties about that mass
# All functions require the following keywords:
#   'value': The mass value itself
#   'function': The name of the function approximate the mass profile: GramCharlier/Gaussian
#   'width': The value/values for the width. A single number specifies a fixed width whereas
#            a list of numbers specifies a default and allowed range for the width as [min,default,max]
# If the function is set to GramCharlier (generally only for the first mass) then there are additional
# keywords:
#    'hermite_coeffs': A list of 1/0 to indicate whether a given coefficient is included
#    'k_free': True/False to indicate whether k is either fixed or remains free
#    'sears_flag': If k is fixed then sears_flag=1 fixes k=sqrt(2)/12 whereas sears_flag=0
#                 fixes k=0
mass1 = {'value': 1.0079, 'function': 'GramCharlier', 'width': [2, 5, 7],
         'hermite_coeffs': [1,0,0], 'k_free': 0, 'sears_flag': 1}
mass2 = {'value': 16.0, 'function': 'Gaussian', 'width': 10}
mass3 = {'value': 27.0, 'function': 'Gaussian', 'width': 13}
mass4 = {'value': 133.0, 'function': 'Gaussian', 'width': 30}
flags['masses'] = [mass1, mass2, mass3, mass4]

# Intensity constraints. Can be None or a list of lists defining the required
# constraints to be imposed between the intensity values for each mass.
# Example 1: list([0, 1, 0, -4]) defines a single constraint that the intensity of mass 2
#            should be 4 times the intensity of mass 4
flags['intensity_constraints'] = list([0, 1, 0, -4])

# Background. Defines a function an associated attributes
# Currently only a Polynomial is supported. To switch off the
# background set
#flags['background'] = None
flags['background'] = {'function': 'Polynomial', 'order': 2}

# --------------------------------------------------------------------------------
# Corrections flags
# --------------------------------------------------------------------------------

# Run(s) for an empty container subtracted before final fitting is performed
# Set to None for no subtraction or a set of run numbers in the same format as runs variable
flags['container_runs'] = None

# Fixes the scale factor for the container, if not specified then the factor is
# calculated form the result of a linear fit
flags['fixed_container_scaling'] = None

# Outputs workspaces containing the correction factors for each correction and
# workspaces containing the input data with the single correction applied.
flags['output_verbose_corrections'] = True

# Enable gamma correction for gamma emissions due to neutron absorption in
# sheilding
flags['gamma_correct'] = True

# Scale factor for gamma correction
# If set to None the scale factor is calculated from a linear fit of
# gamma and total scattering to the input data
flags['fixed_gamma_scaling'] = None

# Holds flags specific to multiple scattering
flags['ms_flags'] = dict()

# Sample size in cm
flags['ms_flags']['SampleWidth'] = 10.0
flags['ms_flags']['SampleHeight'] = 10.0
flags['ms_flags']['SampleDepth'] = 0.5

# Sample density in g/cm^3
flags['ms_flags']['SampleDensity'] = 241

# Optional parameters (default values are given)
# flags['ms_flags']['Seed'] = 123456789
# flags['ms_flags']['NumScatters'] = 3
# flags['ms_flags']['NumRuns'] = 10
# flags['ms_flags']['NumEvents'] = 50000
# flags['ms_flags']['SmoothNeighbours'] = 3
# flags['ms_flags']['BeamRadius'] = 2.5


# --------------------------------------------------------------------------------
# Advanced flags
# --------------------------------------------------------------------------------

# Calibration file specifying the detector positions and parameter values
flags['ip_file'] = 'IP0004_10.par'

# Differencing mode. You should rarely need to modify this. Options are:
#    single
#    double
flags['diff_mode'] = 'single'

# Maximum number of iterations to use in fitting
flags['max_fit_iterations'] = 5000

# Configuration of the minimizer to use in fitting
flags['fit_minimizer'] = 'Levenberg-Marquardt,AbsError=1e-08,RelError=1e-08'

# Number of iterations to perform of the corrections and fitting workflow
# Each new iteration used the previous fitted parameters as the starting parameters
flags['iterations'] = 1

# The maxmimum change in the cost function value between iterations at which the
# parameters are assumed to have converged.
# When this condition is met iterations will be stopped, if the value is set to
# None then the number of iterations defined above will always be performed.
flags['convergence'] = None

# --------------------------------------------------------------------------------
# Run fit
# --------------------------------------------------------------------------------
fit_tof(runs, flags, flags['iterations'], flags['convergence'])
