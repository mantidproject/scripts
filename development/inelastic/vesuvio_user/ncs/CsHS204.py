from mantid.simpleapi import *
import ncs

runs = "15039-15045"
spectra = "135"
diff_type="SingleDifference" # Allowed values=Single,Double,Thick
ip_file = "IP0004_10.par"

## Define fitting options ##
fit_options = ncs.FitOptions()
fit_options.workspace_index = 0
fit_options.smooth_points = None
fit_options.bad_data_error = 1e6
fit_options.background_order = 3 # None to switch off

# Mass options
mass1 = {'value':1.0079, 'widths':[2,5,7], 'function':'GramCharlier',
          'hermite_coeffs':[1,0,0],'k_free':True, 'sears_flag':1}
mass2 = {'value':16.0, 'widths':10, 'function':'Gaussian'}
mass3 = {'value':27.0, 'widths':13, 'function':'Gaussian'}
mass4 = {'value':133.0, 'widths':30, 'function':'Gaussian'}

fit_options.masses = [mass1, mass2, mass3, mass4]

## Intensity constraints
fit_options.constraints = ([0,1,0,-4])


## Load data & preprocess ##

tof_data = LoadVesuvio(Filename=runs, SpectrumList=spectra,
                     Mode=diff_type,InstrumentParFile=ip_file)
tof_data = CropWorkspace(tof_data,XMin=50.0,XMax=562.0)
tof_data = ncs.preprocess(tof_data, fit_options)

## Run fitting ##
reduced_chi_square, params_ws = ncs.run_fit(tof_data, fit_options)

## Print results to screen ##
ncs.display_fit_output(reduced_chi_square, params_ws,fit_options)
