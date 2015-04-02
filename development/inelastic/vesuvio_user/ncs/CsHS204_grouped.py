from mantid.simpleapi import *
import ncs

forward_banks = ((135, 142), (143, 150), (151, 158), (159, 166),
                             (167, 174), (175, 182), (183, 190), (191, 198))


runs = "15039-15045"
spectra = []
for srange in forward_banks:
    spectra.append("{0}-{1}".format(*srange))
spectra = ";".join(spectra)
sum_spectra = True
diff_type="SingleDifference" # Allowed values=Single,Double,Thick
ip_file = "IP0004_10.par"

## Define fitting options ##
fit_options = ncs.FitOptions()
fit_options.smooth_points = None
fit_options.bad_data_error = 1e6
fit_options.background_order = None # None to switch off

# Mass options
mass1 = {'value':1.0079, 'widths':[2,5,7], 'function':'GramCharlier',
          'hermite_coeffs':[1,0,0],'k_free':False, 'sears_flag':1}
mass2 = {'value':16.0, 'widths':10, 'function':'Gaussian'}
mass3 = {'value':27.0, 'widths':13, 'function':'Gaussian'}
mass4 = {'value':133.0, 'widths':30, 'function':'Gaussian'}

fit_options.masses = [mass1, mass2, mass3, mass4]

## Intensity constraints
fit_options.constraints = ([0,1,0,-4])


## Load data & preprocess ##

tof_data = LoadVesuvio(Filename=runs, SpectrumList=spectra,
                     Mode=diff_type,InstrumentParFile=ip_file,SumSpectra=sum_spectra)
tof_data = CropWorkspace(tof_data,XMin=50.0,XMax=562.0)
tof_data = ncs.preprocess(tof_data, fit_options)

## Run fitting ##
fitted_ws, fitted_params = [], []
for idx in range(tof_data.getNumberHistograms()):
    fit_options.workspace_index = idx
    reduced_chi_square, params_ws = ncs.run_fit(tof_data, fit_options)
    fitted_ws.append("fit_%d" % (idx+1))
    RenameWorkspace("fit_Workspace", OutputWorkspace=fitted_ws[-1])
    fitted_params.append("params_%d" % (idx+1))
    RenameWorkspace("fit_Parameters", OutputWorkspace=fitted_params[-1])

# Group
fitted_data= GroupWorkspaces(fitted_ws)
fitted_pars = GroupWorkspaces(fitted_params)

## Print results to screen ##
#ncs.display_fit_output(reduced_chi_square, params_ws,fit_options)
