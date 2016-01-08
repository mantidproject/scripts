import mantid

from mantid.kernel import funcreturns as _funcreturns
from mantid.simpleapi import *

RESOLUTION_MODEL = 'TobyFitResolutionModel'

def run_simulation(spe, **params):
    # Use the name of the variable that this function is assigned to as the output name
    output_workspace = _funcreturns.lhs_info()[1][0]

    # Required log entries, can be taken from real ones by placing an instrument parameter of the same
    # name pointing to the log name
    ei = params["ei"]
    AddSampleLog(Workspace=spe, LogName='Ei',LogText=str(ei), LogType="Number",NumberType='Double')

    # UB matrix
    alatt = params['alatt']
    blatt = params['blatt']
    clatt = params['clatt']
    uvec =  params['uvec']
    vvec = params['vvec']
    SetUB(Workspace=spe,a=alatt,b=blatt,c=clatt,u=uvec,v=vvec)

    # Sample rotation. Simulate 1 run at zero degrees psi
    psi = params['psi']
    omega = params['omega']
    alpha = params['alpha']
    beta = params['beta']
    gamma = params['gamma']
    AddSampleLog(Workspace=spe,LogName='psi',LogText=str(psi),LogType='Number',NumberType='Double')
    SetGoniometer(Workspace=spe,Axis0="psi,0,1,0,1")

    # Create the MD workspace in the lab frame as it cuts down on the number of transformations we have to make
    qframe = 'Q_lab'
    qscale = 'Q in A^-1'
    spe_q3d = ConvertToMD(InputWorkspace=spe, QDimensions="Q3D", Q3DFrames=qframe, QConversionScales=qscale,
                          OverwriteExisting=True,
                          SplitInto=[1],MaxRecursionDepth=1,TopLevelSplitting=False)

    # Do the simulation
    foreground_model = params['foreground_model']
    all_model_params =  params['resolution_params'] + ',' +  params['foreground_params']
    simulated = SimulateResolutionConvolvedModel(InputWorkspace=spe_q3d,
                                                             ResolutionFunction=RESOLUTION_MODEL,
                                                             ForegroundModel=foreground_model, 
                                                             Parameters=all_model_params,
                                                             OutputWorkspace=output_workspace)
    return simulated
# --------------------------------------------------------------------------------------------------------------------------

def translate_par_file(parfile_in, detector_dat):
    '''Translates a TobyFit .par file to a format Mantid can understand
    TobyFit columns: R,theta,-phi,-,-,spectrum
    First line is skipped
    '''
    parfile = open(parfile_in, 'r')
    datfile = open(detector_dat, 'w')
    for line in parfile:
        columns = line.split()
        if len(columns) != 6:
            continue
        spectrum = columns[5]
        r, theta, neg_phi = columns[:3]
        datfile.write("%s  %s    %s    %.4f\n" % (spectrum, r, theta, -float(neg_phi)))
    
# --------------------------------------------------------------------------------------------------------------------------
