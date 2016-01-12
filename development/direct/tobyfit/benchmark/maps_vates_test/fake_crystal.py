from mantid.simpleapi import *
import tobyfit
reload(tobyfit)

import math
import os
import sys

THIS_DIR = os.path.dirname(__file__)
DATA_ETC = os.path.join(THIS_DIR, 'data_etc')
sys.path.insert(0, DATA_ETC)
    
# --------------------------------------------------------------------------------------------------------------------------

def run_comparision(name, ws1, ws2, rel_errs, keep_mtx=True):
    nbins = 41
    if not hasattr(rel_errs, '__len__'):
        rel_errs = [rel_errs] * 4
    for i in range(4):
        print "Binning dimension %d to %d bins" % (i, nbins)
        as_matrix = to_matrix_ws(ws1,i, nbins)
        ws1_matrix = RenameWorkspace(as_matrix, OutputWorkspace='ws1_matrix_' + name + '_' + str(i))
        del as_matrix
        as_matrix = to_matrix_ws(ws2,i, nbins)
        ws2_matrix = RenameWorkspace(as_matrix, OutputWorkspace='ws2_matrix_' + name  + '_' + str(i))
        status = CompareWorkspaces(ws1_matrix, ws2_matrix, Tolerance=rel_errs[i],ToleranceRelErr=True)
        if not keep_mtx:
            DeleteWorkspace(ws1_matrix)
            DeleteWorkspace(ws2_matrix)
        print 'Dimension %d signals equal (within rel. error of %.3f) = %r'% (i, rel_errs[i], status[0])

def to_matrix_ws(ws, non_integrated, nbins):
    kwargs = {'AxisAligned': True}
    for i in range(4):
        dim = ws.getDimension(i)
        dim_string = '%s,%.5f,%.5f' % (dim.getName(), dim.getMinimum(), dim.getMaximum())
        if i == non_integrated:
            nbins = 10
        else:
            nbins = 1
        dim_string += ',' + str(nbins)
        kwargs['AlignedDim' + str(i)] = dim_string
    #endfor
    __binned = BinMD(InputWorkspace=ws, **kwargs)
    __binned_query = QueryMDWorkspace(__binned)
    xname = __binned_query.getColumnNames()[3+non_integrated]
    yname = __binned_query.getColumnNames()[0]
    __matrix_ws = ConvertTableToMatrixWorkspace(__binned_query,ColumnX=xname,ColumnY=yname)
    DeleteWorkspace(__binned)
    DeleteWorkspace(__binned_query)
    return __matrix_ws


#============================================================
# Reduction
#============================================================
# Use same detector parmeters as tobyfit
parfile_path =  os.path.join(THIS_DIR, '4to1_102.par')
datfile_path =  os.path.join(THIS_DIR, '4to1_102.dat')
if not os.path.exists(datfile_path):
    tobyfit.translate_par_file(parfile_path, datfile_path)

nxs_file = 'map24076_ei50.nxs'
nxs_filepath = os.path.join(DATA_ETC, nxs_file)
if os.path.exists(nxs_filepath):
    if 'spe' in AnalysisDataService:
        spe = mtd['spe']
    else:
        spe = LoadNexusProcessed(Filename=nxs_filepath)
else:
    reduction_script = os.path.join(DATA_ETC, 'reduction.py')
    execfile(reduction_script)
    output_ws = 'converted_to_energy_transfer_ws'
    spe = RenameWorkspace(output_ws)

    # Preparations to match TobyFit
    #   - clear mask flags
    #   - use the same .par file positions
    #   - convert the data to point data so that we get the same MD boundaries (I suspect ConvertToMD should do this)
    #   - move the moderator instrument component back to it's IDF position for the the resolution calculation
    UpdateInstrumentFromFile(spe, Filename=datfile_path,AsciiHeader="spectrum,R,theta,phi")
    spe = ConvertToPointData(InputWorkspace=spe)
    base_moderator = spe.getInstrument().getBaseInstrument().getSource()
    mod_pos = base_moderator.getPos()
    MoveInstrumentComponent(spe,ComponentName=base_moderator.getName(),X=mod_pos.getX(),
                            Y=mod_pos.getY(),Z=mod_pos.getZ(),RelativePosition=False)
    SaveNexusProcessed(InputWorkspace=spe, Filename=nxs_file)
    DeleteWorkspace('SR_CurrentMasking')
    DeleteWorkspace('WB_CurrentMasking')
    DeleteWorkspace('WB_MAP024019_norm_white')
ClearMaskFlag(spe)

#============================================================
# Simulation
#============================================================
params = {}
params['ei'] = 50
params['uvec'] = [0,0,1]
params['vvec'] = [1,0,0]
params['psi'] = 0.0
params['omega'] = 0.0
params['resolution_params'] = 'MCLoopMin=10,MCLoopMax=10,MCType=1,ForegroundOnly=1'
params['foreground_model'] = 'QCoordinate'

def compare_with_tobyfit(spe, params, crystal_system):
    # === H ===
    print 'Comparing %s lattice, Q_H' % crystal_system
    params['foreground_params'] = 'Coord=H'
    fake_name = 'fake_' + crystal_system
    simul = tobyfit.run_simulation(fake_name + '_h_mt', spe, **params)
    # Load tobyfit comparison
    filename = fake_name + '_h.sqw'
    sqw = LoadSQW(os.path.join(THIS_DIR, filename), Q3DFrames='Q_lab', OutputWorkspace=filename)
    run_comparision('h',sqw, simul, rel_errs=1e-3)

    # === K ===
    print 'Comparing %s lattice, Q_H' % crystal_system
    params['foreground_params'] = 'Coord=K'
    fake_name = 'fake_' + crystal_system
    simul =  tobyfit.run_simulation(fake_name + '_k_mt', spe, **params)
    # Load tobyfit comparison
    filename = fake_name + '_k.sqw'
    sqw = LoadSQW(os.path.join(THIS_DIR, filename), Q3DFrames='Q_lab', OutputWorkspace=filename)
    run_comparision('k', sqw, simul, rel_errs=[1e-3,1e-3,1e-2,1e-3])

    # === L ===
    print 'Comparing %s lattice, Q_H' % crystal_system
    params['foreground_params'] = 'Coord=L'
    fake_name = 'fake_' + crystal_system
    simul = tobyfit.run_simulation(fake_name + '_l_mt', spe, **params)
    # Load tobyfit comparison
    filename = fake_name + '_l.sqw'
    sqw = LoadSQW(os.path.join(THIS_DIR, filename), Q3DFrames='Q_lab', OutputWorkspace=filename)
    run_comparision('l',sqw, simul, rel_errs=1e-3)
#end

# ================ CUBIC =================
params['alatt'] = 5.
params['blatt'] = 5.
params['clatt'] = 5.
params['alpha'] = 90.0
params['beta'] = 90.0
params['gamma'] = 90.0
compare_with_tobyfit(spe, params, 'cubic')

# ============== TETRAGONAL =================
params['alatt'] = 3.
params['blatt'] = 3.
params['clatt'] = 10.
params['alpha'] = 90.0
params['beta'] = 90.0
params['gamma'] = 90.0
compare_with_tobyfit(spe, params, 'tetragonal')

# ============== HEXAGONAL =================
params['alatt'] = 3.
params['blatt'] = 3.
params['clatt'] = 10.
params['alpha'] = 90.0
params['beta'] = 90.0
params['gamma'] = 120.0
compare_with_tobyfit(spe, params, 'hexagonal')