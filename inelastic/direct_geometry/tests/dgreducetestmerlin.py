from qtiGenie import *
from PySlice2 import *
inst='mer'
iliad_setup(inst)
ext='.raw'
#mapfile='one2one_094'
mapfile='rings_113'

#det_cal_file must be specified if the reduction sends out put to a workpsace
cal_file='6399'
#load vanadium file
wb_run="6399"
LoadRaw(Filename=wb_run,OutputWorkspace="wb_wksp",LoadLogFiles="0")
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
ei=18
rebin_params='-10,.2,15'

#runs=[16654]# ,16643,16652]
runs=6398
LoadRaw(Filename=runs,OutputWorkspace="run",LoadLogFiles="0")
mono_van=16654
wb_mono=16637
samp_rmm=1
samp_mass=1
monovan_mapfile='mari_res'
#w1=iliad_abs(wb_run,runs,mono_van,wb_mono,samp_rmm,samp_mass,ei,rebin_params,mapfile,monovan_mapfile,use_sam_msk_on_monovan=True)
w1=iliad('wb_wksp','run',ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current',diag_remove_zero=False,detector_van_range=[20,300])
w2=data2D(w1)
w2.rebinproj('0.2,.05,5')
w2.display(10)
w2.CutAlongE(0,3,-10,1,80)