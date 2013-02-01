from qtiGenie import *
from PySlice2 import *
inst='mar'
iliad_setup(inst)
ext='.raw'
mapfile='mari_res'
#det_cal_file must be specified if the reduction sends out put to a workpsace
cal_file='MAR16637.raw'
#load vanadium file
whitebeamfile="16637"
LoadRaw(Filename=whitebeamfile,OutputWorkspace="wb_wksp",LoadLogFiles="0")
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
ei=100
rebin_params='-10,.2,95'
#load run
#runfile="16644"

runs=[16654]# ,16643,16652]
#16642,16643,16644,16645,16648,16649,16652]
#save .spe file
for runfile in runs:
	save_file=inst+str(runfile)+'_norm.spe'
	LoadRaw(Filename=runfile,OutputWorkspace="run_wksp",LoadLogFiles="0")
	#w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method ='current')
	w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current',diag_sigma=2,mask_run=16652,diag_remove_zero=True)
	SaveSPE('w1',save_file)
	
w2=data2D(w1)
w2.rebinproj('0,.07,10')
w2.display(10)
w2.CutAlongE(0,3,-10,1,80)
w2.CutAlongE(3,5,-10,1,80,Over=True)
w2.CutAlongQ(30,50,1,.1,10)
