from qtiGenie import *
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
ei=250
rebin_params='-50,.5,240'
#load run
#runfile="16644"

runs=[16653]
#[16639,16640,16641,16646,16647,16650,
#save .spe file
for runfile in runs:
	save_file=inst+str(runfile)+'_norm.spe'
	LoadRaw(Filename=runfile,OutputWorkspace="run_wksp",LoadLogFiles="0")
	#w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method ='current')
	w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current',bkgd_range=[10000,19000],diag_sigma=2,diag_remove_zero=True)#,mask_run=16639)
	SaveSPE('w1',save_file)
