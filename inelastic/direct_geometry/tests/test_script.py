from qtiGenie import *
from PySlice2 import *
#execfile('PySlice2.py')

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

runs=[16654]#,16643,16652,16642,16643,16644,16645,16648,16649]
#save .spe file
for runfile in runs:
	save_file=inst+str(runfile)+'_norm.spe'
	LoadRaw(Filename=str(runfile),OutputWorkspace="run_wksp",LoadLogFiles="0")
	
	#w2=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,det_cal_file=cal_file,diag_remove_zero=False,norm_method='current')
	w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current')
	#w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current',diag_sigma=2,diag_remove_zero=False,hardmaskPlus='/Users/jon/Work-computing/mantid_test_data/mari/mask.dat')
	#SaveSPE('w2',save_file)
	#DeleteWorkspace('_wksp.spe-white')
	#RenameWorkspace(w2,str(runfile))
print type(w1)
w2=data2D('w1')
w2.rebinproj('0.4,.06,15')
w2.display(10)
w2.ECut(0,1,10,.1,80)
w2.ECut(1,5,10,.1,80)
w2.ECut(6,10,10,1,80,over=True)
w2.ECut(1,12,0,.51,80)
w2.ECut(3,5,-10,1,80,Over=True)
w2.QCut(14.9,15.3,1,.02,15)
#w2.CutAlongQ(10,10.2,1,.08,10,Over=True)
w2.figures()
w2.mergeFigs(5,6)

w2.gofe(300,0.002)
w2.boseFac(250)
w2.smooth(2,10)
