from qtiGenie import *
from PySlice2 import *
from autoEi import *
inst='mar'
iliad_setup(inst)
ext='.raw'
mapfile='mari_res2012'
#det_cal_file must be specified if the reduction sends out put to a workpsace
#cal_file='MAR17578.raw'
cal_file='MAR18040.raw'
#load vanadium file

run='18249'
WB='18040'
ei,rebin_params=autoEi(str(run))
w1=iliad(WB,run,ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current')

run='18313'
WB='18040'
ei,rebin_params=autoEi(str(run))
w1=iliad(WB,run,ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current',save_format='')

run='18322'
WB='18040'
ei,rebin_params=autoEi(str(run))
w1=iliad(WB,run,ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current',save_format='')


run='18069'
WB='18040'
ei,rebin_params=autoEi(str(run))
w1=iliad(WB,run,ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current',save_format='')