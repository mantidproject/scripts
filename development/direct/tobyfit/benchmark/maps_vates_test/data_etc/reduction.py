from mantid import config
from iliad_maps import *
import time
import os
import datetime

config['defaultsave.directory'] = '/home/dmn58364/Scripts/tobyfit/benchmark/maps_vates_test/data_etc'

#white beam vanadium run number
wbvan=24019

#mass of sample and RMM for absolute normalisation
sam_mass=0
sam_rmm=0


#----------------------------------------------------------------------------------
# Ei=80 250Hz
ei=50
rebin_pars=[-5,0.5,45]
bg_range=[15000,19000]
monovan=[]

# T=8K
runno=24076
iliad_maps_crystal(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range)
