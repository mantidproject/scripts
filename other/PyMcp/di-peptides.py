from qtiGenie imort *
iliad_setup('mar')
ala_pro_dry=iliad(12439,12520,250,'-20,.4,240','mari_res')
gly_ala_dry=iliad(12439,12503,250,'-20,.4,240','mari_res')
gly_pro_dry=iliad(12439,12478,250,'-20,.4,240','mari_res')
ala_pro_wet=iliad(12439,12555,250,'-20,.4,240','mari_res')
gly_ala_wet=iliad(12439,12537,250,'-20,.4,240','mari_res')
gly_pro_wet=iliad(12439,12571,250,'-20,.4,240','mari_res')


SofQW(InputWorkspace="ala_pro_dry",OutputWorkspace="bb",QAxisBinning="-1,0.1,15",EMode="Direct",EFixed="250")
Transpose(InputWorkspace="bb",OutputWorkspace="bbb")
Rebin(InputWorkspace="bbb",OutputWorkspace="bbb",Params="-10,1,250",PreserveEvents="0")
Transpose(InputWorkspace="bbb",OutputWorkspace="ala_pro_cut")
