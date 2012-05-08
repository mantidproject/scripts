from numpy import arange, degrees
from mantid import *
from mantid.simpleapi import *
from mantid.kernel import V3D



# if no angles set, it will mask from 0 to 180
# default TTmin = 0
# default TTmax = 180
def MaskAngle(**kwargs):


    workspace=kwargs.get('Workspace',None)
    if (workspace==None):
        raise RuntimeError("Workspace not set for angle mask")
    ttmin = (kwargs.get('twothetamin',0.0))
    ttmax = (kwargs.get('twothetamax',180.0))

#    ttmin = kwargs.get('twothetamin')
#    ttmax = kwargs.get('twohetamax')
    print ttmin, ttmax


    if ttmin== None:
        ttmin = 0.0
    if ttmax == None:
        ttmax = 180.0

    #check for silly angle ranges
    if ttmin < 0 :
        ttmin = 0
    if ttmax > 180 :
        ttmax =180

    if ttmin > ttmax :
        raise ValueError("ttmin > ttmax, please reset angle range for masking")

    detlist=[]

    #get the number of spectra
    ws = mtd[workspace]
    numspec = ws.getNumberHistograms()
    #detlist=[]
    for i in range(numspec):
        det=ws.getDetector(i)
        if not det.isMonitor():
            tt=degrees(det.getTwoTheta(V3D(0,0,0),V3D(0,0,1)))
            if tt>= ttmin and tt<= ttmax:
                detlist.append(det.getID())

    if len(detlist)> 0:
        MaskDetectors(Workspace=workspace,DetectorList=detlist)
    else:
        print "no detectors within this range"
    return detlist

		

