from math import fabs, sqrt
# Work out (roughly) the tof range based on the incident energy
l1 = fabs(input.getInstrument().getSource().getPos().getZ())
l2 = 4.54 # Rough average
ltot = l1 + l2

ei = mtd['Live:TubeTOF'].getRun().getLogData('Ei').value if mtd.workspaceExists('Live:TubeTOF') else GetEi(InputWorkspace=input,FixEi='1').getPropertyValue('IncidentEnergy')

tof = ltot * sqrt(1.675e-27/2/1.602e-22/float(ei)) * 1e6

tofmin = tof-10000 if tof > 10000 else 0
binning = tofmin, 100, tof+10000

Rebin(InputWorkspace=input,OutputWorkspace=input,Params=binning)
GroupDetectors(InputWorkspace=input,OutputWorkspace=output,MapFile='/SNS/HYSA/shared/adara/128x1pixels.xml')
