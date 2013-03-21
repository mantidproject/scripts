from math import sqrt, cos, radians, floor, ceil, fabs

# All this is to work out the Q range to pass to SofQW
Ei = input.getRun().getProperty('Ei').value
S2 = fabs(input.getRun().getProperty('s2').value[0])

xvals = input.readX(0)
dE = [ xvals[0], 0.0, xvals[xvals.size-1] ]

factor = 2.072194

Qlo = []
Qhi = []
for E in dE:
	kf2 = (Ei-E)/factor
	ki2 = Ei/factor
	Qhi.append( sqrt(ki2+kf2-(2*sqrt(ki2)*sqrt(kf2)*cos(radians(S2+30)))) )
	Qlo.append( sqrt(ki2+kf2-(2*sqrt(ki2)*sqrt(kf2)*cos(radians(S2-30)))) )

Qbin = 0.05
Qmin = floor(min(Qlo))+Qbin  # Bringing in the limit from an exact integer means less whitespace around the plot
Qmax = ceil(max(Qhi))-Qbin
Qbinning = Qmin, Qbin, Qmax
# Done working out Q binning

GroupDetectors(InputWorkspace=input,OutputWorkspace=input,MapFile=r'/SNS/HYSA/shared/adara/HYS_Grouping_2x.xml',Behaviour='Average')
SofQW3(InputWorkspace=input,OutputWorkspace=input,QAxisBinning=Qbinning,EMode='Direct')
Transpose(InputWorkspace=input,OutputWorkspace=output)
