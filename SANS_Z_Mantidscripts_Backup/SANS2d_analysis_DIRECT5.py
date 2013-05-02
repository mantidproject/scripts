###########################################################################
# generate direct beam efficiency file for SANS2d       #

# 22/6/10 try new direct beam 5706-add with 1% time bins, at L2=8m
#
# 20/04/10 re-work for M1 direct beam,  first use new add code to make 524-add
# note, changed params to CropWorkspace
###########################################################################
# was   import LOQFunctions
import SANSUtility

#path="//isis/inst$/NDXSANS2D/Instrument/data/cycle_09_1/"
#path = '/home/dmn58364/Desktop/SANS2D_Data/'

#import SANSadd
#   attempt to lose mons did not work, no shape for something, now detector shift not working check instr descr again !
# run 521 is at   19.281 + 6.000 + 0.058 =  23.281 + 2.058
# run 5706 is at   19.281 + 8.000 + 0.058 =  23.281 + 4.058
#SANSadd.add_runs(path, [521,522,523,524])
resultWS="5706-add"
# was front to get rear, swapped id lists in latest instrument description 26/6/9
detbank = 'rear-detector'
# natural bins for M1 are wavbin="0.4412647,-0.01,13.80"
wavbin="0.4412647,-0.01,18.0"
#
monitorspectrum=1 
# natural bins for M2 are                       wavbin="0.4434466,-0.01,13.732"
wavbin="0.4434466,-0.01,18.00"
monitorspectrum=2
#
monitorWS="added_mon"
CropWorkspace(resultWS, OutputWorkspace = monitorWS, StartWorkspaceIndex=str(monitorspectrum - 1), EndWorkspaceIndex=str(monitorspectrum - 1))
ConvertUnits(monitorWS, monitorWS, "Wavelength")
Rebin(monitorWS, monitorWS,wavbin)
#
first=1+8
last=8+192*192
#resultWS="added"
CropWorkspace(resultWS, OutputWorkspace = "lowend", StartWorkspaceIndex=str(first-1), EndWorkspaceIndex=str(last - 1))
resultWS2="lowend"
# pick spectrum for centre, enter number below - beware got these numbers before cropping the monitors off
wksp = mantid.getMatrixWorkspace(resultWS)
instrument = wksp.getInstrument()
det = wksp.getDetector(0)
print det.getPos()
det = wksp.getDetector(1)
print det.getPos()
# M1 says it is at 7.217m,    M2 at 17.937m  (as per RKH's Qranges_tgt2_17Jun09.xls )
det = wksp.getDetector(11485)
print det.getPos()
det = wksp.getDetector(11486)
print det.getPos()
det = wksp.getDetector(11293)
print det.getPos()
det = wksp.getDetector(11294)
print det.getPos()
# then check the coords for the shift  
# mask101Q for 8m has centre at set centre -260.0 -195.0 5.1 5.1 mm
#  pixels above give centre at approx x = -.2754   y = -.1887 m  but then it changed to x  ~ -.3162  after the 8 monitors got cropped off ????
# was -.2342           -.2342+.150=-.0842
xshift = str( +.2754)
yshift = str( +.1887)
zshift = str(4.058)
# DO NOT REPEAT THIS LINE  - IT WILL MOVE IT FURTHER - PERHAPS USE ABSOLUTE COORDINATES ?
MoveInstrumentComponent(resultWS2, detbank, X = xshift, Y = yshift, Z = zshift, RelativePosition="1")
#
# check the result, note use resultWS2 here, and have to subtract 8 from the index!
wksp = mantid.getMatrixWorkspace(resultWS2)
instrument = wksp.getInstrument()
det = wksp.getDetector(11485-8)
print det.getPos()
# final centre coord is ca  -0.150, 0.0, 23.281    (see note above for z coord)
# try convert to wavelength FIRST, before tof is lost !
ConvertUnits(resultWS2, resultWS2, "Wavelength")
# then rebin - as the spectra all have slighlty different paths and hence lambda bins
Rebin(resultWS2,resultWS2+"_wav",wavbin)
#
SumSpectra(resultWS2+"_wav",resultWS2+"_sumbefore")
# decide what radius to include  
#  was 0.035
rmax=str(0.035)
#
# was LOQFunctions.MaskOutsideCylinder(resultWS, rmax)
Rebin(resultWS2,resultWS2+"_mask",wavbin)
SANSUtility.MaskOutsideCylinder(resultWS2+"_mask", rmax)
SumSpectra(resultWS2+"_mask",resultWS2+"_sum")
# take ratio to M2
Divide(resultWS2+"_sum",monitorWS,"dirbeam")

