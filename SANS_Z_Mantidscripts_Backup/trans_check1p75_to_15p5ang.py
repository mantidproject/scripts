# SANS2d    
#  26/7/12 check transmissions with M4
from mantidsimple import *
from ISISCommandInterface import *
path="\\\isis\inst$\NDXSANS2D\Instrument\data\cycle_12_1\\"
#path="U:\\processed\\"
pfix="sans2d"
sfix=".nxs"
#
trans=13391
direct=13388
#
specincid=1
spectrans=4
backincid1=35000
backincid2=98000
backtrans1=10
backtrans2=100
#zshift=-100
#
#trans=7991
#direct=7970
#specincid=2
#spectrans=3
#backincid1=85000
#backincid2=98000
#backtrans1=85000
#backtrans2=98000
wavbins="1.75,0.125,15.5"
#
#	Load the monitors only and divide by uamphr
transrun=str(trans)
directrun=str(direct)
#
# read direct beam   ================================================================
nzeros=8-len(directrun)
fpad=""
for ii in range(nzeros):
   fpad+="0"
filename=path+pfix+fpad+directrun+sfix
print "reading file:   "+filename
check=LoadISISNexus(Filename=filename,OutputWorkspace=directrun,SpectrumMin=1,SpectrumMax=1).workspace()
print "uamphr = " + str(check.getSampleDetails().getProtonCharge())
m1=LoadISISNexus(Filename=filename,OutputWorkspace=directrun,SpectrumMin=1,SpectrumMax=4)

MoveInstrumentComponent(Workspace=directrun,ComponentName='rear-detector',X='0.10639999999999999',Y='0.1691',Z='2.0582880859400001')
MoveInstrumentComponent(Workspace=directrun,ComponentName='some-sample-holder',Z='0.052999999999999999')
MoveInstrumentComponent(Workspace=directrun,ComponentName='monitor4',Z='-4.7607119140599998')
#NormaliseByCurrent(run+"_norm",run+"_norm")
#
# incident spectrum
CropWorkspace(directrun, directrun+"_incid", StartWorkspaceIndex=specincid-1, EndWorkspaceIndex=specincid-1)
ConvertToDistribution(directrun+"_incid")
#                                        integrate long time region to average flat background
Integration(directrun+"_incid","integral",backincid1,backincid2)
flat = mtd["integral"].readY(0)[0]/(backincid2-backincid1)
CreateSingleValuedWorkspace("scalar",flat)
#                                          subtract flat background
Minus(directrun+"_incid","scalar",directrun+"_incid_sub")
#  trans spectrum
CropWorkspace(directrun, directrun+"_trans", StartWorkspaceIndex=spectrans-1, EndWorkspaceIndex=spectrans-1)
Integration(directrun+"_incid","integral",backtrans1,backtrans2)
ConvertToDistribution(directrun+"_trans")
Integration(directrun+"_trans","integral",backtrans1,backtrans2)
flat = mtd["integral"].readY(0)[0]/(backtrans2-backtrans1)
CreateSingleValuedWorkspace("scalar",flat)
Minus(directrun+"_trans","scalar",directrun+"_trans_sub")
#
ConvertUnits(directrun+"_incid_sub",directrun+"_incid_wav","Wavelength")
ConvertUnits(directrun+"_trans_sub",directrun+"_trans_wav","Wavelength")
#
# read trans run ============================================================================
nzeros=8-len(transrun)
fpad=""
for ii in range(nzeros):
   fpad+="0"
filename=path+pfix+fpad+transrun+sfix
print "reading file:   "+filename
check=LoadNexus(Filename=filename,OutputWorkspace=transrun,SpectrumMin=1,SpectrumMax=1).workspace()

print "uamphr = " + str(check.getSampleDetails().getProtonCharge())
m1=LoadISISNexus(Filename=filename,OutputWorkspace=transrun,SpectrumMin=1,SpectrumMax=4)
MoveInstrumentComponent(Workspace=transrun,ComponentName='rear-detector',X='0.10639999999999999',Y='0.1691',Z='2.0582880859400001')
MoveInstrumentComponent(Workspace=transrun,ComponentName='some-sample-holder',Z='0.052999999999999999')
MoveInstrumentComponent(Workspace=transrun,ComponentName='monitor4',Z='-4.7607119140599998')
#NormaliseByCurrent(run+"_norm",run+"_norm")
# incident spectrum
CropWorkspace(transrun, transrun+"_incid", StartWorkspaceIndex=specincid-1, EndWorkspaceIndex=specincid-1)
ConvertToDistribution(transrun+"_incid")
Integration(transrun+"_incid","integral",backincid1,backincid2)
flat = mtd["integral"].readY(0)[0]/(backincid2-backincid1)
CreateSingleValuedWorkspace("scalar",flat)
Minus(transrun+"_incid","scalar",transrun+"_incid_sub")
# trans spectrum
CropWorkspace(transrun, transrun+"_trans", StartWorkspaceIndex=spectrans-1, EndWorkspaceIndex=spectrans-1)
ConvertToDistribution(transrun+"_trans")
Integration(transrun+"_trans","integral",backtrans1,backtrans2)
flat = mtd["integral"].readY(0)[0]/(backtrans2-backtrans1)
#
# TEST
#flat=0.0
#
CreateSingleValuedWorkspace("scalar",flat)
Minus(transrun+"_trans","scalar",transrun+"_trans_sub")
# convert to wavelength
ConvertUnits(transrun+"_incid",transrun+"_incid_wav","Wavelength")
ConvertUnits(transrun+"_trans",transrun+"_trans_wav","Wavelength")
#
# rebin everything
InterpolatingRebin(directrun+"_incid_wav",directrun+"_incid_wav",wavbins)
InterpolatingRebin(directrun+"_trans_wav",directrun+"_trans_wav",wavbins)
InterpolatingRebin(transrun+"_incid_wav",transrun+"_incid_wav",wavbins)
InterpolatingRebin(transrun+"_trans_wav",transrun+"_trans_wav",wavbins)
#
Divide(directrun+"_trans_wav",directrun+"_incid_wav","direct_trans_by_incid")
Divide(transrun+"_trans_wav",transrun+"_incid_wav","trans_trans_by_incid")
#
Divide("trans_trans_by_incid","direct_trans_by_incid",transrun+"_by_"+directrun)

# FITH INTEGRATION ========================================================================================
#	Integration(run+"_m3_div","integral",str(7.0),str(12.0))
#	try:
#				cr5 = mtd["integral"].readY(0)[0]
