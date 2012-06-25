autows='__auto_ws'
CloneWorkspace(InputWorkspace=input, OutputWorkspace=autows)
MaskDetectors(Workspace=autows, MaskedWorkspace='BASIS_MASK')
ModeratorTzero(InputWorkspace=autows,OutputWorkspace=autows)
ConvertUnits(InputWorkspace=autows, OutputWorkspace=autows, Target='DeltaE', EMode='Indirect')
CorrectKiKf(InputWorkspace=autows, OutputWorkspace=autows,EMode='Indirect')
Rebin(InputWorkspace=autows, OutputWorkspace=autows, Params='-0.12,0.0004,0.12')
GroupDetectors(InputWorkspace=autows, OutputWorkspace=output, MapFile='/SNS/BSS/shared/autoreduce/BASIS_Grouping.xml', Behaviour='Sum')
