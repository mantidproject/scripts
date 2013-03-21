# The incident energy needs to be given to DgsReduction because the 'EnergyRequest' log will not contain a value for every chunk after the first one
ei = mtd['__dgs-interim'].getRun().getLogData('Ei').value if mtd.workspaceExists('__dgs-interim') else GetEi(InputWorkspace=input,FixEi='1').getPropertyValue('IncidentEnergy')

DgsReduction(SampleInputWorkspace=input,OutputWorkspace=output,IncidentEnergyGuess=ei,UseIncidentEnergyGuess=1,HardMaskFile=r'/SNS/HYSA/shared/adara/MonsterMask.xml',SofPhiEIsDistribution=0,RejectZeroBackground=0)
