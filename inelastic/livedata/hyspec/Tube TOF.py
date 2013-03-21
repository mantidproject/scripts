wsName = 'Live:TubeTOF'
mtd.remove(wsName)
prescript = open('/SNS/HYSA/shared/adara/TubeTOFPreProcessing.py').read()

StartLiveData(UpdateEvery='5',Instrument='HYSPECA',ProcessingScript=prescript,EndRunBehavior='Stop',OutputWorkspace=wsName)

workspace_mtx = importMatrixWorkspace(wsName)
graph2D = workspace_mtx.plotGraph2D().activeLayer()
graph2D.logColor()
