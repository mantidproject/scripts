wsName = 'Live:TotalCounts'
mtd.remove(wsName)
prescript = open('/SNS/HYSA/shared/adara/TotalCountsPreProcessing.py').read()

StartLiveData(UpdateEvery='5',Instrument='HYSPECA',ProcessingScript=prescript,PreserveEvents='0',EndRunBehavior='Stop',OutputWorkspace=wsName)

total_counts = plot(wsName,0).activeLayer().legend().hide()
