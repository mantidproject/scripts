def simpleRebin(ppz,ppznew,data):
	for i in range(len(ppznew)-1):
		print i
		tmp=where(ppz>ppznew[i])#need the first element of this
		tmp2=where(ppz<=ppznew[i+1])#need the last element of this
		
		lolim=int(tmp[0][0])
		hilim=int(tmp2[0][size(tmp2)-1])
		print lolim, hilim
		sum(data[lolim:hilim])/(hilim-lolim);