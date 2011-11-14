def loadSpring8(fname,scans,dets):
	
	updat=[1,4,6,7]
	dndat=[2,3,5,8]
	upmon=0.0

	dnmon=0.0
	updata=numpy.zeros((8192,10))
	dndata=numpy.zeros((8192,10))
	upDatNorm=numpy.zeros((8192,10))
	dnDatNorm=numpy.zeros((8192,10))
	upDatNormE=numpy.zeros((8192,10))
	dnDatNormE=numpy.zeros((8192,10))
	tmpvec=numpy.zeros((10))
	tmpdat=numpy.zeros((8192,10))
	usecols=[2,3,4,5,6,7,8,9,10,11]
	ch=range(1,8192+1)
	for j  in range(scans+1):
		print 'reading scan ',j
		for i in updat:
			fnameIn=fname+'_'+str.rjust(str(j+1),3,'0')+'_'+str(i)+'_A.dat'
			#dat=numpy.loadtxt(fnameIn)
			f = open(fnameIn, 'r')
			monline=f.readline()
			upmon=upmon+float(monline.split()[0])
			temp1=float(monline.split()[1])
			tempsamp=float(monline.split()[2])
			
			lline=0
			for line in f:
				nextline=line
				tmpstr= nextline.split()
				#print nextline
				j=0
				for i in usecols:
					tmpvec[j]=float(tmpstr[i])
					j=j+1
				tmpdat[lline,:]=tmpvec
				lline=lline+1
			updata=updata+tmpdat
		for i in dndat:
			fnameIn=fname+'_'+str.rjust(str(j+1),3,'0')+'_'+str(i)+'_B.dat'
			#dat=numpy.loadtxt(fnameIn)
			f = open(fnameIn, 'r')
			monline=f.readline()
			dnmon=dnmon+float(monline.split()[0])
			temp1=float(monline.split()[1])
			tempsamp=float(monline.split()[2])
			
			lline=0
			for line in f:
				nextline=line
				tmpstr= nextline.split()
				#print nextline
				j=0
				for i in usecols:
					tmpvec[j]=float(tmpstr[i])
					j=j+1
				tmpdat[lline,:]=tmpvec
				lline=lline+1
			dndata=dndata+tmpdat
			

		#for i in dndat:
		#	fnameIn=fname+'_'+str.rjust(str(j+1),3,'0')+'_'+str(i)+'_B.dat'
		#	dat=numpy.loadtxt(fnameIn)
		#	dnmon=dnmon+dat[0,0]
		#	temp1=dat[0,1]
		#	tempsamp=dat[0,2]
		#	dat=numpy.loadtxt(fnameIn,skiprows=1,usecols=(2,3,4,5,6,7,8,9,10,11))

			#for det in range(10):
			#	dndata[:,det]=dndata[:,det]+dat[:,det]
				
	for det in dets:
		upDatNorm[:,det]=updata[:,det]*(updata[1,det]/updata[0,det]/upmon)
		dnDatNorm[:,det]=dndata[:,det]*(dndata[1,det]/dndata[0,det]/dnmon)
		upDatNormE[:,det]=numpy.sqrt(updata[:,det])*(updata[1,det]/updata[0,det]/upmon)
		dnDatNormE[:,det]=numpy.sqrt(dndata[:,det])*(dndata[1,det]/dndata[0,det]/dnmon)

	

	#for f in (dets):
	#	CreateWorkspace(fname+'_Updata_Det_'+str(f),ch,upDatNorm[:,f],upDatNormE[:,f])
	#for f in (dets):
	#		ConjoinWorkspaces(InputWorkspace1=fname+'_Updata_Det_'+str(1),InputWorkspace2=fname+'_Updata_Det_'+str(f),CheckOverlapping="0")		
	#for f in (dets):
	#	CreateWorkspace(fname+'_Dndata_Det_'+str(f),ch,dnDatNorm[:,f],dnDatNormE[:,f])
	#for f in (dets):
	#		ConjoinWorkspaces(InputWorkspace1=fname+'_Dndata_Det_'+str(1),InputWorkspace2=fname+'_Dndata_Det_'+str(f),CheckOverlapping="0")


fname='/Users/jon/Desktop/bl08_Nov2011/111104/TbMnO3/TbMnO3_b_13_K_2p5_T'

dets=[1,2,3,4,5,6]

loadSpring8(fname,40,dets)