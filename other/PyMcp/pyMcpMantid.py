#load a pointer to a file
import specfile as sp
import numpy as np
#from pylab import *
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
#import MantidFramework 
#MantidFramework.mtd.initialise()
#from mantidsimple import *
sqrt=np.sqrt
cos=np.cos
arccos=np.arccos
ones=np.ones
zeros=np.zeros
array=np.array
exp=np.exp
where=np.where
size=np.size
pi=np.pi
flipud=np.flipud
class spectrum2D():
	def __init__(self):
		self.initilise()
	
	def initilise(self):
		self.updata=0
		self.dndata=0
		self.uperr=0
		self.dnerr=0
		self.mcpraw=0
		self.compton_up=[]
		self.compton_dn=[]
		self.diode_up=[]
		self.diode_dn=[]
		self.fluo_up=[]
		self.fluo_dn=[]
		self.detectors=[]
		self.prefix=''
		self.xlabel='Channel'
		self.xscale=[]
		self.numDets=0
		self.mca_size=0
	
	def plotMCA(self,dir,type):
		if type=='single':
			if dir =='up':
				for i in range(len(self.updata)):
					plot(self.xscale,self.updata[i][:],'o')
					waitforbuttonpress()
					cla()
			if dir =='dn':
				for i in range(len(self.dndata)):
					plot(self.xscale,self.dndata[i][:],'o')
					waitforbuttonpress()
					cla()
		if type=='all':
			if dir =='up':
				for i in range(len(self.updata)):
					plot(self.xscale,self.updata[i][:],'o')
			if dir =='dn':
				for i in range(len(self.dndata)):
					plot(self.xscale,self.dndata[i][:],'o')
	
class spectrum1D():
	def __init__(self):
		self.initilise()
	
	def initilise(self):
		self.data=0
		self.err=0
		self.detectors=[]
		self.prefix=''
		self.xlabel=''
		self.direction=''
		self.xscale=[]
		self.numDets=0
		self.mca_size=0
		self.mcp=[]
		
	def plot(self,dir,type):
			plot(self.xscale,self.data,'o')
	
	def fitGauss(self):
		y= self.data
		x = self.xscale
		## Parametric function: 'v' is the parameter vector, 'x' the independent varible

		fp = lambda v, x: (v[0]+(v[1]**2)*exp(-(x-v[2])**2/(2*v[3]**2))) 
		
		
		## Error function
		e = lambda v, x, y: (fp(v,x)-y)
		
		
		## Initial parameter value guess from moments
		
		cen = sum(x*y)/sum(y)
		width = sqrt(abs(sum((x-cen)**2*y)/sum(y)))
		max = y.max()
		min=y.min()
		v0 = [min,max, cen, width]
		
		## Fitting
		v, success = leastsq(e, v0, args=(x,y), maxfev=10000)
		
		#plot(x, fp(v,x))
		print 'Background= ', v[0], '; peak to background = ',v[1]/v[0]
		print 'Centre = ',v[2]
		print 'FWHM = ',v[3]
		return v

class mcp1D():
	def __init__(self):
		self.initilise()
	
	def initilise(self):
		
		self.flipdata=[]
		self.fliperr=[]
		self.unflipdata=[]
		self.unflippederr=[]
		self.rawerr=[]
		self.datraw=[]
		self.detectors=[]
		self.prefix=[]
		self.xlabel=[]
		self.rebin=[]
		self.pz=[]
		self.rawPz=[]
		self.numDets=[]
		self.mca_size=[]
		
	def plot(self):
		errorbar(self.pz,self.flipdata,yerr=self.fliperr,xerr=None,marker='o')
		title('Flipped & Rebined MCP from '+self.prefix)
	
	def plotDetectors(self):
		for f in (self.detectors):
			figure()
			errorbar(self.pz,self.datraw[f-1][:],yerr=self.rawerr[f-1][:],xerr=None,marker='o')
			title('Detector element'+str(f))
		
	def scale(self,factor):
		self.flipdata=self.flipdata*factor
		self.fliperr=self.fliperr*factor
		
	def integrate(self):
		inte=np.trapz(self.pz,self.flipdata)
		print 'integral = ',inte
	
	def McpToMantid(self,name):
		CreateWorkspace(name,self.pz,self.flipdata,self.fliperr)
	
	def RawMcpToMantid(self,name):
		#print 'PZ', self.pz.shape
		#print 'data', self.datraw.shape
		#print 'err', self.rawerr.shape
		rawPz=numpy.nan_to_num(self.rawPz)
		rawERR=numpy.nan_to_num(self.rawerr)
		datraw=numpy.nan_to_num(self.datraw)
		
		for f in (self.detectors):
			CreateWorkspace(name+'Det_'+str(f),rawPz[f-1],datraw[f-1],rawERR[f-1])
			Rebin(name+'Det_'+str(f),name+'Det_'+str(f),'-20,.05,20')
		for f in (self.detectors):
			ConjoinWorkspaces(InputWorkspace1=name+'Det_'+str(self.detectors[0]),InputWorkspace2=name+'Det_'+str(f),CheckOverlapping="0")
		
		
	def fitGauss(self):
		y= self.data
		x = self.xscale
		## Parametric function: 'v' is the parameter vector, 'x' the independent varible

		fp = lambda v, x: (v[0]+(v[1]**2)*exp(-(x-v[2])**2/(2*v[3]**2))) 
		
		
		## Error function
		e = lambda v, x, y: (fp(v,x)-y)
		
		
		## Initial parameter value guess from moments
		
		cen = sum(x*y)/sum(y)
		width = sqrt(abs(sum((x-cen)**2*y)/sum(y)))
		max = y.max()
		min=y.min()
		v0 = [min,max, cen, width]
		
		## Fitting
		v, success = leastsq(e, v0, args=(x,y), maxfev=10000)
		
		#plot(x, fp(v,x))
		print 'Background= ', v[0], '; peak to background = ',v[1]/v[0]
		print 'Centre = ',v[2]
		print 'FWHM = ',v[3]
		return v

class mcp:
	def __init__(self,location):
		self.location=location
		self.initilise()
	
	def initilise(self):
		if self.location =='esrf':
			print 'A python class for reading ESRF spec files and perfoming analysis'
			print 'of magnetic compton scattering data'
			self.updata=0
			self.dndata=0
			self.uperr=0
			self.dnerr=0
			self.mcpraw=0
			self.mcpInterp=[]
			self.compton_up=[]
			self.compton_dn=[]
			self.diode_up=[]
			self.diode_dn=[]
			self.fluo_up=[]
			self.fluo_dn=[]
			self.prefix=''
			self.xlabel='Channel'
			print 'Setting xlabel to ',self.xlabel
			self.xscale=[]
			self.ch=[]
			self.en=[]
			self.numDets=13
			print 'Setting number of detectors to ',self.numDets
			self.mca_size=4096
			print 'Setting default mca size to ',self.mca_size, ' channels'
			self.detectors=[1,3,4,5,7,8,9,10,11,12,13]
			print 'Setting detectors to be includes to ',self.detectors
			self.diode_chan='diode2'
			print 'Setting diode normalisation column to ',self.diode_chan
			self.current='srcur'
			print 'Setting Ring current column to ',self.current
			self.compton='compton'
			print 'Setting compton counter column to ',self.compton
			self.fluo_chan='snlines'
			print 'Setting fluo counter column to ',self.fluo_chan
			self.real_time='Seconds'
			print 'Setting real time column to ',self.real_time
			self.live_time='DeadTime'
			print 'Setting Live time column to ',self.live_time
			self.num_dets=13
		elif self.location =='spring8':
			print 'A python class for perfoming analysis'
			print 'of magnetic compton scattering data taken on Bl08W at spring8'
			self.updata=0
			self.dndata=0
			self.uperr=0
			self.dnerr=0
			self.mcpraw=0
			self.mcpInterp=[]
			self.compton_up=[]
			self.compton_dn=[]
			self.diode_up=[]
			self.diode_dn=[]
			self.fluo_up=[]
			self.fluo_dn=[]
			self.prefix=''
			self.xlabel='Channel'
			print 'Setting xlabel to ',self.xlabel
			self.xscale=[]
			self.ch=[]
			self.en=[]
			self.numDets=10
			print 'Setting number of detectors to ',self.numDets
			self.mca_size=8192
			print 'Setting default mca size to ',self.mca_size, ' channels'
			self.detectors=[1,2,3,4,5,6,7,8,9,10]
			print 'Setting detectors to be includes to ',self.detectors
			self.diode_chan='diode2'
			print 'Setting diode normalisation column to ',self.diode_chan
			self.current='srcur'
			print 'Setting Ring current column to ',self.current
			self.compton='compton'
			print 'Setting compton counter column to ',self.compton
			self.fluo_chan='snlines'
			print 'Setting fluo counter column to ',self.fluo_chan
			self.real_time='Seconds'
			print 'Setting real time column to ',self.real_time
			self.live_time='DeadTime'
			print 'Setting Live time column to ',self.live_time
			self.num_dets=10
		else:
			print 'where were the data taken'
	def load():
		if self.location=='esrf':
			self.loadspecfile()
		if self.location=='spring8':
			self.loadspring8data

	def RawDataToMantid(self,name):
		for f in (self.detectors):
			CreateWorkspace(name+'_Updata_Det_'+str(f),self.ch[f-1],np.nan_to_num(self.updata[f-1]),np.nan_to_num(self.uperr[f-1]))
				
		for f in (self.detectors):
			ConjoinWorkspaces(InputWorkspace1=name+'_Updata_Det_'+str(self.detectors[0]),InputWorkspace2=name+'_Updata_Det_'+str(f),CheckOverlapping="0")
		
		for f in (self.detectors):
			CreateWorkspace(name+'_Dndata_Det_'+str(f),self.ch[f-1],np.nan_to_num(self.dndata[f-1]),np.nan_to_num(self.dnerr[f-1]))
				
		for f in (self.detectors):
			ConjoinWorkspaces(InputWorkspace1=name+'_Dndata_Det_'+str(self.detectors[0]),InputWorkspace2=name+'_Dndata_Det_'+str(f),CheckOverlapping="0")
		
			
	def plotMCA(self,dir,type):
		if type=='single':
			if dir =='up':
				for i in range(len(self.updata)):
					plot(self.xscale[i][:],self.updata[i][:],'o')
					waitforbuttonpress()
					cla()
			if dir =='dn':
				for i in range(len(self.dndata)):
					plot(self.xscale[i][:],self.dndata[i][:],'o')
					waitforbuttonpress()
					cla()
		if type=='all':
			if dir =='up':
				for i in range(len(self.updata)):
					plot(self.xscale[i][:],self.updata[i][:],'o')
			if dir =='dn':
				for i in range(len(self.dndata)):
					plot(self.xscale[i][:],self.dndata[i][:],'o')

	def plotMCP(self,type):
		
		if type=='single':
			for i in range(len(self.mcpraw)):
				plot(self.xscale[i][:],self.mcpraw[i][:],'o')
				waitforbuttonpress()
				cla()
		if type=='all':
			for i in range(len(self.mcpraw)):
				plot(self.xscale[i][:],self.mcpraw[i][:],'o')
				
	def plotdiode(self,dir):
		if dir =='up':
			plot(self.diode_up,'o')
			show()
		if dir =='dn':
			plot(self.diode_dn,'o')
					
	def plotcompton(self,dir):
		if dir =='up':
			plot(self.compton_up,'o')
			show()
		if dir =='dn':
			plot(self.compton_dn,'o')
			show()
			
	def plotfluo(self,dir):
		if dir =='up':
			plot(self.fluo_up,'o')
			show()
		if dir =='dn':
			plot(self.fluo_dn,'o')
			show()

	def calibrate_energy(self,e1,e2,p1,p2,window):
		#gives the calibration coefficients a and b where
		#	Energy = a*Channel + b
		#
		# b1 and b2 are calculated respectively using e1,c1 and e2,c2 and their values 
		# should be identical.
		#
		# e1 = first energy value
		# e2 = second energy value
		# c1 = channel corresponding to first energy peak
		# c2 = channel corresponding to second energy peak
		#
		if self.xlabel != 'Channel':
			print 'Xscale should be channel for this operation use setXunits to change'
			return
				
		pos1=ones(self.numDets)
		pos2=ones(self.numDets)
		ch=np.arange(1,self.mca_size+1)
		en=np.zeros([self.numDets,self.mca_size])
		for f in (self.detectors):
			print '################################################################'
			print 'Detector-- ', f
			#fit first peak
			tmp=self.spec1D([p1-window,p1+window],f,'up')
			v=tmp.fitGauss()
			pos1=v[2]
			#fit second peak
			tmp=self.spec1D([p2-window,p2+window],f,'up')
			v=tmp.fitGauss()
			pos2=v[2]	
			a = (e2-e1)/(pos2-pos1)
			b = e1-a*pos1
			print a,b
			e = a*ch+b
			en[f-1][:]=e
			print '################################################################'
		self.en=en
		self.xscale=en
		self.xlabel='Energy [keV]'
	
	def ConvertToPz(self,eC,Ei,window,**kwargs):
	
		if self.xlabel != 'Energy [keV]':
			print 'Xscale should be Energy for this operation use setXunits to change'
			return
		
		if kwargs.has_key('fixei'):
			Ei = kwargs.get('fixei')
			print 'Setting fixei to ', Ei
		
		#create a n detetector pzscale
		Pz=np.zeros([self.numDets,self.mca_size])
		print '########Convert to Pz##############'
		for f in (self.detectors):
			print '################################################################'
			print 'Detector-- ', f
			#fit Compton peak
			tmp=self.spec1D([eC-window,eC+window],f,'up')
			v=tmp.fitGauss()
			PosCompton=v[2]
			#fit second peak
			tmp=self.spec1D([Ei-window,Ei+window],f,'up')
			v=tmp.fitGauss()
			if kwargs.has_key('fixei'):
				PosEi=Ei
			else:
				PosEi=v[2]	
			ScatAngle=arccos(1.0+511.0*((1.0/PosEi)-(1.0/PosCompton)))/pi*180
			print 'Scattering angle for detector ',f,' = ',ScatAngle,'deg'
			#(w1=Ei,w2=energy_vec,th=theta)
			w1=PosEi
			w2=self.xscale[f-1]
			th=ScatAngle
			P=0.0
			m=511.003		# Mass in natural units
			th=th/180.0*pi		# Angle to radians
			alp=1.0/137.036		# Fine structure constant
			r0=2.8179e-15		# Electron radius
			q=sqrt(w1**2+w2**2-2*w1*w2*cos(th))			# Momentum transfer
			pz=q/2-(w1-w2)*sqrt(1.0/4.0+m**2/(2.0*w1*w2*(1.0-cos(th))))	# In natural units
			E=sqrt(m**2+pz**2)
			A=((w1-w2)*E-w1*w2*(1-cos(th)))/q
			D=(w1-w2*cos(th))*A/q
			R=w1*(E-D)
			R2=R-w1*w2*(1-cos(th))
			chi=R/R2+R2/R+(2*m**2*(1/R-1/R2)+m**4*(1/R-1/R2)**2)*(1-P)
			cf=2*w1*q*E/(m**2*r0**2*w2*chi)
			cf=cf*(1e-28*(m*alp))		# Cross section now in barns/atom/keV/srad
			pz=pz/(m*alp)	       		# pz to atomic units (a.u.)
			
			Pz[f-1][:]=pz
			#correct data for xsection
			self.updata[f-1][:]=self.updata[f-1][:]*cf
			self.dndata[f-1][:]=self.dndata[f-1][:]*cf
			self.uperr[f-1][:]=self.uperr[f-1][:]*cf
			self.dnerr[f-1][:]=self.dnerr[f-1][:]*cf
			print '########################################################'
		self.Pz=Pz
		self.xscale=Pz
		self.xlabel='Pz [au]'
	
	def interpMCP(self,rebin):
		if self.xlabel != 'Pz [au]':
			print 'Xscale should be Pz for this operation use setXunits to change'
			return
		#rebin onto a constant pz scale using interpolation
		PPZ=np.linspace(rebin[0],rebin[2],rebin[2]-rebin[0]/rebin[1])
		mcpInterp=np.zeros([self.numDets,len(PPZ)])
		mcpInterpERR=np.zeros([self.numDets,len(PPZ)])
		for f in (self.detectors):
			x=self.xscale[f-1]
			y=self.updata[f-1]-self.dndata[f-1]
			yerr=sqrt((self.uperr[f-1]**2)+(self.dnerr[f-1]**2))
			iterp=interp1d(x,y)
			interpTMP=iterp(PPZ)
			
			iterpERR=interp1d(x,yerr)
			interpTMPERR=iterpERR(PPZ)
			mcpInterp[f-1][:]=interpTMP
			mcpInterpERR[f-1][:]=interpTMPERR
		mcpSum=sum(mcpInterp,0)
		mcpSumERR=sqrt(sum(mcpInterpERR**2,0))
		errorbar(PPZ,mcpSum,yerr=mcpSumERR,xerr=None,marker='o')
		title('Interpolated rebin MCP from '+self.prefix)
		xlabel=self.xlabel
		mcpSym=(mcpSum+flipud(mcpSum))/2
		mcpSymERR=sqrt(mcpSumERR**2+flipud(mcpSumERR**2))/2
		figure()
		errorbar(PPZ,mcpSym,yerr=mcpSymERR,xerr=None,marker='o')
		title('Flipped Interpolated rebin MCP from '+self.prefix)
		xlabel=self.xlabel
		
	def rebinMCP(self,rebin):
		if self.xlabel != 'Pz [au]':
			print 'Xscale should be Pz for this operation use setXunits to change'
			return
		#rebin onto a constant pz scale using a basic rebin
		#PPZ=np.linspace(rebin[0],rebin[2],rebin[2]-rebin[0]/rebin[1])
		PPZ=numpy.arange(rebin[0],rebin[2],rebin[1])
		mcpReb=np.zeros([self.numDets,len(PPZ)-1])
		mcpRebERR=np.zeros([self.numDets,len(PPZ)-1])
		print mcpReb.shape
		tmpdat=np.zeros(len(PPZ)-1)
		tmpdatERR=np.zeros(len(PPZ)-1)
		for f in (self.detectors):
			data=self.updata[f-1]-self.dndata[f-1]
			ERR=sqrt((self.uperr[f-1]**2)+(self.dnerr[f-1]**2))
			for i in range(len(PPZ)-1):
			
				tmp=where(self.Pz[f-1]>PPZ[i])#need the first element of this
				tmp2=where(self.Pz[f-1]<=PPZ[i+1])#need the last element of this
				lolim=int(tmp[0][0])
				hilim=int(tmp2[0][size(tmp2)-1])
				tmpdat[i]=sum(data[lolim:hilim])/(hilim-lolim);
				tmpdatERR[i]=sqrt(sum(ERR[lolim:hilim]**2))/(hilim-lolim);
			mcpReb[f-1][:]=tmpdat
			mcpRebERR[f-1][:]=tmpdatERR
		
		rawmcpdata=np.zeros([self.num_dets,self.mca_size])
		rawmcpERR=np.zeros([self.num_dets,self.mca_size])

		for f in (self.detectors):
			rawmcpdata[f-1,:]=self.updata[f-1]-self.dndata[f-1]
			rawmcpERR[f-1,:]=sqrt((self.uperr[f-1]**2)+(self.dnerr[f-1]**2))
		
		dpz=(PPZ[1]-PPZ[0])/2.0
		PPZnew=PPZ[0:len(PPZ)-1]+dpz
		self.PPZreb=PPZnew
		mcpRebSum=sum(mcpReb,0)
		mcpRebSumERR=sqrt(sum(mcpRebERR**2,0))
		
		#errorbar(PPZnew,mcpRebSum,yerr=mcpRebSumERR,xerr=None,marker='o')
		#title('Rebined MCP from '+self.prefix)
		#xlabel=self.xlabel
		mcpSym=(mcpRebSum+flipud(mcpRebSum))/2
		mcpSymERR=sqrt(mcpRebSumERR**2+flipud(mcpRebSumERR**2))/2
		#figure()
		#errorbar(PPZnew,mcpSym,yerr=mcpSymERR,xerr=None,marker='o')
		#title('Flipped & Rebined MCP from '+self.prefix)
		#xlabel=self.xlabel
		#CreateWorkspace('test',PPZnew,mcpSym,mcpSymERR)
		
		dat=mcp1D()
		print mcpReb.shape
		dat.flipdata=mcpSym
		dat.fliperr=mcpSymERR
		dat.unflipdata=mcpRebSum
		dat.unflippederr=mcpRebSumERR
		dat.rawerr=rawmcpERR
		dat.datraw=rawmcpdata
		dat.detectors=self.detectors
		dat.prefix=self.prefix
		dat.xlabel=self.xlabel
		dat.rebin=rebin
		dat.pz=PPZnew
		dat.rawPz=self.Pz
		dat.numDets=self.numDets
		dat.mca_size=self.mca_size
		
		return dat
		
	def setXUnits(self,set):
		if set =='Channel'or set=='channel':
			print 'setting x Units to channel'
			self.xscale=self.ch
			self.xlabel='Channel'
		if set =='Energy' or set =='energy':
			print 'setting x Units to energy'
			try:
				self.xscale=self.en
				#self.mca_size=4096
				self.xlabel='Energy [keV]'
			except NameError:
				print 'No Energy scale'
			
		if set =='Pz' or set =='pz':
			print 'setting x Units to Pz'
			try:
				self.xscale=self.Pz
				self.xlabel='Pz [au]'
			except NameError:
				print 'No momentum scale'
	def replaceEnergyScale(self,newEn):
		#replaces the self.en with newEn save haveing to recalibrate all the time
		print 'setting ',newEn ,'as the energy scale for this object'
		self.en=newEn
	def getScanNum(self):
		f=1
		fname=self.prefix+'%02d'%f
		print fname
		file=sp.Specfile(fname)
		num_scan=(file.scanno()/4)*4
		print 'number of scans', num_scan,'; with ', num_scan/4,' sets of ABBA cycles'
		return num_scan, num_scan/4
		
	def spec1D(self,lims,det,dir):
		out_dat=spectrum1D()
		
		if self.xlabel=='Channel':
			if dir=='up':
				out_dat.data=self.updata[det-1][lims[0]:lims[1]]
				out_dat.xscale=array(self.xscale[det-1][lims[0]:lims[1]])
				out_dat.mca_size=len(self.xscale[det-1][lims[0]:lims[1]])
			if dir=='dn':
				out_dat.data=self.dndata[det-1][lims[0]:lims[1]]
				out_dat.xscale=array(self.xscale[det-1][lims[0]:lims[1]])
				out_dat.mca_size=len(self.xscale[det-1][lims[0]:lims[1]])
			
		if self.xlabel=='Energy [keV]':
			tmp=where(self.xscale[det-1]>=lims[0])#need the first element of this
			tmp2=where(self.xscale[det-1]<=lims[1])#need the last element of this
			lolim=int(tmp[0][0])
			hilim=int(tmp2[0][size(tmp2)-1])
			if dir=='up':
				out_dat.data=self.updata[det-1][lolim:hilim]
				out_dat.xscale=array(self.xscale[det-1][lolim:hilim])
				out_dat.mca_size=len(self.xscale[det-1][lolim:hilim])
			if dir=='dn':
				out_dat.data=self.dndata[det-1][lolim:hilim]
				out_dat.xscale=array(self.xscale[det-1][lolim:hilim])
				out_dat.mca_size=len(self.xscale[det-1][lolim:hilim])
				
		if self.xlabel=='Pz [au]':
			tmp=where(self.xscale[det-1]>=lims[0])#need the first element of this
			tmp2=where(self.xscale[det-1]<=lims[1])#need the last element of this
			lolim=int(tmp[0][0])
			hilim=int(tmp2[0][size(tmp2)-1])
			
			out_dat.mcp=self.updata[det-1][lolim:hilim]-self.dndata[det-1][lolim:hilim]
			out_dat.xscale=array(self.xscale[det-1][lolim:hilim])
			out_dat.mca_size=len(self.xscale[det-1][lolim:hilim])
			
				
			
		out_dat.err=0
		out_dat.detectors=det
		out_dat.prefix=self.prefix
		out_dat.xlabel=self.xlabel
		out_dat.direction=dir
		
		out_dat.numDets=0
		#if self.xlabel=='Pz [au]':
			#plot(out_dat.xscale,out_dat.mcp)
		#else:
			#plot(out_dat.xscale,out_dat.data)

		return out_dat
	
	def cycleFromScans(self,scans_included):
		list(scans_included)
		sc1=array([1,2,3,4])
		listscan=[]
		for sc in scans_included:
			sc=list(sc1+(sc-1)*4)
			listscan.extend(sc)
			
       		print 'Scans ',str(listscan),' correspond to ABBA cycles' ,str(scans_included)
       		return listscan
	def McpsPerScan(self,rebin):
		#assumes a preloaded complete dataset then runs through ABBA cycle by ABBA cycle
		#to give the mcps there must be an existing pz scale for all channels
		numscan,numcycle=self.getScanNum()
		for cycle in range(1,numcycle+1):
			print cycle
			ss='['+str(cycle)+']'
			cur_scan=eval(ss)
			print type(cur_scan)
			
			self.loadspecfile(scans=cur_scan)
			self.setXUnits('pz')
			self.rebinMCP(rebin)
		
	
	def loadspecfile(self,*args,**kwargs):
		#out_dat=mcp()
		#print type(out_dat)
		#out_dat.numDets=13
		#out_dat.mca_size=4096
		#out_dat.prefix=prefix
		#num_dets=13
		#mca_size=4096
		perscan=0
		if len(args)==1 and type(args[0])==str:
			self.prefix=args[0]
		
		if kwargs.has_key('scans'):
			scans_included=kwargs.get('scans')
			print scans_included
			
			perscan=1
			listscan=self.cycleFromScans(scans_included)
			

			
		updat=np.zeros([self.num_dets,self.mca_size])
		dndat=np.zeros([self.num_dets,self.mca_size])
		uperr=np.zeros([self.num_dets,self.mca_size])
		dnerr=np.zeros([self.num_dets,self.mca_size])
		ch=np.zeros([self.num_dets,self.mca_size])
		xscale=np.zeros([self.num_dets,self.mca_size])
		monitor_up=[]
		monitor_dn=[]
		current_up=[]
		current_dn=[]
		compton_up=[]
		compton_dn=[]
		fluo_up=[]
		fluo_dn=[]
		#number of scans
		
		#prefix='/Users/jon/Desktop/HE3579/id15/EuFeCoAs_2K_2T_b_'
		detectors=self.detectors
		diode_chan=self.diode_chan
		current=self.current#'srcur'
		compton=self.compton#'compton'
		fluo_chan=self.fluo_chan#'snlines'
		real_time=self.real_time#'Seconds'
		live_time=self.live_time#'DeadTime'
		for f in (detectors):
			fname=self.prefix+'%02d'%f
			print fname
			file=sp.Specfile(fname)
			num_scan=(file.scanno()/4)*4
			print 'number of scans', num_scan,'; with ', num_scan/4,' sets of ABBA cycles'
			#select scan
			if perscan==0:
				scans_run=range(1,num_scan+1)
			if perscan == 1:
				scans_run=listscan
				
			for i in scans_run:
				scan=file.select(str(i))
				#number of lines in scan data
				#print 'number of lines', scan.lines()
				#number of mcas in scan
				#print 'number of MCAs', scan.nbmca()
				#note that the current macro will have twice the number of mcas as scans
				# only half are useful depending on the magnet direction
				#names of the counters in the scan returns a list
				counters=scan.alllabels()
				
				magnet=scan.datacol('Magnet')
				if magnet.mean()==1:
					print 'magnet is up'
					dir=1
						
				if magnet.mean()==-1:
					print 'magnet is down'
					dir=-1
				
				if dir == 1:
					#returns the scan data for column labled diode2
					diode=scan.datacol(diode_chan)
					monitor_up.extend(diode)
					scan.datacol(current)
					current_up.extend(scan.datacol(current))
					compton_up.extend(scan.datacol(compton))
					fluo_up.extend(scan.datacol(fluo_chan))
					RT=scan.datacol(real_time)
					LT=scan.datacol(live_time)
					for i in range(1,scan.nbmca()+1,2):
						mon_ind=0
						tmpup=(scan.mca(i)/diode[mon_ind])*(RT[mon_ind]/LT[mon_ind])
						tmpupERR=(sqrt(scan.mca(i))/diode[mon_ind])*(RT[mon_ind]/LT[mon_ind])
						try:
							updat[f-1][:]=updat[f-1][:]+tmpup
							uperr[f-1][:]=sqrt(uperr[f-1][:]**2+tmpupERR**2)
						except:
							#catches instances where the mca is read as 4095 length!!
							updat[f-1][0:len(tmpup)]=updat[f-1][0:len(tmpup)]+tmpup[0:len(tmpup)]
							uperr[f-1][0:len(tmpup)]=sqrt(uperr[f-1][0:len(tmpup)]**2+tmpupERR[0:len(tmpupERR)]**2)
						mon_ind=mon_ind+1
				if dir == -1:
					#returns the scan data for column labled diode2
					diode=scan.datacol(diode_chan)
					monitor_dn.extend(diode)
					current_dn.extend(scan.datacol(current))
					compton_dn.extend(scan.datacol(compton))
					fluo_dn.extend(scan.datacol(fluo_chan))
					RT=scan.datacol(real_time)
					LT=scan.datacol(live_time)
					for i in range(2,scan.nbmca()+1,2):
						mon_ind=0
						tmpdn=(scan.mca(i)/diode[mon_ind])*(RT[mon_ind]/LT[mon_ind])
						tmpdnERR=(sqrt(scan.mca(i))/diode[mon_ind])*(RT[mon_ind]/LT[mon_ind])
						try:
							dndat[f-1][:]=dndat[f-1][:]+tmpdn
							dnerr[f-1][:]=sqrt(dnerr[f-1][:]**2+tmpdnERR**2)
						except:
							dndat[f-1][0:len(tmpdn)]=dndat[f-1][0:len(tmpdn)]+tmpdn[0:len(tmpdn)]
							dnerr[f-1][0:len(tmpdn)]=sqrt(dnerr[f-1][0:len(tmpdn)]**2+tmpdnERR[0:len(tmpdnERR)]**2)
						
						mon_ind=mon_ind+1

		#package and return
		self.updata=updat
		self.dndata=dndat
		self.uperr=uperr
		self.dnerr=dnerr
		self.mcpraw=updat-dndat
		self.compton_up=compton_up
		self.compton_dn=compton_dn
		self.diode_up=monitor_up
		self.diode_dn=monitor_dn
		self.fluo_up=fluo_up
		self.fluo_dn=fluo_dn
		self.detectors=detectors
		for f in (detectors):
			ch[f-1][:]=np.arange(1,self.mca_size+1)
			xscale[f-1][:]=np.arange(1,self.mca_size+1)
		self.xscale=xscale
		self.xlabel='Channel'
		self.ch=ch

	def loadspring8data(self,*args,**kwargs):
		#out_dat=mcp()
		#print type(out_dat)
		#out_dat.numDets=13
		#out_dat.mca_size=4096
		#out_dat.prefix=prefix
		#num_dets=13
		#mca_size=4096
		perscan=0
		
		perscan=1
		
		updat=np.zeros([self.num_dets,self.mca_size])
		dndat=np.zeros([self.num_dets,self.mca_size])
		uperr=np.zeros([self.num_dets,self.mca_size])
		dnerr=np.zeros([self.num_dets,self.mca_size])
		ch=np.zeros([self.num_dets,self.mca_size])
		xscale=np.zeros([self.num_dets,self.mca_size])
		upmag=[1,4,6,7]
		dnmag=[2,3,5,8]
		upmon=0.0
		dnmon=0.0
		updata=np.zeros([self.mca_size,self.num_dets])
		dndata=np.zeros([self.mca_size,self.num_dets])
		updataerr=np.zeros([self.mca_size,self.num_dets])
		dndataerr=np.zeros([self.mca_size,self.num_dets])
		tmpdatnorm=np.zeros([self.mca_size,self.num_dets])
		tmpdnornERR=np.zeros([self.mca_size,self.num_dets])
		ch=numpy.zeros((self.mca_size,self.num_dets))
		xscale=np.zeros((self.num_dets,self.mca_size))
		ch=numpy.zeros((self.num_dets,self.mca_size))
		tmpvec=np.zeros((10))
		tmpdat=np.zeros((8192,10))
		usecols=[2,3,4,5,6,7,8,9,10,11]
		updatfilesuffix=[1,4,6,7]
		dndatfilesuffix=[2,3,5,8]
		for j  in range(self.scans):
			print 'reading scan ',j+1
			for i in updatfilesuffix:
				fnameIn=self.prefix+'_'+str.rjust(str(j+1),3,'0')+'_'+str(i)+'_A.dat'
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
					jj=0
					for i in usecols:
						tmpvec[jj]=float(tmpstr[i])
						jj=jj+1
					tmpdat[lline,:]=tmpvec
					lline=lline+1
						
				for det in (self.detectors):
					dd=det-1
					tmpdatnorm[:,dd]=tmpdat[:,dd]/(float(monline.split()[0])* (tmpdat[1,dd]/tmpdat[0,dd]))
					tmpdnornERR[:,dd]=sqrt(tmpdat[:,dd])/(float(monline.split()[0])* (tmpdat[1,dd]/tmpdat[0,dd]))
					updata[:,dd]=updata[:,dd]+((tmpdatnorm[:,dd]*tmpdnornERR[:,dd])/tmpdnornERR[:,dd])
					updataerr[:,dd] = sqrt(updataerr[:,dd]**2+tmpdnornERR[:,dd]**2)
					
			for i in dndatfilesuffix:
				fnameIn=self.prefix+'_'+str.rjust(str(j+1),3,'0')+'_'+str(i)+'_B.dat'
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
					jj=0
					for i in usecols:
						tmpvec[jj]=float(tmpstr[i])
						jj=jj+1
					tmpdat[lline,:]=tmpvec
					lline=lline+1
				for det in (self.detectors):
					dd=det-1
					tmpdatnorm[:,dd]=tmpdat[:,dd]/(float(monline.split()[0])*(tmpdat[1,dd]/tmpdat[0,dd]))
					tmpdnornERR[:,dd]=sqrt(tmpdat[:,dd])/(float(monline.split()[0])*(tmpdat[1,dd]/tmpdat[0,dd]))
					dndata[:,dd]=dndata[:,dd]+((tmpdatnorm[:,dd]*tmpdnornERR[:,dd])/tmpdnornERR[:,dd])
					dndataerr[:,dd]=sqrt(dndataerr[:,dd]**2+tmpdnornERR[:,dd]**2)
	
		updat=updata.transpose()
		dndat=dndata.transpose()
		uperr=updataerr.transpose()
		dnerr=dndataerr.transpose()
		#package and return
		self.updata=updat
		self.dndata=dndat
		self.uperr=uperr
		self.dnerr=dnerr
		self.mcpraw=updat-dndat
		self.compton_up=0
		self.compton_dn=0
		self.diode_up=upmon
		self.diode_dn=dnmon
		self.fluo_up=0
		self.fluo_dn=0
		self.detectors=self.detectors
		for f in (self.detectors):
			ch[f-1,]=np.arange(1,self.mca_size+1)
			xscale[f-1,:]=np.arange(1,self.mca_size+1)
		self.xscale=xscale
		self.xlabel='Channel'
		self.ch=ch

			
			
			








				