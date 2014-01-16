#!/usr/bin/env python

from numpy import *
from pylab import *
from scipy import integrate
from scipy import optimize
from constants import *
import scipy.special

# Subtracts a baseline from the spectrum, fitting within the range low-hi. DOES exclude signal region by default, but this can be turned off. Works on both 1d and multi-dimensional y-matrices.

# Defines a channel mask based on the signal region and the min and max channels of interest
# Note: If want full spectrum (except signal region) use xlolim=min(x) etc...

# Define a mask over the signal region
def mask_sigreg(ary,siglolim,siguplim):
	msklo=where(ary<siglolim)[0]
	mskhi=where(ary>siguplim)[0]
	msk=concatenate([msklo,mskhi])
	return msk

def x_mask(x1,xlolim,xuplim,siglolim,siguplim,excl_sig='y',multid_x='n'):

	if multid_x=='n':
		x=ravel(x1)
	else:
		x=x1
	
	if excl_sig=='y':	
		msklo=where(logical_and(x>xlolim,x<siglolim))[0]
		mskhi=where(logical_and(x>siguplim,x<xuplim))[0]
		xmsk=concatenate([msklo,mskhi])
	elif excl_sig=='n':
		xmsk=where(logical_and(x>xlolim,x<xuplim))[0]
		
	return xmsk	
	
#4th order
# NOTE: The multi_x='y' doesn't work
def basefit(x1,y,xlolim,xuplim,siglolim,siguplim,excl_sig='y',multid_y='n',plt='n',multid_x='n'):
# 	for key,value in d.iteritems():
# 	exec(str(key)+"="+str(value))	
		
	msk=x_mask(x1,xlolim,xuplim,siglolim,siguplim,excl_sig,multid_x=multid_x)
	if multid_x=='n':
		x=ravel(x1)
	else:
		x=x1
	
	if multid_y=='n':
		a,b,c,d,f=ma.polyfit(x[msk],y[msk],4)
	
	# Make sure new fit has same dimensions as y if y has more than 1 dimension
	elif multid_y=='y':
	
		a=[]
		b=[]
		c=[]
		d=[]
		f=[]
		num=shape(y)[1]
		for i in range(num):
			zi=ma.polyfit(x[msk],y[msk,i],4)
			a=a+[zi[0]]
			b=b+[zi[1]]
			c=c+[zi[2]]
			d=d+[zi[3]]
			f=f+[zi[4]]
		
		x=x1
		
	newy=y-(a*x**4+b*x**3+c*x**2+d*x+f)
	
	if plt=='y':
		figure()
		plot(x,y,'k')
		plot(x,a*x**4+b*x**3+c*x**2+d*x+f,'r')
	
	return newy,msk

# Third order
# def basefit(x1,y,xlolim,xuplim,siglolim,siguplim,excl_sig='y',multid_y='n',plt='n'):
# # 	for key,value in d.iteritems():
# # 	exec(str(key)+"="+str(value))
# 		
# 	msk=x_mask(x1,xlolim,xuplim,siglolim,siguplim,excl_sig='y')
# 	x=ravel(x1)
# 	
# 	if multid_y=='n':
# 		a,b,c,d=ma.polyfit(x[msk],y[msk],3)
# 	
# 	# Make sure new fit has same dimensions as y if y has more than 1 dimension
# 	elif multid_y=='y':
# 	
# 		a=[]
# 		b=[]
# 		c=[]
# 		d=[]
# 		num=shape(y)[1]
# 		for i in range(num):
# 			zi=ma.polyfit(x[msk],y[msk,i],3)
# 			a=a+[zi[0]]
# 			b=b+[zi[1]]
# 			c=c+[zi[2]]
# 			d=d+[zi[3]]
# 		
# 		x=x1
# 		
# 	newy=y-(a*x**3+b*x**2+c*x+d)
# 	
# 	if plt=='y':
# 		figure()
# 		plot(x,y,'k')
# 		plot(x,a*x**3+b*x**2+c*x+d,'r')
# 	
# 	return newy,msk

	
# Fit cubic polynomial to whole baseline (using masked arrays)
def cubicfit(x,y):
	p=ma.polyfit(x,y,4)
	newy=y-polyval(p,x)
	return newy

# Plots a flux spectrum
def stack_var(ary,usecol,colour,lb):
	dat=genfromtxt(ary,usecols=usecol)
	x=dat[:,0]
	y=dat[:,1]*1000
	plot(x,y,colour,label=lb)
	return x,y
	
# Integrate a spectrum between the signal region
def integrate_spectrum(ary,sigreg,spec_res):
	integ = ma.sum(ary[sigreg],axis=0)*spec_res
	return integ

# Plots vertical lines at edges of signal region	
def siglims_plot(xlo,xhi,Colour='m',Line='--',lw=1):
	ymin=min(ylim())
	ymax=max(ylim())
	vlines(xlo,ymin,ymax,color=Colour,linestyle=Line,lw=lw)
	vlines(xhi,ymin,ymax,color=Colour,linestyle=Line,lw=lw)

# Print value in scientific notation. Type should be 'g' (single value) or 's' (array)
def print_sci(val,type):
	if type=='2df':
		out='%.3e' %val
	if type=='g':
		out=map(lambda x: '%.3e' %x,val)
	return out
	
# Same source code as masked_inside() but with the 'all' added.
def masked_between(x, v1, v2, copy=True):
   if all(v2 < v1):
	   (v1, v2) = (v2, v1)
   xf = ma.filled(x)
   condition = (xf >= v1) & (xf <= v2)
   return ma.masked_where(condition, x, copy=copy)
   
# My function to blank BETWEEN bad channels
# 	def mask_between_me(ary,bad_lo,bad_hi):
# 		ary_ma1=ma.masked_less(ary,bad_hi)
# 		ary_ma2=ma.masked_greater(ary,bad_lo)
# 		ary_ma3=ma.masked_where(ary_ma1.mask==ary_ma2.mask,ary)
# 		return ary_ma3
   
# Find the angular separation between two sources. Coordinates in radians
def angular_sep(ra1,ra2,dec1,dec2):
# 	if ra1==ra2 and dec1==dec2:
# 	   ang_sep_rad=0.0
# 	else:
	ang_sep_rad=arccos(cos(dec1)*cos(ra1)*cos(dec2)*cos(ra2)+cos(dec1)*sin(ra1)*cos(dec2)*sin(ra2)+sin(dec1)*sin(dec2))
	
	ang_sep_deg=ang_sep_rad*180/pi
	return ang_sep_deg
	   
# Convolve a value with a gaussian. Outputs scale factor based on beam response away from centre.
def convolve_gaussian(angsep,beamwidth):
	r0=beamwidth/2.355
	weight=e**(-angsep**2/(2*r0**2))
	return weight
	
#***	
	
# Error on the luminosity given the redshift error and the absolute magnitude

def Merr_to_Lerr(Merr, L):
	lum_err = 0.4*log(10)*Merr*L
	return lum_err

def lum_error(z,zerr, L, m_err):

	dl_fracerr_2 = (zerr/z)**2 + (zerr/(1+z))**2	
	Merr = ( m_err**2 + (5/log(10))**2 * dl_fracerr_2 )**0.5	
	#lum_err = 0.4*Merr*L/M
	lum_err = Merr_to_Lerr(Merr, L)
	return lum_err
	
# Error on the frequency, given the error on the redshift	
def freq_error(f,z, zerr):
	freq_err = (f*zerr)/(1+z)
	return freq_err

# Given a set of values with errors, find the error on the average.	Note: axis=0 means it can work over a ugriz array
def error_on_av(err):
	er_on_av = (ma.sum(err**2,axis=0)/ma.size(err,axis=0)**2)**0.5
	return er_on_av
	   
# Test if file exists. 'name' should be a string. Outputs True or False.
def existence(name):	   
   import os.path
   return os.path.exists(name)
   
# Convert ra from hexadecimal to degrees
def ra_hextodeg(h,m,s):
	ra=(h+m/60.+s/3600.)*15.
	return ra
	
# Convert dec from hexadecimal to degrees
def dec_hextodeg(d,m,s):
	dst=str(d)
	if dst[0]=='-':
		dec=float(d)-float(m)/60.-float(s)/3600.
	elif float(d)>=0.0:
		dec=float(d)+float(m)/60.+float(s)/3600.
	elif float(d)<0.0:
		dec=float(d)-float(m)/60.-float(s)/3600.
	return dec
# Caution!!! Check d=-00 case!

# Make an annotation file from a list of RA and Dec values	   
def make_ann_file(ra,dec):
	size=['0.05']*len(ra)
	shape=['CROSS']*len(ra)
	string=vstack((shape,ra,dec,size,size)).T
	return string

# Convert redshift to frequency	
def z_to_freq(z,f0):
	f=f0/(1.0+z)
	return f
	
# Convert frequency to redshift
def freq_to_z(f,f0):
	z=(f0/f)-1
	return z

# Convert redshift to NON-RELATIVISTIC velocity
def z_to_vel(z):
	v=c*z
	return v
	
def vel_to_z(v):
	z=v/c
	return z	

# Convert NON-RELATIVISTIC velocity to frequency
def vel_to_freq(vel,f0):
	z=vel_to_z(vel)
	f=z_to_freq(z,f0)
	return f
	
# Convert NON-RELATIVISTIC frequency to velocity
def freq_to_vel(f,f0):
	z=freq_to_z(f,f0)
	v=z_to_vel(z)
	return v
	
def delfreq_to_delvel(delf, f0):
	delv = delf * c / f0
	return delv

def delvel_to_delfreq(delv, f0):
	delf = delv * f0 / c
	return delf

	
# Calculate luminosity distance from redshift and cosmology parameters
# NOTE: Uses cosmology parameters from constants.py
def z_to_DL(z):
	E=lambda z: sqrt(Omega_M*(1.+z)**3.+Omega_k*(1.+z)**2.+Omega_lambda)
	integral=[]
	# Note: This is for an array of z values
	for m in z:
	  ary=integrate.quad(lambda x: 1./E(x),0.,m)
	  integral.append(ary[0])
	  
	DC=DH*array(integral)
	if Omega_k==0:
	  DM=DC
	else:
	  print 'Omega_k does not equal zero'
	 
	DA=DM/(1.+z)
	DL=(1+z)**2*DA
	DL=DL/(MpctoMet) #convert to megaparsec
	return DL

def z_to_DC_ary(z):
	E=lambda z: sqrt(Omega_M*(1.+z)**3.+Omega_k*(1.+z)**2.+Omega_lambda)
	integral=[]
	# Note: This is for an array of z values
	for m in z:
	  ary=integrate.quad(lambda x: 1./E(x),0.,m)
	  integral.append(ary[0])
	  
	DC=DH*array(integral)
	DC=DC/(MpctoMet) #convert to megaparsec
	return DC

def z_to_DL_nonary(z):
	#Note: This is for a single z
	E=lambda z: sqrt(Omega_M*(1.+z)**3.+Omega_k*(1.+z)**2.+Omega_lambda)
	integral=integrate.quad(lambda x: 1./E(x),0.,z)[0]
	  
	DC=DH*array(integral)
	if Omega_k==0:
	  DM=DC
	else:
	  print 'Omega_k does not equal zero'
	 
	DA=DM/(1.+z)
	DL=(1+z)**2*DA
	DL=DL/(MpctoMet) #convert to megaparsec
	return DL
	
def z_to_DA_Mpc(z):
	#Note: This is for a single z
	E=lambda z: sqrt(Omega_M*(1.+z)**3.+Omega_k*(1.+z)**2.+Omega_lambda)
	integral=integrate.quad(lambda x: 1./E(x),0.,z)[0]
	  
	DC=DH*array(integral)
	if Omega_k==0:
	  DM=DC
	else:
	  print 'Omega_k does not equal zero'
	 
	DA=DM/(1.+z)
	DA=DA/(MpctoMet) #convert to megaparsec
	return DA
	
def z_to_DC_Mpc(z):
	#Note: This is for a single z
	E=lambda z: sqrt(Omega_M*(1.+z)**3.+Omega_k*(1.+z)**2.+Omega_lambda)
	integral=integrate.quad(lambda x: 1./E(x),0.,z)[0]
	  
	DC=DH*array(integral)
	DC=DC/(MpctoMet) #convert to megaparsec
	return DC

# Separately define the function E(z) as defined in Hogg (2000)
def Ez(z):
	E=sqrt(Omega_M*(1.+z)**3.+Omega_k*(1.+z)**2.+Omega_lambda)
	return E

# The evolving Hubble constant	
def Hz(z):
	Hz=H0*Ez(z)
	return Hz

# Convert degrees on the sky to Mpc
def deg_to_Mpc(deg,z):
	mpc=z_to_DC_Mpc(z)*deg*pi/180.0
	return mpc

# Convert degrees on the sky to Mpc (array)
def deg_to_Mpc_ary(deg,z):
	mpc=z_to_DC_ary(z)*deg*pi/180.0
	return mpc

# Find the sky volume from an ra and dec range (in degrees), and the maximum redshift
def sky_vol(ra,dec,z):
	V=z_to_DC_Mpc(z)*deg_to_Mpc(ra,z)*deg_to_Mpc(dec,z)/3
	return V

# Find the sky volume from an ra and dec range (in degrees), and an array of redshift
def sky_vol_ary(ra,dec,z):
	V=z_to_DC_ary(z)*deg_to_Mpc_ary(ra,z)*deg_to_Mpc_ary(dec,z)/3
	return V

# Find the weighted sky volume
def sky_vol_weight(ra,dec,z,w):
	V=sum(w*z_to_DC_ary(z)*deg_to_Mpc_ary(ra,z)*deg_to_Mpc_ary(dec,z)/3)/sum(w)
	return V
	
# Calculate rhoHI from a mass-to-light ratio. 	
def rhoHI_calc(MassLumRatio,rhoL):
	rhoHI=MassLumRatio*rhoL
	return rhoHI
   
# Convert rhoHI to OmegaHI. Note constants defined.
#def OmegaHI_calc(rhoHI):
def rho_to_OmegaHI(rhoHI,z=0):	
	Omega_HI=(8*pi*G*rhoHI)/(3*Hz(z)**2)
	return Omega_HI
	
def OmegaHI_to_rho(OmegaHI,z=0):	
	rhoHI=OmegaHI / ((8*pi*G)/(3*Hz(z)**2))
	return rhoHI
	
# Calculate error on rhoHI and OmegaHI
def OHI_error(OHI,MHI,MHI_err,rhoL,rhoL_err,Lav,z=0,Lav_err=0):

	rhoHI_err = ( (rhoL*MHI_err/Lav)**2+(MHI*rhoL_err/Lav)**2+(-1*MHI*rhoL*Lav_err/Lav**2)**2 )**0.5
	OHI_err = (8*pi*G*rhoHI_err)/(3*Hz(z)**2)
	
	return OHI_err, rhoHI_err	
	
# Calculate error on rhoHI and OmegaHI when OHI = <M/L> x rho_L
def OHI_error_v2(ML,ML_err,rhoL,rhoL_err,z=0):

	rhoHI_err = ( (rhoL*ML_err)**2 + (ML*rhoL_err)**2 )**0.5
	OHI_err = (8*pi*G*rhoHI_err)/(3*Hz(z)**2)
	
	return OHI_err, rhoHI_err	
	
# Convert apparent magnitude to absolute magnitude	
def app_to_absmag(DL,appmag):
	absmag=appmag-5.*log10(DL*10.**5.)
	return absmag

# Convert absolute magnitude to luminosity (careful of AB vs Vega)
def absmag_to_lum(absmag,SolarMag):
	lum=10.**((SolarMag-absmag)*0.4)
	return lum
	
def lum_to_absmag(lum,SolarMag):
	absmag = SolarMag - 2.5*log10(lum)
	return absmag

def z_to_L( z, appmag, SolarMag, ary='y'):
	if ary=='y':
		dl = z_to_DL(z)
	if ary=='n':
		dl = z_to_DL_nonary(z)
	absmag = app_to_absmag(dl, appmag)
	l = absmag_to_lum(absmag, SolarMag)
	return l
	
# Save a spectrum into a text file
def save_spect(x,y,file):
	a=y
	a.shape=(len(a),1)
	b=concatenate((x,a),axis=1)
	savetxt(file,b)

# Find the nearest value in an array
def find_nearest(array,value):
	idx=(abs(array-value)).argmin()
	return idx,array[idx]
	
def find_level(perc,x,y,peak,xpeak):
	lev=peak*perc/100
	idx_closest,val_closest=find_nearest(y,lev)
	return x[idx_closest],val_closest
	
def find_sigrange(x,y,perc):
	peak,xpeak=max(y),x[where(y==max(y))[0]]
	msklo=where(x<=freq0)[0]
	mskhi=where(x>freq0)[0]
	xlo,ylo=x[msklo],y[msklo]
	xhi,yhi=x[mskhi],y[mskhi]
	xlimlo,vallo=find_level(perc,xlo,ylo,peak,xpeak)
	xlimhi,valhi=find_level(perc,xhi,yhi,peak,xpeak)
	return {'xlo':xlimlo,'xhi':xlimhi,'ylo':vallo,'yhi':valhi}

# Find the 1D fft of a spectrum and plot the magnitude
def fft_spec(ary):
	ft_spec=fft.fft(ary)
	mag=abs(ft_spec)
	figure()
	plot(mag)
	ift_spec=fft.ifft(ft_spec)
	return ft_spec,mag,ift_spec
	
# Convert a string to latex regular text string
def texstr(str):
	outstr=r'$\rm{'+str+'}$'
	return outstr

# Fit a function to a spectrum. See OmegaHI_sin_fit on how to define fitfunc	
def fit_function(sx,sy,p0,fitfunc):

	# Difference between target function and actual data
	errfunc = lambda p, sx, sy: (fitfunc(p, sx) - sy)
	# Minimising the difference (least squares) to find the coefficients
	p1, success = optimize.leastsq(errfunc, p0[:], args=(sx, sy))
	#print success
	# The resulting fit
	sy_fit=fitfunc(p1,sx)
	
	return p1,sy_fit
	
# Fit a function to a spectrum, weighting by the y-errors
def fit_function_yweight(sx, sy, sz, p0, fitfunc):
	errfunc = lambda p, sx, sy, sz: (fitfunc(p, sx) - sy) /sz
	p1, success = optimize.leastsq(errfunc, p0[:], args=(sx, sy, sz))
	sy_fit=fitfunc(p1,sx)	
	return p1,sy_fit

# Fit a straight line, accounting for y-errors
def fit_line(logx,logy,logyerr,pinit): # Note: values don't have to be logs

	fitfunc = lambda p, x: p[0] + p[1] * x
	errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
	
	#pinit = [10, -.1]
	out = optimize.leastsq(errfunc, pinit, args=(logx, logy, logyerr), full_output=1)
	
	pfinal = out[0]
	yint, slope = pfinal[0], pfinal[1]
	covar = out[1]
	
	return yint, slope, covar
	
def lumfn_2dfgrs(L, SolarMag=5.25):

	
	phi = 1.61*10**-2 * h**3
	alpha = -1.21
	Mstar = -19.66 + 5*log10(h)
	Lstar = absmag_to_lum( Mstar, SolarMag )
	
	N = 0.921 * phi * (L/Lstar)**(alpha+1) * e**(-1*L/Lstar)
	#N = 0.921 * phi * (L/Lstar)**(alpha) * e**(-1*L/Lstar) # test version only
	
	return N

def lumfn_error_2dfgrs(L, SolarMag=5.25):

 	phistar = 1.61*10**-2 * h**3
	alpha = -1.21
	Mstar = -19.66 + 5*log10(h)
	Lstar = absmag_to_lum( Mstar, SolarMag )
	
	phistarerr = 0.08*10**-2 * h**3
	alphaerr = 0.03
	Mstarerr = 0.07
	Lstarerr = ( (-0.4 * 10**(0.4*SolarMag * 5*log10(h)) * log(10) * 10**(-0.4*Mstar) * Mstarerr)**2 )**0.5
	
	const = 0.4 * log(10)
	
	dphi_dphistar = const * (L/Lstar)**(alpha+1) * e**(-1*L/Lstar)
	dphi_dalpha = const * phistar * log(L/Lstar) * (L/Lstar)**(alpha+1) * e**(-1*L/Lstar)
	dphi_dLstar = const * phistar * L**(alpha+2)/Lstar**(alpha+3) * e**(-1*L/Lstar) \
- const * phistar * (alpha+1) * L**(alpha+1)/Lstar**(alpha+2) *	e**(-1*L/Lstar)

	phierr = ( (dphi_dphistar * phistarerr)**2 + (dphi_dalpha * alphaerr)**2 + (dphi_dLstar * Lstarerr)**2 )**0.5 
	
	return phierr

def schecter( phi, alpha, xstar, x ):
	y = phi * (x/xstar)**(alpha+1) * e**(-1*x/xstar)
	#y = (phi/xstar) * (x/xstar)**(alpha) * e**(-1*x/xstar) # Test only!
	return y
	
def schecter_logx( phi, alpha, xstar, x ):
	y = phi * log(10) * (x/xstar)**(alpha+1) * e**(-1*x/xstar)
	return y
	
def schecter_logxy( phi, alpha, xstar, x ):
	y = log10(phi*log(10)) +(alpha+1)* log10(x)-log10(xstar) - log10(e) * 10**(log10(x)-log10(xstar))
	return y
	
def dbl_schecter(phi1, phi2, alpha1, alpha2, xstar, x):
	y = log(10) * (phi1 * (x/xstar)**(alpha1+1) + phi2 * (x/xstar)**(alpha2+1)) * e**(-1*x/xstar)
	return y
	
def dbl_schecter_logx(phi1, phi2, alpha1, alpha2, xstar, logx):
	y = log(10) * ( phi1 * (10**logx/xstar)**(alpha1+1) + phi2 * (10**logx/xstar)**(alpha2+1)) * e**(-1*10**logx/xstar)
	return y
	
def dbl_schecter_dens_logx(phi1, phi2, alpha1, alpha2, xstar, logx):
	y = xstar*log(10) * ( phi1 * (10**logx/xstar)**(alpha1+2) + phi2 * (10**logx/xstar)**(alpha2+2)) * e**(-1*10**logx/xstar)
	return y

def schecter_fn_to_integrate(phi, alpha, xstar, x):
	y = phi * log(10) * (x/xstar)**(alpha+1) * e**(-1*x/xstar)
	return y
	
def schecter_absmag( phi, alpha, Mstar, M):
	y = 0.4 * log(10) * phi * (10**(0.4*(Mstar-M)))**(alpha+1) * e**(-1*10**(0.4*(Mstar-M)))
	return y

def schecter_absmag_logy(phi, alpha, Mstar,M):
	y = log10(0.4*log(10)) + log10(phi) + (alpha+1)*0.4*(Mstar-M) - 10**(0.4*(Mstar-M)) * log10(e)
	return y

def schecter_integral(phistar, alpha, Lstar, xmin, xmax):
	import scipy.special
	integ = phistar*Lstar* (scipy.special.gamma(alpha+2) * scipy.special.gammainc(alpha+2, xmax/Lstar) - scipy.special.gamma(alpha+2) *scipy.special.gammainc(alpha+2, xmin/Lstar))
	return integ
	
# def dbl_schecter_integral(phi1, phi2, alpha1, alpha2, xstar, xmin, xmax):
# 	integ_part = lambda phi,alpha: phi*xstar**2* (scipy.special.gamma(alpha+1) * scipy.special.gammainc(alpha+1, xmax/xstar) - scipy.special.gamma(alpha+1) *scipy.special.gammainc(alpha+1, xmin/xstar))
# 	integ = integ_part(phi1,alpha1) + integ_part(phi2,alpha2)
# 	return integ

def dbl_schecter_integral(phi1, phi2, alpha1, alpha2, xstar, xmin, xmax):
	integ_part = lambda phi,alpha: phi*xstar* (scipy.special.gamma(alpha+2) * scipy.special.gammainc(alpha+2, xmax/xstar) - scipy.special.gamma(alpha+2) *scipy.special.gammainc(alpha+2, xmin/xstar))
	integ = integ_part(phi1,alpha1) + integ_part(phi2,alpha2)
	return integ
	
# def dbl_schecter_integral(phi1, phi2, alpha1, alpha2, xstar, xmin, xmax):
# 	integ_part = lambda phi,alpha: log(10)*phi*xstar**2* (scipy.special.gamma(alpha+2) * scipy.special.gammainc(alpha+2, xmax/xstar) - scipy.special.gamma(alpha+2) *scipy.special.gammainc(alpha+2, xmin/xstar))
# 	integ = integ_part(phi1,alpha1) + integ_part(phi2,alpha2)
# 	return integ
	
def schecter_integral_absmag(phistar, alpha, Mstar, xmin, xmax):
	import scipy.special
	integ = -1*phistar* (scipy.special.gamma(alpha+1) * scipy.special.gammainc(alpha+1, absmag_to_lum(xmax,Mstar)) - scipy.special.gamma(alpha+1) *scipy.special.gammainc(alpha+1, absmag_to_lum(xmin,Mstar)))
	return integ
	
def schecter_integral_absmag_logy(phi, alpha, Mstar, xmin, xmax):
	integ = lambda M: M*log10(0.4*log(10)*phi) + 0.4*M*Mstar*(alpha+1) - (0.4*(alpha+1)*M**2)/2.0 + log10(e)*10**(0.4*(Mstar-M))/(0.4*log(10))
	a_to_b = integ(xmax) - integ(xmin)
	print 'integ(xmax), integ(xmin) = ',integ(xmax), integ(xmin)
	return a_to_b
	
def schecter_integral_logx(phistar, alpha, Lstar, xmin, xmax):
	import scipy.special
	integ = phistar * Lstar**2 * log(10) * (scipy.special.gamma(alpha+3) * scipy.special.gammainc(alpha+3, xmax/Lstar) - scipy.special.gamma(alpha+3) *scipy.special.gammainc(alpha+3, xmin/Lstar))
	return integ

# The definition of the bj luminosity function according to Norberg 2001	
def L_density_fn( L ):
	(phi, alpha, Lstar) = (1.61e-2, -1.21, 9.54993e9)
	y = 0.921 * phi * (L/Lstar)**(alpha+1) * e**(-1*L/Lstar) * L
	return y

# Same as above but in log-log space	
def L_density_loglog_fn( logl ):
	(phi, alpha, Lstar) = (1.61e-2, -1.21, 9.54993e9)
	logy = log10(0.921) + log10(phi) + (alpha+1) * log10( logl/Lstar ) \
- (10**logl/Lstar) * log10(e) + logl
	return logy
	
# Find the inclination of a galaxy from its semi major and minor axes
def find_incl(a, b, eos=0.12):

	#eos = 0.12 # Width of Edge-on-spiral. Defined in Meyer 08.
	cosi2= ( (b/a)**2 - eos ) / (1-eos)
	i = arccos( sqrt(cosi2) )
	
	bd=where(cosi2<0)[0] # These are the edge-on galaxies (are flatter than eos)
	i[bd] = pi/2  # So set their inclination to 90deg
	
	return i

# Use the Tully-Fisher relation to convert absolute magnitude to rotational velocity
def TF_vel( absmag, offset, slope ):
	v = 10**( (absmag-offset) /slope )
	return v

# Convert the rotational velocity to line-width:
def rotvel_to_W(rotvel, incl):
	W=2*rotvel*sin(incl)
	return W

# Convert between Supercosmos/APM bj and rf bands into B band (probably Johnson B).	
def bj_to_B(bj, rf):
	B = bj - 0.005 + 0.236 * (bj - rf)
	return B
	
# At restframe:
def dv_to_df(dv):
	df = dv * freq0 / c
	return df
	
# Test if there are any spectra that are *completely* blanked
def all_blank_test(y_ma):
	for i in range(len(y_ma[0])):
		numblnk = all(y_ma[:,i].mask)
		blnk = where(y_ma[:,i].mask==True)
		num = len(blnk[0])
		#print i, num
		if numblnk==True:
			print i

# Identify G15 sources
def which_gama(ra):		
	if ra>211:
		gtyp='g15'
	else:
		gtype='g9'
	return gtype
	
def deg_to_rad(deg):
	rad = deg*math.pi/180.
	return rad
	
def rad_to_deg(rad):
	deg = rad * 180. / math.pi
	return deg
	
def gama_absmag_fluxscale(M, fluxscale):
	M_fix = M - 2.5*log10(fluxscale) + 5*log10(h/0.7)
	return M_fix

def gama_stellarmass_fluxscale(logM, fluxscale):
	logM_fix = logM + log10(fluxscale) - 2*log10(h/0.7)
	return logM_fix

def restrict_fluxscale(fluxscale):
	fluxscale[where(fluxscale<1.0)] = 1.0
	fluxscale[where(fluxscale>1.3)] = 1.3
	return fluxscale
	
def lum_completeness_lim(z_lim, appmag_lim, SolarMag):
	L_lim = z_to_L( z_lim, appmag_lim, SolarMag)
	return L_lim
	
def sm_completeness_lim(DLmax, gmin):
	SMmin = 10**gmin * DLmax**2
	return SMmin
	
# Combine errors in quadrature
def quad_error(erary):
	comb=0
	for i in erary:
		comb+=i**2
	comb = (comb)**0.5
	return comb

#********************************************
#Find lookback time in Gyrs for proch,rao etc
def z_to_LBT_1d(z):
	#**Calculate the luminosity distance**
	E=lambda zi: sqrt(Omega_M*(1.+zi)**3.+Omega_k*(1.+zi)**2.+Omega_lambda)
	integral=[]
	for m in z:
		ary=integrate.quad(lambda x: 1./((1+x)*E(x)),0.,m)
		integral.append(ary[0])
	integral=array(integral)
	tl=integral*tH
	Gyr=tl/10**9 #convert from Mpc to Gly
	return Gyr	

#********************************************

#Find lookback time in Gyrs for zwaan, lah etc.
def z_to_LBT_2d(z):
	#**Calculate the luminosity distance**
	E=lambda zi: sqrt(Omega_M*(1.+zi)**3.+Omega_k*(1.+zi)**2.+Omega_lambda)
	integral=integrate.quad(lambda x: 1./((1+x)*E(x)),0.,z)
	tl=integral[0]*tH
	Gyr=tl/10**9 #convert from Mpc to Gly
	#print 'z =',z,'t=',Gyr,'Gyr'
	return Gyr
	
#********************************************