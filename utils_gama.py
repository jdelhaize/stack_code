from numpy import *
from pylab import *
from scipy import integrate
from scipy import optimize
from scipy.interpolate import interp1d
from constants import *
import scipy.special
from utils import *

def SolarMag_gama(band):
	sm = {'u':6.38, 'g':5.15, 'r':4.71, 'i':4.56, 'z':4.54}
	return sm[band]

def gama_lumfn_params(band):
	schecd = {'u':{'Mstar':-18.60,'phi':2.03, 'alpha':-1.03},\
'g':{'Mstar':-20.09,'phi':1.47, 'alpha':-1.10},\
'r':{'Mstar':-20.86,'phi':1.24, 'alpha':-1.12},\
'i':{'Mstar':-21.30,'phi':1.00, 'alpha':-1.17},\
'z':{'Mstar':-21.52,'phi':1.02, 'alpha':-1.14}}

	phi = schecd[band]['phi']  * 10**(-2)
	Mstar = schecd[band]['Mstar'] - 5*log10(h)
	alpha = schecd[band]['alpha']
	
	return phi, alpha, Mstar

def gama_lumfn_params_error(band):
	schecd = {'u':{'Mstar':0.03,'phi':0.065, 'alpha':0.01},\
'g':{'Mstar':0.03,'phi':0.04, 'alpha':0.01},\
'r':{'Mstar':0.03,'phi':0.035, 'alpha':0.01},\
'i':{'Mstar':0.03,'phi':0.035, 'alpha':0.01},\
'z':{'Mstar':0.035,'phi':0.03, 'alpha':0.01}}

	phierr = schecd[band]['phi']  * 10**(-2)
	Mstarerr = schecd[band]['Mstar'] - 5*log10(h)
	alphaerr = schecd[band]['alpha']
	
	return phierr, alphaerr, Mstarerr
	
def gama_smfn_params():
	phi1 = 3.96E-3 * (h/0.7)**3
	phi2 = 0.79E-3 * (h/0.7)**3
	alpha1 = -0.35
	alpha2 = -1.47
	xstar = 10**10.66 * (h/0.7)**(-2)
	return phi1, phi2, alpha1, alpha2, xstar

def lumfn_error_gama(band, L):

	phistar, alpha, Mstar = gama_lumfn_params(band)
	Lstar = absmag_to_lum( Mstar, SolarMag_gama(band) )
	phistarerr, alphaerr, Mstarerr = gama_lumfn_params(band)
	Lstarerr = absmag_to_lum( Mstarerr, SolarMag_gama(band) ) # Not sure about this!!!
	
	const=1.0
	
	dphi_dphistar = const * (L/Lstar)**(alpha+1) * e**(-1*L/Lstar)
	dphi_dalpha = const * phistar * log(L/Lstar) * (L/Lstar)**(alpha+1) * e**(-1*L/Lstar)
	dphi_dLstar = const * phistar * L**(alpha+2)/Lstar**(alpha+3) * e**(-1*L/Lstar) \
- const * phistar * (alpha+1) * L**(alpha+1)/Lstar**(alpha+2) *	e**(-1*L/Lstar)

	phierr = ( (dphi_dphistar * phistarerr)**2 + (dphi_dalpha * alphaerr)**2 + (dphi_dLstar * Lstarerr)**2 )**0.5 
	
	return phierr


def lumfn_gama_L(band, L):
	phi, alpha, Mstar = gama_lumfn_params(band)
	Lstar = absmag_to_lum( Mstar, SolarMag_gama(band) )
	y = schecter( phi, alpha, Lstar, L)
	return y

def lumfn_gama(band, M):
	phi, alpha, Mstar = gama_lumfn_params(band)
	y = schecter_absmag( phi, alpha, Mstar, M) / (0.4*log(10))
	return y
	
def lumdensfn_gama(band, M):
	phi, alpha, Mstar = gama_lumfn_params(band)
	y = 10**(0.4*(SolarMag_gama(band)-M)) * lumfn_gama( band, M)
	return y

def lumdens_gama(band, lolim, hilim):
	rho = integrate.quad( lambda x: lumdensfn_gama(band, x), lolim, hilim)[0]
	print 'rho(L) = ', '%.3e' %rho
	return rho
	
def lumdens_gama_gammaform(band):
	phi, alpha, Mstar = gama_lumfn_params(band)
	integ = phi * 10**(0.4*(SolarMag_gama(band)-Mstar)) * scipy.special.gamma(alpha+2)
	return integ

def lumdens_gama_gammaform_bounded(band, xmin, xmax):
	phi, alpha, Mstar = gama_lumfn_params(band)
	integ = -1 * phi * 10**(0.4*(SolarMag_gama(band)-Mstar)) * \
(scipy.special.gamma(alpha+2) * scipy.special.gammainc(alpha+2, absmag_to_lum(xmax,Mstar)) \
- scipy.special.gamma(alpha+2) *scipy.special.gammainc(alpha+2,absmag_to_lum(xmin,Mstar)))

	return integ
	
def plot_lumfn_absmag(band, xmin, xmax):
	x=linspace(xmin,xmax,num=100)
	figure()
	plot(x,lumfn_gama(band, x))
	yscale('log')
	ylabel(r'$\phi_{L}(L)\ \rm{Mpc}^{-3}\,\rm{0.5\ mag}^{-1}$')
	xlabel(r'$M$')
	
def plot_lumdensfn_absmag(band, xmin, xmax):
	x=linspace(xmin,xmax,num=100)
	figure()
	plot(x,lumdensfn_gama(band, x))
	yscale('log')
	ylabel(r'$\rho_{L}(L)\ L_{\odot}\,\rm{Mpc}^{-3}\,\rm{0.5\ mag}^{-1}$')
	xlabel(r'$M$')
	
	
def SMfn_gama(logx):
	phi1, phi2, alpha1, alpha2, xstar = gama_smfn_params()
	y = dbl_schecter_logx(phi1, phi2, alpha1, alpha2, xstar, logx)
	return y
	
def SMdensfn_gama(logx):
	y = 10**logx * SMfn_gama(logx)
	return y
	
	
def plot_SMdensfn(xmin, xmax):
	logx=linspace(xmin,xmax,num=100)
	figure()
	plot(logx,SMdensfn_gama(logx))
	yscale('log')
	ylabel(r'$\rho_{M*}(M*)\ M_{\odot}\,\rm{Mpc}^{-3}$')
	xlabel(r'$\log(M*)$')
	
# Now define the SM density integral in terms of the gamma function.

def SMdensfn_gama_gammaform_bounded(xmin, xmax):
	phi1, phi2, alpha1, alpha2, xstar = gama_smfn_params()
	integ = dbl_schecter_integral(phi1, phi2, alpha1, alpha2, xstar, xmin, xmax)
	return integ
	
def SMdensfn_gama_gammaform():
	phi1, phi2, alpha1, alpha2, xstar = gama_smfn_params()
	integ_part = lambda phi,alpha: phi*xstar* (scipy.special.gamma(alpha+2))
	integ = integ_part(phi1,alpha1) + integ_part(phi2,alpha2)
	return integ

def SMdens_gammaform_params(phi1, phi2, alpha1, alpha2, xstar):
	integ_part = lambda phi,alpha: phi*xstar* (scipy.special.gamma(alpha+2))
	integ = integ_part(phi1,alpha1) + integ_part(phi2,alpha2)
	return integ

def completenessF_gama_r(zmax, appmag_lim, schec_lo, SolarMag):
	#appmag_lim = 19.4 # Limiting magnitude of the survey
	#schec_lo=-24
	#SolarMag_r = 4.71
	
	Lmin = z_to_L( zmax, appmag_lim, SolarMag, ary='y')
	Mmax = lum_to_absmag(Lmin, SolarMag)

	Lfull = lumdens_gama_gammaform('r')
	Lpart = lumdens_gama_gammaform_bounded('r', schec_lo, Mmax)
	Lpercen = Lpart/Lfull
	LcompleteF = 1.0 / Lpercen

	return LcompleteF

# Based on a manually-defined SM 'flux' survey limit (gmin)
def SMcompletenessF_gama_JD(zmax, gmin, schec_hi):

	#gmin = 3.63 # Where g = SM/DL**2
	#schec_hi = 10**14
	
	DLmax = z_to_DL(zmax)
	SMmin = 10**gmin * DLmax**2

	Mfull = SMdensfn_gama_gammaform()
	Mpart = SMdensfn_gama_gammaform_bounded(SMmin, schec_hi)
	Mpercen = Mpart/Mfull
	McompleteF = 1.0 / Mpercen
	
	return McompleteF, SMmin
	

# Calculate the SM completeness limit at z, defined by Aaron Robotham's non colour-biased cut
def SMcompleteness_AR(z):
	ar=genfromtxt('../catalogues/zmax_19p4_afo_logmstar.dat') #'/Users/jacinta/workdir/stacking/catalogues/zmax_19p4_afo_logmstar.dat'
	zar = ar[:,1]
	smar = ar[:,0]
	smar100 =  10**smar *(h/0.7)**(-2) 
	arfit = interp1d(zar,smar100)
	SMlim = arfit(z)
	return SMlim
	
	
def SMcompletenessF_gama(z, schec_hi):
	SMmin = SMcompleteness_AR(z)
	Mfull = SMdensfn_gama_gammaform()
	Mpart = SMdensfn_gama_gammaform_bounded(SMmin, schec_hi)
	Mpercen = Mpart/Mfull
	McompleteF = 1.0 / Mpercen
	return McompleteF, SMmin






	