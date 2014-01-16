#!/usr/bin/env python

from numpy import *
from pylab import *
from scipy import *
from scipy import optimize
from OmegaHI_vars import *
from utils import *
#ion()

#########################################################################

# Define useful parameters to make the plot pretty
def graph_props(typ):
	hlines(0,min(xlim()),max(xlim()), colors='k', linestyles='dotted')
	vlines(freq0,min(ylim()),max(ylim()),linestyle='dotted',color='k')
	xlabel(texstr('Frequency\ (MHz)'))
	if typ=='f':
		ylabel(texstr('Flux\ (mJy)'))
	if typ=='m':
		ylabel(r'$\rm{Average\ Mass}\ (M_{\odot})$')
	legend(loc=0,frameon=False)
	xlim(xlolim,xuplim)

#########################################################################

# Function to fit and subtract a sinusoid from a spectrum
def remove_sinusoid(sx,sy,typ):

	#Plot the original data
	figure()
	plot(sx,sy,'g',label='Original spectrum')
	
	# Choosing only the fitted region
	fitreg=x_mask(sx,xlolim,xuplim,siglolim,siguplim,excl_sig='n')
	sx_fitreg=sx[fitreg]
	sy_fitreg=sy[fitreg]
	
	# Choosing only the x-values between the fitted region and the [lower/upper] signal boundary
	msklo=x_mask(sx,xlolim,siglolim,siglolim,siguplim,excl_sig='n')
	msklo_all=x_mask(sx,min(sx)-spec_res,siglolim,siglolim,siguplim,excl_sig='n')
	
	mskhi=x_mask(sx,siguplim,xuplim,siglolim,siguplim,excl_sig='n')
	mskhi_all=x_mask(sx,siguplim,max(sx)+spec_res,siglolim,siguplim,excl_sig='n')
	
	sxlo=sx[msklo]
	sylo=sy[msklo]
	sxhi=sx[mskhi]
	syhi=sy[mskhi]
	
	
	# If blanking of the signal region is required
	# sx_ma=masked_between(sx,siglolim,siguplim)				
	# sy_ma=ma.MaskedArray(sy,sx_ma.mask)
	# sy_ma2=ma.filled(sy_ma,nan)
	# Note, wasn't fitting if there were 'nan's involved
	
	# The function to fit (sinusoid)
	fitfunc = lambda p, sx: p[0]*cos(2*pi/p[1]*sx+p[2]) + p[3]*sx
	
	# Initial guess for the parameters
	if typ=='f':
		p0 = [0.001,5.0,5.0,0.0]
	if typ=='m':
		p0 = [1E8,5.0,5.0,0.0]
	
	# Fit sinusoid to low and high freq part separately
	plo,fitlo=fit_function(sxlo,sylo,p0,fitfunc)
	print 'lofit:', plo
	phi,fithi=fit_function(sxhi,syhi,p0,fitfunc)
	print 'hifit:',phi
	
	# Define a sinusoid with the above fitted parameters across the whole freq range
	fitlo_all=fitfunc(plo,sx)
	fithi_all=fitfunc(phi,sx)
	plot(sx,fitlo_all,'m',label='Lower fit')
	plot(sx,fithi_all,'r',label='Upper fit')
	
	# Sum the fits within the signal region
	sigreg=x_mask(sx,siglolim,siguplim,siglolim,siguplim,excl_sig='n')
	fit_sigreg=fitlo_all[sigreg]+fithi_all[sigreg]
	plot(sx[sigreg],fit_sigreg,'b',label='Sum')
	
	fit_full=concatenate((fitlo_all[msklo_all],fit_sigreg,fithi_all[mskhi_all]),axis=0)
	
	sy_new=sy-fit_full
	
	plot(sx,fit_full,'k',label='Final fit')
	graph_props(typ)
	
	# Plot the original spectrum and the original minus the final fit
	figure()
	plot(sx,sy,'g',label='Original')
	plot(sx,sy_new,'k',label='Original-fit')
	
	graph_props(typ)
	
	return sy_new

#########################################################################

# Function to fit a gaussian profile to a spectrum	
def fit_gauss(sx,sy,typ):

	#Plot the original data
	figure()
	plot(sx,sy,'k',label='Original spectrum')
	
	# The function to fit (gaussian)
	fitfunc = lambda p, sx: p[0]*exp(-(sx-p[1])**2/(2.0*p[2]**2)) # +p[3]
	
	# Initial guess for the parameters
	if typ=='f':
		p0 = [0.07,1420.4,1.0]
	if typ=='m':
		p0 = [1E9,1420.4,1.0]
	
	# Fit sinusoid
	pout,yfit=fit_function(sx,sy,p0,fitfunc)
	print 'Amp=%.3e' %pout[0],'Mean=%.3f' %pout[1],'Sigma=%.4f' %pout[2]
	
	plot(sx,yfit,'r',label='Fit')
	
	#yout=sy-yfit
	#plot(sx,yout,'g',label='Original - Fit')
	
	graph_props(typ)
	
	return yfit,pout #yout
	

################# END ############################################

# Eg. of how to use this program independently

#spectra=genfromtxt('../hipass/mn_hicat_full_stack_flux.txt')
# ion()
# spectra=genfromtxt('../hipass/mn_hicat_full_stack_mass.txt')
# sx=spectra[:,0]
# sy=spectra[:,1]
# #sy_out=remove_sinusoid(sx,sy,'m')
# #y_gfit,pout=fit_gauss(sx,sy,'m')
# mean_stack,[ga,gm,gs]=fit_gauss(sx,sy,'m')
# for ns in [0.6744897501960817, 0.8416212335729142, 1.0364333894937896, \
# 1.2815515655446004, 1.6448536269514726, 1.9599639845400543, \
# 2.575829303548901]:
# 	f_min=gm-ns*gs
# 	f_max=gm+ns*gs
# 	gaus_mask=x_mask(sx,gm-ns*gs,gm+ns*gs,siglolim,siguplim,excl_sig='n')
# 	# Double check that correct percentage of flux/mass included
# 	perc=sum(mean_stack[gaus_mask])/sum(mean_stack)*100
# 	
# 	gausrng_lo=freq_to_vel(gm-ns*gs,freq0)
# 	gausrng_hi=freq_to_vel(gm+ns*gs,freq0)
# 	dv=gausrng_lo-gausrng_hi
# 	print 'N=%.2f' %ns, 'frac=%.3f' %perc, 'f_min=%.2f' %f_min,'f_max=%.2f' %f_max, 'DV=%.2f' %dv
# 	siglims_plot(f_min,f_max,Colour='y')
# 	siglims_plot(siglolim,siguplim,Colour='b')