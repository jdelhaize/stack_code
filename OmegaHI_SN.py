#!/usr/bin/env python

# This code takes an input spectrum and calculates its S/N based on the integral/peak flux in the signal region and the rms of the rest of the spectrum

from numpy import *
from pylab import *
from utils import siglims_plot
from utils import integrate_spectrum
from OmegaHI_plots import get_plot_props

#from OmegaHI_control import *
#from OmegaHI_vars import * #siglolim,siguplim,xlolim,xuplim,spec_res

def SN_calc(x,y,siglolim,siguplim,xlolim,xuplim,spec_res,type,indiv='n',shuf_dir=None,nosig_msk=None): #note that shuf_dir should be a lambda function (can replace x,y by det_stack if opening array)

# 	#****************************
# 	# Function to integrate a spectrum over a range
# 	def integrate_spect(spect,sigreg):
# 		MHI=ma.sum(spect[sigreg])*spec_res #Check if 0 needed
# 		return MHI
# 	#****************************
	
	#****************************
	# Load the normal stacked spectrum
	#****************************
# 	det=genfromtxt(det_stack)
# 	
# 	x,y=det[:,0],det[:,1]


	scal,ylab = get_plot_props(type)

	figure()
	#plot(x,y*scal)

	sigreg=where(logical_and(x<siguplim,x>siglolim)) # The integration range
	
	if indiv=='n':
		# Find the integrated mass/flux of detection stack
		integ_sig=integrate_spectrum(y,sigreg[0],spec_res)
	
	# Find the peak flux of the stack
	peak_flux=ma.max(y[sigreg[0]],axis=0)
	#print peak_flux
	
	# Load the shuffled stacks
	mhiary=[]
	stdary=[]
	
	if indiv=='n':
		for n in range(1,11):
		
			#print n
			n=str(n)
			shuf=genfromtxt(shuf_dir(n))
			
			# Make sure spectral axis same as for detection stack
			xs,ys=shuf[:,0],shuf[:,1]
			if all(xs==x)=='False':
				print 'WARNING: x-axis mismatch!'
			
			plot(xs,ys*scal)
				
			# Integrated mass/flux in signal region
			mhi_n=integrate_spectrum(ys,sigreg[0],spec_res)  # Notice the zero!
			mhiary+=[mhi_n]
			
			# rms in signal region
			std_n=std(ys[sigreg[0]])
			stdary+=[std_n]
			
		# The average rms of the shufflez stacks
		std_noise=mean(stdary)  # Should this be mean or std???
		#print 'Got here'
		std_noise_err=std(stdary)
			
	elif indiv=='y':
	
		nonsigreg=where(logical_or(x<siglolim,x>siguplim))
		std_noise=ma.std(y[nonsigreg[0]],axis=0) # NOTE: To find a fair noise level, calculate from the full spectrum, not just xlolim to xuplim
		#std_noise=ma.std(y[nosig_msk],axis=0) #ie. only between xlolim and xuplim
		#print std_noise
	
	#****************************
	# Calculate the noise and S/N levels and print
	#****************************
	
	SN_peak=peak_flux/std_noise
	#print SN_peak
	
	if indiv=='n':
	
		# The integrated noise is the stdev of the integrated flux/masses of the shuffled spectra
		integ_noise=std(mhiary)
		SN_integ=integ_sig/integ_noise
		
		print 'Integrated signal =', '%e' %integ_sig
		print 'Integrated noise =', '%e' %integ_noise
		print 'Integrated S/N =', "%.2f" %SN_integ
		
		print 'Peak signal =', '%e' %peak_flux
		print 'Stdv noise =', '%e' %std_noise
		print 'Peak S/N =', "%.2f" %SN_peak

		print
	
	#****************************
	# Make the plot pretty
	#****************************
	plot(x,y*scal,'k',lw=2)
	xlim(xlolim,xuplim)
	xlabel(r'$\rm{Frequency\ (MHz)}$')
	ylabel(ylab)
	
	#grid()
	
	if indiv=='n':
		return integ_sig,integ_noise,SN_integ,peak_flux,std_noise,std_noise_err,SN_peak
	if indiv=='y':
		return peak_flux,std_noise,SN_peak



