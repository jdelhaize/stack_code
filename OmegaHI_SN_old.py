#!/usr/bin/env python

# This code takes an input spectrum and calculates its S/N based on a series of random stacks

from numpy import *
from pylab import *
from utils import siglims_plot

#from OmegaHI_control import *
#from OmegaHI_vars import * #siglolim,siguplim,xlolim,xuplim,spec_res

def SN_calc(x,y,shuf_dir,siglolim,siguplim,xlolim,xuplim,spec_res): #note that shuf_dir should be a lambda function (can replace x,y by det_stack if opening array)

	#****************************
	# Function to integrate a spectrum over a range
	def integrate_spect(spect,sigreg):
		MHI=ma.sum(spect[sigreg])*spec_res #Check if 0 needed
		return MHI
	#****************************
	
	#****************************
	# Load the normal stacked spectrum
	#****************************
# 	det=genfromtxt(det_stack)
# 	
# 	x,y=det[:,0],det[:,1]
	figure()
	plot(x,y*10**-9)
	sigreg=where(logical_and(x<siguplim,x>siglolim)) # The integration range
	
	# Find the integrated mass of detection stack
	integ_sig=integrate_spect(y,sigreg[0])
	
	# Find the peak flux of the stack
	peak_flux=ma.max(y)
	
	# Load the shuffled stacks
	mhiary=[]
	stdary=[]
	for n in range(1,11):
	
		#print n
		n=str(n)
		shuf=genfromtxt(shuf_dir(n))
		
		# Make sure spectral axis same as for detection stack
		xs,ys=shuf[:,0],shuf[:,1]
		if all(xs==x)=='False':
			print 'WARNING: x-axis mismatch!'
		
		plot(xs,ys*10**-9)
			
		mhi_n=integrate_spect(ys,sigreg[0])  # Notice the zero!
		mhiary+=[mhi_n]
		
		std_n=std(ys[sigreg[0]])
		stdary+=[std_n]
		
	#****************************
	# Calculate the noise and S/N levels
	#****************************
	print mhiary
	integ_noise=std(mhiary)
	std_noise=mean(stdary) # Should this be mean or std???
	std_noise_err=std(stdary)
	
	SN_integ=integ_sig/integ_noise
	SN_peak=peak_flux/std_noise
	
	
	#****************************
	#Print the output values
	#****************************
	print 'Integrated signal =', '%e' %integ_sig
	print 'Integrated noise =', '%e' %integ_noise
	print 'Integrated S/N =', SN_integ
	
	print 'Peak signal =', '%e' %peak_flux
	print 'Stdv noise =', '%e' %std_noise
	print 'Peak S/N =', SN_peak
	
	#****************************
	# Make the plot pretty
	#****************************
	xlim(xlolim,xuplim)
	xlabel('Frequency (MHz)')
	ylabel(r'Average mass (M$_\odot$)$\times 10^{9}$')
	grid()
	
	
	return integ_sig,integ_noise,SN_integ,peak_flux,std_noise,std_noise_err,SN_peak



