#!/usr/bin/env python

from numpy import *
from pylab import *
from constants import *
from OmegaHI_vars import *
from utils import *

# NOTE: This is only configured for f_or_m=='f'

#Use this if individual bounds have been defined with TF, otherwise use SN_calc in OmegaHI_SN.py
def calc_indivsn(x_ma, y_ma, indiv_siglolim, indiv_siguplim):

	# Mask out the signal region
	x_mask_sig = masked_between(x_ma,indiv_siglolim,indiv_siguplim)
	y_mask_sig = ma.MaskedArray(y_ma,x_mask_sig.mask)

	# Find the noise level (stdv)
	# NOTE: To find a fair noise level, calculate from the full spectrum, not just xlolim to xuplim (but excluding signal region)
	#std_noise=ma.std(y_mask_sig, axis=0) # Only excludes 'linewidth' region
	
	nonsigreg=where(logical_or(x_ma[:,0]<siglolim,x_ma[:,0]>siguplim))
	std_noise=ma.std(y_ma[nonsigreg[0]],axis=0) # Excludes full 'stacked profile width' region 
	
	# Mask everything except the signal region
	x_mask_notsig = masked_between(x_ma, min(x_ma[:,0]), indiv_siglolim)
	x_mask_notsig = masked_between(x_mask_notsig, indiv_siguplim, max(x_ma[:,0]))
	y_mask_notsig = ma.MaskedArray(y_ma,x_mask_notsig.mask)
	
	# Find the peak flux of the stack
	peak_flux=ma.max(y_mask_notsig,axis=0)
	
	
	# Calculate the noise and S/N levels	
	SN_peak=peak_flux/std_noise

	return peak_flux,std_noise,SN_peak



def indiv_SN_investig(x,yorig,y_fit,ids,indiv_S,indiv_N,indiv_SN,z, indiv_siglolim, indiv_siguplim, cut=3.7):

	notmaskd=where(indiv_S.mask==False)
	print 'No. unmasked S/N values = ', shape(notmaskd), 'out of',len(indiv_S)
	figure()
	hist(indiv_SN[notmaskd[0]],bins=20)
	xlabel('S/N')
	ylabel('# spectra')
	
	dubious=where(indiv_SN[notmaskd[0]]>=cut)
	print 'Num S/N>=',cut,' = ',shape(dubious[0])
	dub_ids=ids[notmaskd][dubious]
	print dub_ids
	
	f=z_to_freq(z,freq0)
	print shape(ids),shape(z)
	
	for idx in dub_ids:
		figure()
		ind=where(ids==idx)
		
		plot(x,yorig[:,ind[0]]*1000,'b')
		plot(x,y_fit[:,ind[0]]*1000,'g')
		
		ylim(-30,30)
		xlim(xlolim,xuplim)
		nm="%.1f" %(indiv_SN.data[ind][0])
		title('ID='+str(idx)+'  S/N='+nm+'  freq='+str(f[ind][0]))
		vlines(freq0,min(ylim()),max(ylim()),linestyles='dashed')
		vlines([siglolim,siguplim],min(ylim()),max(ylim()),linestyles='dashed',colors='y')
		vlines(indiv_siglolim[ind[0]],min(ylim()),max(ylim()),linestyles='dashed',colors='m')
		vlines(indiv_siguplim[ind[0]],min(ylim()),max(ylim()),linestyles='dashed',colors='m')
		hlines(0,min(xlim()),max(xlim()),linestyles='dotted',colors='k')
		xlabel('Frequency (MHz)')
		ylabel('Average flux density (mJy)')