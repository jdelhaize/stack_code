#!/usr/bin/env python

# Need to make this multidimensional for processing gama

from numpy import *
from pylab import *
from utils import *
from constants import *
from OmegaHI_params import *
from OmegaHI_vars import *
if OHIvDV=='y':
	from OmegaHI_control import *

#print siglolim,siguplim

def spectral_confusion( nlw, nv, vref, dv=600 ):

#Uncomment this if OHIvDL=='y'
	if OHIvDV=='y':
		siguplim,siglolim = choose_siglims(dv_to_df(dv))
	else:
		print 'Nothing'

	#print 'vref, nv, nlw =',vref, nv, nlw

	DLO = ( freq_to_vel(freq0,freq0) - freq_to_vel(siguplim,freq0) ) *1
	DHI = ( freq_to_vel(siglolim,freq0) - freq_to_vel(freq0,freq0) ) *1
	
	# Separation between the velocity of the central galaxy in question (vref) and the lower/upper edge of the potentially-confused galaxy profile
	sep1 = (nv - nlw/2) - vref
	sep2 = vref - (nv + nlw/2)
	
	
	# Find the fraction of the galaxy profile that sits within the integration bounds (assumes a boxcar bj lum profile)
	if nv > vref and any(sep1<DHI):
		lum_fact = ( (vref+DHI) - (nv - nlw/2) )/nlw
	elif nv < vref and any(sep2<DLO):
		lum_fact = ( (nv + nlw/2) - (vref-DLO) )/nlw
	else:
		lum_fact = zeros(shape(nlw))
		
	if survey=='2df':
		if lum_fact > 1:
			lum_fact = 1
	if survey=='g':
		lum_fact[ where(lum_fact > 1) ] = 1
		lum_fact[ where(lum_fact < 0) ] = 0 # Because some bands satisfy sep1<DHI or sep2<DLO criteria, and others don't. Ones that don't are negative.

	return lum_fact


# Plot an individual spectrum, eg. to show an example of confusion
def plot_indiv_spec(id):

	ind=where(ID==id)#300593.0)
	
	figure()
	plot(x,y[:,ind][:,0,0]*1000,'k')
	
	ylim(-200,300)
	xlim(1410,1430)
	
	ymin=min(ylim())
	ymax=max(ylim())
	vlines(siglolim,ymin,ymax,color='r',linestyle='--',lw=1.2)
	vlines(siguplim,ymin,ymax,color='r',linestyle='--',lw=1.2)
	
	hlines(0,min(xlim()),max(xlim()), colors='k', linestyles='dotted')
	vlines(freq0,min(ylim()),max(ylim()),linestyle='dotted',color='k')
	
	ylabel(r'$\rm{Flux\ density\ (mJy)}$')
	xlabel(r'$\rm{Frequency\ (MHz)}$')
	
	
# Plot the redshift distribution of the input sample and the full 2dfgrs Nz (scaled to the same field size)
def nz_n_full2df( z_in ):

	figure()

	#sgp = genfromtxt('../catalogues/p669+2df_zcat+i.txt',usecols = (24))
	sgp = z_in
	full = genfromtxt('../catalogues/best.observations_pksrange.txt',usecols = (24))
	
	rng=[min(sgp), max(sgp)]
	
	outs = hist(sgp, histtype='step', bins=20, range=rng, ec='k', lw=2)
	[Ns,zs] = outs[0], outs[1]
	
	
	outf = hist(full, histtype='step', bins=20, range=rng, ec='g', lw=2)
	[Nf,zf] = outf[0], outf[1]
	
	wdth = zf[1] - zf[0]
	scale = 2000/42.0
	binl=list(zf[0:-1])
	
	figure()
	#bar(zf[0:-1], Nf/scale, width = wdth, ec='g', fc='None')
	#bar(zf[0:-1], Ns, width = wdth, ec='k', fc='None')
	
	x=[]
	for i in range(len(binl)):
		x+=[binl[i]]*(Nf[i]/scale)
	x=array(x)
	
	hist(x, histtype='step', bins=20, range=rng, ec='g', lw=2)
	hist(sgp, histtype='step', bins=20, range=rng, ec='k', lw=2)
	
	xlabel(texstr('Redshift'))
	ylabel(texstr('Count'))
	






		
	