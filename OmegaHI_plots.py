#!/usr/bin/env python

from numpy import *
from pylab import *
from utils import *
from OmegaHI_vars import *
import utils

freqlab=texstr('Frequency\ (MHz)')

#************
# Flux stack
#************

def flux_stack(x,mn_flux,md_flux,siglolim,siguplim):
	
	# Normal plots	(Comment this section out if using movie plots section below)
	#plot(x,md_flux*1000,'g',label='Median')	
	plot(x,mn_flux*1000,'k',label='Mean')
	xlabel(freqlab)
	ylabel(r'$\rm{Average\ flux\ density\ (mJy)}$')
	#title('Flux stack')
	#legend()
	#grid()
	if noshift!='y':
		siglims_plot(siglolim,siguplim)

def get_plot_props(type):
	if type=='f':
		scal=1000
		ylab=texstr('Average\ flux\ density\ (mJy)')
	if type=='m':
		scal=1E-9
		ylab=r'$\rm{\left\langle M_{HI} \right\rangle}\ (\times 10^{9}\ h^{-2}\ M_{\odot}\ \rm{MHz}^{-1})$'
		#ylab=r'$\rm{Average\ Mass}\ (\times 10^{9}\ h^{-2}\ M_{\odot}\ \rm{MHz}^{-1})$'
	if type=='mmstar':
		ylab=r'$\left\langle \rm{M_{HI}}/\mathcal{M}\right\rangle\ (\rm{MHz}^{-1})$'
		scal=1.0
	if type=='ml':
		ylab=r'$\rm{\left\langle M_{HI}/L_{r}}\right\rangle\ (M_{\odot}\, L_{\odot}^{-1}\, \rm{MHz}^{-1})$'
		scal=1.0
	return scal,ylab		
		
def plot_stack(x,y,type,siglolim,siguplim):	
	#plot(x,md_flux*1000,'g',label='Median')
	scal,ylab=get_plot_props(type)
	
	plot(x,y*scal,'k',label='Mean',lw=1.5)
	xlabel(freqlab)
	ylabel(ylab)
	if noshift!='y':
		siglims_plot(siglolim,siguplim,Colour='r',lw=1.5)
	if fit=='y':
		xlim(xlolim,xuplim)
	hlines(0,min(xlim()),max(xlim()), colors='k', linestyles='dotted')
	vlines(freq0,min(ylim()),max(ylim()),linestyle='dotted',color='k')
	#title('Flux stack')
	#legend()
	#grid()
	
	
def plot_stack_2xax(x,y,type,siglolim,siguplim):	

	fig = figure()
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twiny()
	
	scal,ylab=get_plot_props(type)
	
	ax1.plot(x,y*scal,'k',label='Mean')
	ax1.set_xlabel(freqlab)
	ax1.set_ylabel(ylab)
	if noshift!='y':
		siglims_plot(siglolim,siguplim,Colour='r')
	if fit=='y':
		ax1.set_xlim(xlolim,xuplim)
	hlines(0,min(xlim()),max(xlim()), colors='k', linestyles='dotted')
	vlines(freq0,min(ylim()),max(ylim()),linestyle='dotted',color='k')
	
	ax2.set_xlim(ax1.get_xlim())
	origx = ax1.get_xticks()
	
	def tick_function(X):
		V = freq_to_vel(X,freq0)
		return [r"$%.1f$" % z for z in V]

	ax2.set_xticklabels(tick_function(origx))

	ax2.set_xlabel(r"$\rm(v\ (km/s))$")
	
	#title('Flux stack')
	#legend()
	#grid()
	
def plot_indiv(x,y,type,siglolim,siguplim):
	
	scal,ylab=get_plot_props(type)
	#plot(x,y[:,0:700]*1000)
	plot(x,y*scal)
	xlabel(freqlab)
	ylabel(ylab)
	siglims_plot(siglolim,siguplim,Colour='k')
	
	
def plot_movie_stack(x,y,type,siglolim,siguplim):	
	scal,ylab=get_plot_props(type)
	
	plot(x,y*scal,'k',label='Mean')
	xlabel(freqlab)
	ylabel(ylab)
	ylim(-4,4)
	if fit=='y':
		xlim(xlolim,xuplim)
	hlines(0,min(xlim()),max(xlim()), colors='k', linestyles='dotted')
	vlines(freq0,min(ylim()),max(ylim()),linestyle='dotted',color='k')
	#title('Flux stack')
	#legend()
	#grid()
	
	
# #***************************
# # If making plots for movie:
# #***************************
	# plot(x,md_flux*1000,'r',label='Median')
	# yrng=6.0
	# ylim(-1*yrng,yrng)
	# xlim(1400,1440)
	# title('N={0}'.format(int(Num)))
	# vlines(freq0,-1*yrng,yrng,linestyle='dashed')
	# figname='../noise_movie/g9_blank/stack_g9_{0}.png'.format(flag)
	# savefig(figname)
	
	
#*****************
# Weight function
#*****************
def plot_Nweight(x,weightfn_flux,ylab):	
	#plot(x,weightfn_mass,'k',label=r'$\rm{Mass}$')
	plot(x,weightfn_flux,'m',label=r'$\rm{Flux}$',lw=1.75)
	xlabel(freqlab)
	#ylabel('Number of galaxies stacked')
	ylabel(ylab)
	#title('Weight function of stacks')
	#legend(loc=0)
	grid()
	
#************
# Mass stack
#************
def mass_stack(x,mean_stack,med_stack,siglolim,siguplim):		
	plot(x,mean_stack,'k',label=r'$\rm{Mean}$')
	plot(x,med_stack,'r',label=r'$\rm{Median}$')
	if noshift!='y':
		siglims_plot(siglolim,siguplim)
	xlabel(freqlab)
	ylabel(r'$\rm{Average\ Mass}\ (M_{\odot})$')
	#title('Mass stacks')
	legend()
	grid()

#***********************
# Plot noise scale factor
#************************
def plot_noise(x,noisefac,ylab,col):	
	plot(x,noisefac,col,lw=1.75)
	grid()
	xlabel(freqlab)
	ylabel(ylab)

def plot_noise2(x,noisefac1,noisefac2,ylab,col1,col2,l1,l2):	
	plot(x,noisefac1,col1,lw=1.75,label=l1)
	plot(x,noisefac2,col2,lw=1.75,label=l2)
	grid()
	xlabel(freqlab)
	ylabel(ylab)
	legend(loc=0)

	
#***********************
# Plot binned MHI
#************************	
def bin_plot(midbin_ary,height,bin_incr,nobin,ylab,xlab):
	# Plot bar chart
	errorbar(midbin_ary,height,xerr=bin_incr/2.,ecolor='k',fmt='o',elinewidth=1.5)	
	
	# Make it pretty and include labels
	xlim(min(midbin_ary)-bin_incr/2.-0.001,max(midbin_ary)+bin_incr/2.+0.001)
	xlabel(xlab)
	ylabel(ylab)
	#title('SGP with blanking')
	grid()	
	
	# Include reference value (eg. no binning)
	hlines(nobin,min(xlim()),max(xlim()),'r',linestyle='dashed')
	
	
	
	
def bin_plot2(meanbin,height,rmsbin,OHI_err,nobin,ylab,xlab):
	# Plot bar chart
	errorbar(meanbin,height,yerr=OHI_err,xerr=rmsbin,ecolor='k',color='k',fmt='o',elinewidth=1.5)	
	
	# Make it pretty and include labels
	#xlim(min(midbin_ary)-bin_incr/2.-0.001,max(midbin_ary)+bin_incr/2.+0.001)
	xlabel(xlab)
	ylabel(ylab)
	#grid()	
	
	# Include reference value (eg. no binning)
	#hlines(nobin,min(xlim()),max(xlim()),'r',linestyle='dashed')
	
def stack_overlay(binnum,x,y,rnd_ary,plotpar):

	scal,ylab=get_plot_props(plotpar)
	
	def stack_var(ary,colour,lb,scal):
		dat=genfromtxt(ary)
		x=dat[:,0]
		y=dat[:,1]*scal
		plot(x,y,colour,label=lb,lw=1.5)
	
	plot(x,y*scal,'k',label=r'$\rm{Original}$',lw=1.5)#r'$3\sigma$')	
	stack_var(rnd_ary,'g',r'$\rm{Control}$',scal)#r'$1\sigma$',scal)
	
	#grid()
	xlabel(freqlab)
	ylabel(ylab)
	#legend(loc=0)
	xlim(xlolim,xuplim)
	
	siglims_plot(siglolim_orig,siguplim_orig,Colour='r',Line='--',lw=1.5)
	hlines(0,min(xlim()),max(xlim()), colors='k', linestyles='dotted')
	vlines(freq0,min(ylim()),max(ylim()),linestyle='dotted',color='k')




#***********************
# Luminosity and z error plots
#************************

def plot_z_v_zerr(z, zerr):
	#figure()
	plot(z,zerr*100/z,'ko')
	xlabel(r'$z$')
	ylabel(r'$\%\rm{\ error\ on\ z}$')

def plot_L_v_Ler(L_err, Lum):
	#figure()
	plot(Lum,L_err*100/Lum,'ko')
	xlabel(r'$\rm{Luminosity}$')
	ylabel(r'$\%\rm{\ error\ on\ L}$')
	xscale('log')


def Ler_hist(L_err, Lum):
	#figure()
	hist(L_err*100/Lum,bins=20)
	xlabel(r'$\%\rm{\ error\ on\ L}$')
	ylabel(r'$\rm{Count}$')
	
def plot_z_v_Ler(L_err, Lum, z):
	#figure()
	plot(z,L_err*100/Lum,'ko')
	xlabel(r'$z$')
	ylabel(r'$\%\rm{\ error\ on\ L}$')


#***********************
# Ldist comparison plots
#************************
# Compare the corrected and uncorrected luminosity distributions. Moved here from OHI_newlum.py	
def lum_compare(newL,lum_raw):	
	figure()
	print 'min(newL), max(newL) = ', ma.min(newL,axis=0), ma.max(newL,axis=0)
	#lmsk = where(log10(lum_raw)>4)[0] # To temporarily mask bad L to see the distbn of just the good values
	#lum_raw, newL = lum_raw[lmsk], newL[lmsk]
	hist(log10(newL),bins=20,label=utils.texstr('Corrected'),histtype='step',lw=3,ec='r')
	hist(log10(lum_raw),bins=20,label=utils.texstr('Original'),histtype='step',lw=3,ec='k')
	legend(loc=0,frameon=False)
	xlabel(r'$\log(L/(L_{\odot}\,h^{-2}))$')
	ylabel(r'$\rm{Count}$')
	
def sm_compare(new_sm,sm_raw):	
	figure()
	print 'min(new_sm), max(new_sm) = ', ma.min(new_sm,axis=0), ma.max(new_sm,axis=0)
	print 'shape(new_sm), shape(sm_raw)', shape(new_sm), shape(sm_raw)
	#lmsk = where(log10(sm_raw)>4)[0] # To temporarily mask bad L to see the distbn of just the good values
	#sm_raw, new_sm = sm_raw[lmsk], new_sm[lmsk]
	hist(log10(new_sm),bins=20,label=utils.texstr('Corrected'),histtype='step',lw=3,ec='r')
	hist(log10(sm_raw),bins=20,label=utils.texstr('Original'),histtype='step',lw=3,ec='k')
	legend(loc=0,frameon=False)
	xlabel(r'$\log(M^*/(M_{\odot}\,h^{-2}))$')
	ylabel(r'$\rm{Count}$')
		
# Plot original and adjusted luminosity error histogram
def lum_err_compare(newL_err,L_err):
	figure()
	hist(log10(newL_err),bins=20,label=utils.texstr('Corrected'),histtype='step',lw=3,ec='r')
	hist(log10(L_err),bins=20,label=utils.texstr('Original'),histtype='step',lw=3,ec='k')
	legend(loc=0)
	xlabel(r'$log(\Delta L/L_{\odot})$')
	ylabel(r'$\rm{Count}$')
	
#***********************
# Completeness scale factor plots
#************************	
def boostf_plot(boost_fact, SMboost_fact):
	#figure()
	#hist(SMboost_fact,bins=30,ec='r',lw=3,histtype='step',range=(min(boost_fact),max(boost_fact)),label=r'$\mathcal{M}$')
	hist(boost_fact,bins=20,ec='k',lw=3,histtype='step',label=r'$L_r$')
	hist(SMboost_fact,bins=20,ec='r',lw=3,histtype='step',label=r'$\mathcal{M}$')
	legend(loc=0,frameon=False)
	xlabel(texstr('Completeness\ scaling\ factor'))
	ylabel(texstr('Count'))

#***********************
# Plot source positions with z on colour axis
#************************	
def plot_position(ra,dec,z):	
	fig=figure()
	ax1=subplot(111)
	sc = scatter(ra,dec,c=z,s=15,hold=True,marker='o',edgecolor='None')#,vmin=0.2,vmax=1.2)
	xlabel(r'$\rm{Right\ Ascension}$')
	ylabel(r'$\rm{Declination}$')
	cax=colorbar(sc)
	cax.ax.set_ylabel(r'$\rm{Redshift}$')
	