#!/usr/bin/env python

# Calculates the average confused MHI and L(Bj) to get OmegaHI for the SGP and GAMA data
# Note: requires 'jacinta' user on Endymion to be mounted (not if bigcat code used first)

#*******************************************************************************
# 									SET UP									   *
#*******************************************************************************

#****************************
# Load packages and functions
#****************************
from numpy import *
from scipy import *
from scipy import integrate
from pylab import *
from utils import *
from pyfn_find_confused import *


#*****************
# Ask what to run
#*****************

survey=raw_input('Use gama or sgp? (g or s): ')
nocontin=raw_input('Exclude spectra near continuum sources? (y or n): ')
fit=raw_input('Fit baseline before stacking? (y or n): ')
blank=raw_input('Blank bad freq channels? (y or n): ')
weighting=raw_input('Use weighting? (y or n): ')

confuse=raw_input('Account for confusion? (y or n): ')
if confuse=='y':
	vel_confused=float(raw_input('--> Use which velocity separation?: '))
	nonconf='n'
else:
	nonconf=raw_input('--> Use only unconfused sources? (y or n):')


#*****************
# Define constants
#*****************
h=1.0
H0=100.*h
G=4.3*10**-9
c=299792.458
freq0=1420.405751786 # Rest freq MHz
DH=(9.26*10.**25.)/h #In metres
Omega_M=0.3
Omega_lambda=0.7
Omega_k=1.-Omega_M-Omega_lambda
MpctoMet=3.0857*10**22
SolarMag_b=5.27
SolarMag_ugriz=array([6.39,5.07,4.62,4.52,4.48])

#**************************
# Define luminosity density
#**************************
if survey=='s':
	rhoL=1.82*10**8*h #Bj luminosity density
if survey=='g':
	rhoL=array([0.8846,1.0068,1.3325,1.6130,1.9083])*10**8 #ugriz. r-band from most recent Driver paper.  [0.8846,1.0068,2.37,1.6130,1.9083]
	
#******************************
# Define radio survey parameters
#******************************
beamwidth=(15./60.) #in degrees
spat_confused=beamwidth
#vel_confused=250.0 #in km/s Changed to raw input
spec_res=0.0625 #Spectral channel resolution in MHz

# *************************
# Define integration ranges
# *************************

xlolim=1425  #1400. # For fitting
xuplim=1440. # For fitting

if survey=='s':
	siglolim=1418.77 # For integrating
	siguplim=1421.45 # For integrating

if survey=='g':
	siglolim=1416.4 # For integrating
	siguplim=1421.7 # For integrating
	
	#Temp only!!!
	#siglolim=1418.77 # For integrating
	#siguplim=1421.45 # For integrating

# ***********************************
# Non-confused input list and spectra
# ***********************************

if nonconf=='y':
    
	if survey=='s':
		unconfused=genfromtxt('../Omega_HI/sgp_250_unconfused.txt',dtype=int32)
		spectra=genfromtxt('../Omega_HI/bigcat_sgp.txt',usecols=spect_cols)
	if survey=='g':
		unconfused=genfromtxt('../Omega_HI/g9_250_unconfused.txt',dtype=int32)
		spectra=genfromtxt('../Omega_HI/bigcat_g9.txt',usecols=spect_cols)
	spect_cols=tuple(insert(unconfused+1,0,0))
	
# **************************************
# Input source list and spectra (all)
# **************************************
if survey=='s':
	sources=genfromtxt('/Users/jacintadelhaize/workdir/stacking/catalogues/pyselection2.txt.combined') 
	#sources=genfromtxt('/home/jdelhaize/workdir/DATA/pyselection2.txt.combined')
	#sources=genfromtxt('/Users/jacintadelhaize/workdir/stacking/Evenf_trial/pyselection2_evenf3.txt')
	# Epeius

	spectra=genfromtxt('../Omega_HI/bigcat_sgp.txt')
	#spectra=genfromtxt('/home/jdelhaize/workdir/DATA/bigcat_sgp.txt')
	#spectra=genfromtxt('/Users/jacintadelhaize/workdir/stacking/Evenf_trial/bigcat_sgp_evenf.txt')
	
	#contin=genfromtxt('../Omega_HI/sgp_continuum_pos.csv',delimiter=',')
	#contin=genfromtxt('/home/jdelhaize/workdir/DATA/sgp_continuum_pos.csv',delimiter=',')

	
	ra=ra_hextodeg(sources[:,10],sources[:,11],sources[:,12])
	dec=sources[:,13]-sources[:,14]/60.-sources[:,15]/3600.
	
	
if survey=='g':

	#sources=genfromtxt('/home/jdelhaize/workdir/DATA/gama9v6.txt.combined')
	sources=genfromtxt('/Users/jacintadelhaize/workdir/stacking/catalogues/gama9v6.txt.combined')
	#sources=genfromtxt('../binning/g9_vhihalomass_match.txt') 
	
	#spectra=genfromtxt('../binning/g9_vhihalomass_bigcat.txt')
	#spectra=genfromtxt('/home/jdelhaize/workdir/DATA/bigcat_g9.txt')
	spectra=genfromtxt('../Omega_HI/bigcat_g9.txt')
	
	#contin=genfromtxt('../Omega_HI/g9_continuum_pos.csv',delimiter=',')
	#contin=genfromtxt('../Omega_HI/g9_NVSS_contin_200.txt')
	#contin=genfromtxt('/home/jdelhaize/workdir/DATA/g9_NVSS_contin_200.txt')
	
	ra=sources[:,3]
	dec=sources[:,4]

# NOTE: Put any source selection conditions here (and in spectra section).

print 'Sources and spectra loaded.'
# **************************
# Select N sources to stack
# **************************
sources_all=sources

flag=0
#flag=-1
#for logN in range(100,110,5):
#for logN in range(100,365,5):
for logN in range(100,380,5): #gama

	flag+=1
	
	logNflt=logN/100.
	Num=round(10**logNflt)

	from random import sample
	Nselection=sample(range(0,len(sources_all[:,0])),int(Num))
	sources=sources_all[Nselection]
	
	# **************************
	# Exclude continuum sources
	# **************************
	if nocontin=='y':
		ra_r=ra*pi/180.0
		dec_r=dec*pi/180.0
		
		#Cra=ra_hextodeg(contin[:,0],contin[:,1],contin[:,2])
		#Cdec=contin[:,3]+contin[:,4]/60.+contin[:,5]/3600.
		Cra=contin[:,0]
		Cdec=contin[:,1]
		
		Cra_r=Cra*pi/180.0
		Cdec_r=Cdec*pi/180.0
		
		ra_big=(ra_r*ones((len(Cra_r),len(ra_r)))).T
		dec_big=(dec_r*ones((len(Cdec_r),len(dec_r)))).T
		angsep=angular_sep(Cra_r,ra_big,Cdec_r,dec_big)
		del ra_big,dec_big
		
		angsep_msk=ma.masked_less(angsep,beamwidth)
		msk=ma.any(angsep_msk.mask,axis=1)
		non_cont=where(msk==False) 
		sources=sources[non_cont[0]]
	
	if nonconf=='y':
		sources=sources[unconfused]
	
	# NOTE: Comment out spectra loaded in here if using 'unconfused'
	
	# ****************************
	# Remove low quality redshifts
	# ****************************
	
	# if survey=='g':
	# 	sources=sources[where(sources[:,7]>2.)]
	
	# NOTE: Make sure source selection columns ABOVE this point.
	
	# ******************
	# Define the columns
	# ******************
	
	if survey=='s':
		ID=sources[:,0]
		z=sources[:,24]
		appmag=sources[:,17]
		ra=ra_hextodeg(sources[:,10],sources[:,11],sources[:,12])
		#dec=dec_hextodeg(sources[:,13],sources[:,14],sources[:,15])
		dec=sources[:,13]-sources[:,14]/60.-sources[:,15]/3600.
	
	if survey=='g':
		ID=sources[:,2]
		z=sources[:,6]
		absmag_ugriz=sources[:,83:88]
		ra=sources[:,3]
		dec=sources[:,4]
	
	
	
	#*******************************************************************************
	# 						CALCULATIONS ON CATALOGUE COLUMNS					   *
	#*******************************************************************************
	
	
	# ***************
	# Define velocity	
	# ***************
	vel=z*c
	
	# ******************************
	# Convert ra and dec to radians
	# ******************************
	ra_r=ra*pi/180.0
	dec_r=dec*pi/180.0
	
	
	# *********************************
	# Calculate the luminosity distance
	# *********************************
	
	E=lambda z: sqrt(Omega_M*(1.+z)**3.+Omega_k*(1.+z)**2.+Omega_lambda)
	integral=[]
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
	
	# Will need additional DL corrections here (eg. Virgo etc.)
	
	# ************************
	# Calculate the luminosity
	# ************************
	
	if survey=='s':
		absmag=appmag-5.*log10(DL*10.**5.)
		lum_b=10.**((SolarMag_b-absmag)*0.4)
		
	if survey=='g':
		#Blank sources that don't have a magnitude entry
		absmag_ugriz=ma.masked_greater(absmag_ugriz,9998.)
		lum_ugriz=10.**((SolarMag_ugriz-absmag_ugriz)*0.4)
		
	#	lum=concatenate([lum_u,lum_g,lum_r,lum_i,lum_z],axis=1)
	
	#**Prepare empty arrays**
	lums_conf=[]
	
	
	
	# ******************************************************************************
	# 							WORK ON SPECTRA									   *
	# ******************************************************************************
	
	
	# ****************************
	# Define freq (x) and flux (y)
	# ****************************
	x=spectra[1::,0] # Frequency (the first column, except the first value - 9999)
	x.shape=(len(x),1)
	y=spectra[1::,1::]
	
	# Put source selection conditions here too (and in column section above).
	y=y[:,Nselection]
	if nocontin=='y':
		y=y[:,non_cont[0]]
	
	# ****************************
	# mask default flux values
	# ****************************
	y_ma=ma.masked_equal(y,-999.)
	
	# Remove 'spectra' from memory to free up some space
	#del spectra
	
	# *************************************************
	# Fit and subtract polynomial to each spectra
	# *************************************************
	if fit=='y':
		y_0=ma.filled(y_ma,fill_value=0.0) # Set blank values to 0. Need to do this?
		p=ma.polyfit(x[:,0],y_0,4)
		newy=y-polyval(p,x)
		#newy=cubicfit(x[:,0],y_0)
		y_ma=ma.MaskedArray(newy,y_ma.mask)
	
	# ****************************
	# Blank bad frequency channels
	# ****************************
	
	if blank=='y':
	
		# Define bad channels
		bad1_lo=1262.000*(1+z)		# Broad-band RFI
		bad1_hi=1275.250*(1+z)
		bad2_lo=1279.875*(1+z)
		bad2_hi=1280.250*(1+z)
		bad3_lo=1312.375*(1+z)
		bad3_hi=1312.813*(1+z)
		bad4_lo=1315.813*(1+z)
		bad4_hi=1316.313*(1+z)
		bad5=1253.05*(1+z)			# low freq cube edge
		#bad5=1300.00*(1+z)
		bad6=1366.90*(1+z)    		# high freq cube edge
		
		# Expand x array to same shape as y (to use as blanking template)
		x_ma=ones(shape(y))*x
		
		# Masking broad-band RFI, satellite RFI and edges of freq range
		x_ma=masked_between(x_ma,bad1_lo,bad1_hi)
		x_ma=masked_between(x_ma,bad4_lo,bad4_hi)
		x_ma=ma.masked_less(x_ma,bad5)
		x_ma=ma.masked_greater(x_ma,bad6)
		
		# Additional satellite RFI regions in GAMA data
		if survey=='g':
			x_ma=masked_between(x_ma,bad2_lo,bad2_hi)
			x_ma=masked_between(x_ma,bad3_lo,bad3_hi)
		
		# Applying these masked regions to y
		y_ma=ma.MaskedArray(y_ma,x_ma.mask)
	
	
	# *****************************************
	# Plot some individual shifted flux spectra
	# *****************************************
	# figure()
	# plot(x,y_ma[:,0:10])
	
	
	# ***********************
	# Weighting and scaling
	# ***********************
# 	
# 	# Finding the noise level of each spectra
# 	sigma=ma.std(y_ma,axis=0)
# 	num_stack=ma.count(y_ma,axis=1)
# 	num_stack.shape=(len(num_stack),1)
# 	# y_ma=y_ma/sigma         # Optional direct scaling by sigma
# 	
# 	# Choose weight factor!!!!
# 	#weight=DL**(-2)
# 	#weight=sigma
# 	weight=(num_stack)**0.5
# 	
# 	
	
	# ***********************
	# Convert to HI mass axis
	# ***********************
	
	# massy=y_ma*(2.356*10**5)*DL**2 		#This is for y in Jy*km/s
# 	
# 	if weighting=='y':
# 		y_ma=y_ma*weight
# 	
# 	mass_array=y_ma*(5.0*10**7)*DL**2
# 	
# 	# if weighting=='y':
# 	# 	mass_array=mass_array*weight
# 	
	
	# ******************************
	# Calculate MHI for each spectra
	# ******************************
# 	
# 	sigreg=where(logical_and(x<siguplim,x>siglolim))
# 	MHI=ma.sum(mass_array[sigreg[0],:],axis=0)*spec_res # Array of indiv HI masses.
# 	
# 	# NOTE: the DL could be multiplied here instead!!! Check whether [0] needs to be here!!! 
# 	
	
	# *********************
	# Find confused sources
	# *********************
	if confuse=='y':
	
		# Call find_confused from pyfn_find_confused.py
		angsep,count_confused = find_confused(vel,ra_r,dec_r,spat_confused,vel_confused)
		
		# Call confused_luminosity from pyfn_find_confused.py to calculate confused luminosities. 
		 # NOTE: Variables lum_b and lum_ugriz are overwritten here
		 
		if survey=='s':
			lum_b=confused_luminosity(lum_b,angsep,beamwidth)
			
		if survey=='g':
			bands={0:'u',1:'g',2:'r',3:'i',4:'z'}
			
			for iter in range(0,5):
				name='lum_{0}'.format(bands[iter])
				vars()[name]=confused_luminosity(lum_ugriz[:,iter],angsep,beamwidth)
				# Outputs lum_u, lum_g etc.
		
	# 		lum_u=confused_luminosity(lum_ugriz[:,0],angsep,beamwidth)
	# 		lum_g=confused_luminosity(lum_ugriz[:,1],angsep,beamwidth)
	# 		lum_r=confused_luminosity(lum_ugriz[:,2],angsep,beamwidth)
	# 		lum_i=confused_luminosity(lum_ugriz[:,3],angsep,beamwidth)
	# 		lum_z=confused_luminosity(lum_ugriz[:,4],angsep,beamwidth)
	
	
			lum_ugriz=vstack((lum_u,lum_g,lum_r,lum_i,lum_z)).T
			#del lum_u,lum_g,lum_r,lum_i,lum_z
			#lum_ugriz=lum_r
			
		print 'Average number of confused sources:',ma.mean(count_confused)
		
		# Plot histogram of confusion
		figure()
		hist(count_confused,bins=30) #accounting for double-up
		xlabel('Number of confused sources')
		ylabel('Count')
	
	#*******************************************************************************
	# 							CALCULATING FINAL VALUES						   *
	#*******************************************************************************
	
	
	#************************
	# Find average luminosity
	#************************
	# 
	# if survey=='s':
	# 	Lav=mean(lum_b)
	# if survey=='g':
	# 	Lav=ma.mean(lum_ugriz,axis=0)
	# 
	# print 'mean Luminosity =',print_sci(Lav,survey) #'%e' %Lav
	# 
	# #******************
	# # Find average mass
	# #******************
	# 
	# if weighting=='y':
	# 	MHIav=ma.sum(MHI)/sum(weight)
	# else:
	# 	MHIav=ma.mean(MHI)
	# print 'mean MHI =','%e' %MHIav
	# print 'median MHI =','%e' %ma.median(MHI)
	# 
	# #**************
	# # Find OmegaHI
	# #**************
	# 
	# OmegaHI=OmegaHI_calc(MHIav/Lav,rhoL)
	# print "OmegaHI =",print_sci(OmegaHI,survey)
	
	
	#**************
	# <M/L> Method
	#**************
	
	# print "MASS TO LIGHT METHOD:"
	# # Find mass-to-light ratio
	# if survey=='s':
	# 	MassToLight=MHI/lum_b
	# 	Av_MassToLight=ma.mean(MassToLight)
	# 	
	# if survey=='g':
	# 	def ML(lum):
	# 		MasstoLight=MHI/lum
	# 		Av_MasstoLight=ma.mean(MassToLight)
	# 		
	# 	Av_MLu=ML(lum_ugriz[:,0])
	# 	Av_MLg=ML(lum_ugriz[:,1])
	# 	Av_MLr=ML(lum_ugriz[:,2])
	# 	Av_MLi=ML(lum_ugriz[:,3])
	# 	Av_MLz=ML(lum_ugriz[:,4])
	# 
	# 
	# OmegaHI2=OmegaHI_calc(Av_MassToLight,rhoL)
	# print "OmegaHI =",print_sci(OmegaHI2,survey)
	
	
	#*******************************************************************************
	# 									PLOTTING								   *
	#*******************************************************************************
	
	
	#************
	# Flux stack
	#************
	yrng=6.0
	mn_flux=ma.mean(y_ma,axis=1)
	md_flux=ma.median(y_ma,axis=1)
	figure()
	#plot(x,mn_flux*1000,'g',label='Mean')
	plot(x,md_flux*1000,'r',label='Median')
	ylim(-1*yrng,yrng)
	xlim(1400,1440)
	#siglims_plot(siglolim,siguplim)
	xlabel('Frequency (MHz)')
	ylabel('Average Flux (mJy)')
	title('N={0}'.format(int(Num)))
	vlines(freq0,-1*yrng,yrng,linestyle='dashed')
	#legend(loc=0)
	#grid()
	
	figname='../noise_movie/g9_blank/stack_g9_{0}.png'.format(flag)
	savefig(figname)
	
	#*****************
	# Weight function
	#*****************
	
# 	weightfn_mass=ma.count(mass_array,axis=1)
# 	weightfn_flux=ma.count(y_ma,axis=1)
# 	
# 	figure() # Plot in new figure window
# 	plot(x,weightfn_mass,'k',label='Mass')
# 	plot(x,weightfn_mass,'m',label='Flux')
# 	xlabel('Frequency (MHz)')
# 	ylabel('Number of galaxies stacked')
# 	title('Weight function of stacks')
# 	legend(loc=0)
# 	grid()
	
	#************
	# Mass stack
	#************
	
	# if weighting=='y':
	# 	mean_stack=ma.sum(mass_array,axis=1)/sum(weight) # #ie. no division by n
	# 	med_stack=ma.median(mass_array,axis=1)/sum(weight)
	# else:
	# 	mean_stack=ma.mean(mass_array,axis=1)
	# 	med_stack=ma.median(mass_array,axis=1)
	# 
	# figure()
	# plot(x,mean_stack,'r',label='Mean')
	# plot(x,med_stack,'k',label='Median')
	# siglims_plot(siglolim,siguplim)
	# xlabel('Frequency (MHz)')
	# ylabel(r'Average Mass (M$_\odot$)')
	# #title('Mass stacks')
	# legend(loc=0)
	# grid()
	
	#*********************
	#How to save an array
	#*********************
	
	# a=md_flux
	# a.shape=(len(a),1)
	# b=concatenate((x,a),axis=1)
	# savetxt('../Evenf_trial/overlays/sgp_evenf_md_nob.txt',b)
	
	#*******************************************************************************
	# 							SUM THEN AVERAGE METHOD							   *
	#*******************************************************************************
	
	# MHI_stack_av=ma.sum(mean_stack[sigreg[0]])*spec_res #Check if 0 needed
	# #if weighting=='y':
	# 	#MHI_stack_av=MHI_stack_av/sum(weight)
	# print "Stacked mean mass =",'%e' %MHI_stack_av
	# 
	# OmegaHI_stack=OmegaHI_calc(MHI_stack_av/Lav,rhoL)
	# print "Stacked OmegaHI =",print_sci(OmegaHI_stack,survey)
	# 
	# print 'Max no. stacked =', max(weightfn_mass)
	
	
	#*********************
	# Clip before stacking
	#*********************
	
	# stdvs=ma.std(mass_array,axis=0)
	# STDV=mean(stdvs)
	# mns=ma.mean(mass_array,axis=0)
	# MN=mean(mns)
	# clip=ma.masked_outside(mass_array,MN-2.5*STDV,MN+2.5*STDV)
	# clip_mnstk=ma.mean(clip,axis=1)
	# clip_mdstk=ma.median(clip,axis=1)
	# plot(x,mean_stack,'b')
	# plot(x,med_stack,'g')
	

#*******************************************************************************
# 								END PROGRAM!								   *
#*******************************************************************************


#savetxt('./output_arrays/mass_array_blanked.txt',mass_array,delimiter=',')

#Watch out for units! DH is in meters!!!!
#What's the z-flow thing?







#Fitting polynomial, subtracting and getting OmegaHI estimate

# xlolim=1400. # For fitting
# xuplim=1440. # For fitting
# siglolim2=1418
# siguplim2=1422.5
# y=md_flux
# mska=where(logical_and(x>xlolim,x<siglolim2))
# mskb=where(logical_and(x>siguplim2,x<xuplim))
# msk=concatenate([mska[0],mskb[0]])
# #plot(x[msk],y[msk])
# a,b,c,d=polyfit(x[msk][:,0],y[msk],3)
# newy=y-(a*x[:,0]**3+b*x[:,0]**2+c[:,0]+d)
# plot(x[:,0],newy)
# # vlines(siglolim3,4,-1,color='r',linestyle='--')
# 
# #siglolim3=1418.77
# 
# plot(x[:,0],newy*1000,'k')
# xlim(xlolim,xuplim)
# ylim(-2,2)
# grid()
# xlabel('Frequency (MHz)')
# ylabel('Flux (mJy)')









