#!/usr/bin/env python

from numpy import *
from pylab import *
import pickle
from utils import *
from OmegaHI_vars import *
from OmegaHI_control import *
#from OmegaHI_parent import cat_num
from OmegaHI_newlum import *
import OmegaHI_TF as ohitf
import OmegaHI_plots as ohiplots
import utils_gama

import random

# cat_num_tst=pickle.load(open('cat_num.txt','r'))

# *******************************
# Function to define the columns
# *******************************

if survey=='2df' and hicatq=='n':
	def def_col_2df(sources):
		#global ID,z,zQ,appmag,ra,dec
		ID=sources[:,0]
		z=sources[:,24]
		zQ=sources[:,26]
		appmag=sources[:,17]
		ra=sources[:,37]
		dec=sources[:,38]
		colour=sources[:,21]-sources[:,22]
		semimajor= sources[:,40] #sources[:,37]+0.1*sources[:,37] #sources[:,40]
		semiminor= sources[:,41] #sources[:,37] #sources[:,41]
		sb_bj=sources[:,21]
		sb_r=sources[:,22]
		return ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r
		
if survey=='2df' and hicatq=='y':
	def def_col_2df(sources):
		#global ID,z,zQ,appmag,ra,dec
		ID=sources[:,0]#[:,-1]
		z=sources[:,36]
		zQ=sources[:,27]*0+5 # Dummy!!!
		appmag=sources[:,28] # Dummy!!!
		ra=sources[:,34]
		dec=sources[:,35]
		colour=sources[:,28]+0.1-sources[:,28] #Dummy!! =0.1
		semimajor=sources[:,8]+30 #Dummy!!
		semiminor=sources[:,8] #Dummy!!
		sb_bj=sources[:,28] #Dummy!!
		sb_r=sources[:,28]+0.1 #Dummy!!
		return ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r


# def def_col_gama(sources):
# 	#global ID,z,zQ,absmag_ugriz,ra,dec,absmagerr
# 	ID=sources[:,0]
# 	z=sources[:,12]
# 	zQ=sources[:,13]
# 	absmag_ugriz=sources[:,6]
# 	ra=sources[:,2]
# 	dec=sources[:,3]
# 	absmagerr_ugriz=sources[:,6] # Columns needed for error propagation
# 	#absmagerr_ugriz=APerr_ugriz-A_ugriz-k_ugriz-DM_cp #Don't need to adjust an error by a constant. For these other columns see OmegaHI_bigcat3.py
# 	semimajor=sources[:,0]+100+30 #Dummy!!
# 	semiminor=sources[:,0]+100 #Dummy!!
# 	colour=sources[:,6]-sources[:,7]
# 	return ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor

def def_col_gama(sources):
	#global ID,z,zQ,absmag_ugriz,ra,dec,absmagerr
	ID=sources[:,0]
	if zversion=='31':
		z=sources[:,34] #34 (autoz) or 38 (tc16)
		zQ=sources[:,13]
	if zversion=='16':
		z=sources[:,38] #34 (autoz) or 38 (tc16) #99 for photoz
		zQ=sources[:,39]
	ra=sources[:,2]
	dec=sources[:,3]

	colour=sources[:,6]-sources[:,7] # Dummy
	
# 	if genshz=='y': # Note: These are all temporary dummy values
# 		semimajor=sources[:,2]*0+20 #From ApMatchedCatv02
# 		semiminor=sources[:,2]*0+10
# 		absmag_u,absmag_g,absmag_r,absmag_i,absmag_z=sources[:,5],sources[:,6],sources[:,7],sources[:,8],sources[:,9]
# 		absmag_ugriz=array([absmag_u,absmag_g,absmag_r,absmag_i,absmag_z]).T
# 		absmagerr_u,absmagerr_g,absmagerr_r,absmagerr_i,absmagerr_z=sources[:,5],sources[:,6],sources[:,7],sources[:,8],sources[:,9]
# 		absmagerr_ugriz=array([absmagerr_u,absmagerr_g,absmagerr_r,absmagerr_i,absmagerr_z]).T	
# 	
# 	else:
	semimajor=sources[:,96]*0.339 #From ApMatchedCatv02
	semiminor=sources[:,97]*0.339 
	fluxscale = sources[:,42]
	fluxscale = restrict_fluxscale(fluxscale)
	#fluxscale = ma.masked_greater(fluxscale,2)
	fluxscale.shape=(len(fluxscale),1)
	
	absmag_u,absmag_g,absmag_r,absmag_i,absmag_z=sources[:,72],sources[:,77],sources[:,82],sources[:,87],sources[:,92]
	absmag_ugriz=array([absmag_u,absmag_g,absmag_r,absmag_i,absmag_z]).T
	absmag_ugriz = gama_absmag_fluxscale(absmag_ugriz, fluxscale)

	absmagerr_u,absmagerr_g,absmagerr_r,absmagerr_i,absmagerr_z=sources[:,73],sources[:,78],sources[:,83],sources[:,88],sources[:,93] # Columns needed for error propagation
	absmagerr_ugriz=array([absmagerr_u,absmagerr_g,absmagerr_r,absmagerr_i,absmagerr_z]).T
	
	return ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor
	
# NOTE: Use this only for Laura's stuff!
# def def_col_2df(sources):
# 	#global ID,z,zQ,appmag,ra,dec
# 	ID=sources[:,0]
# 	z=sources[:,5]
# 	zQ=sources[:,6]
# 	appmag=sources[:,6]
# 	ra=sources[:,2]
# 	dec=sources[:,4]
# 	return ID,z,zQ,appmag,ra,dec


# *************************************
# Function for uniform source selection
# *************************************
# def source_select(msk_ary,sources,Lum,DL,y):
# 	#global sources,Lum,DL,y
# 	sources=sources[msk_ary]
# 	Lum=Lum[msk_ary]
# 	DL=DL[msk_ary]
# 	y=y[:,msk_ary]
# 	return sources,Lum,DL,y

# def source_select(msk_ary, pdict, y):
# 
# 	pdict2={}
# 	
# 	for i in pdict:
# 	
# 		print shape(pdict[i])
# 		i2=pdict[i][msk_ary]
# 		print shape(i2)
# 		
# 		pdict2[i]=i2
# 		print i, shape(pdict2[i])
# 
# 		
# 	y=y[:,msk_ary]
# 	
# 	return pdict2

def source_select(msk_ary, pdict, y):
	
	pdict2={}
	
	for i in pdict:
	
		i2=pdict[i][msk_ary]
		pdict2[i]=i2

		
	y=y[:,msk_ary]
	
	return pdict2,y

if survey=='2df':
	params_lst = [ 'sources', 'Lum', 'DL', 'L_err', 'z_orig', 'zQ_orig', 'semimajor', 'semiminor', 'absmag_uncorr', 'Babsmag','lum_raw' ]
else:
	params_lst = [ 'sources', 'Lum', 'DL', 'L_err', 'z_orig', 'zQ_orig', 'semimajor', 'semiminor', 'absmag_uncorr','lum_raw', 'stellar_mass', 'stellar_mass_raw', 'beamsum' ]

# Could maybe just concatenate these parameters into one big array? Then it would be easy to deselect the relevant rows all at the same time. Could redefine them at the very end to get them out of the array and in to individual variables.

# ******************************************************************************
# ******************************************************************************
# Main function 


# *************************************
# Function for all source selection processes
# *************************************
def sourceselect(sources,spectra,cat_num,dv):
	
	
	if nocontin=='y':
		contin=genfromtxt(contincat)#,delimiter=cdelim)
	
	
	# ****************************
	# Define freq (x) and flux (y)
	# ****************************
	x=spectra[1::,0] # Frequency (the first column, except the first value - 9999)
	x.shape=(len(x),1)
	ids=spectra[0,1::]
	y=spectra[1::,1::] # The array of intensities for each spectrum
	
	if opcat_subset=='y': # I suspect this has never worked!
		id_subset=genfromtxt(subsetcat)
		indexop=array(id_subset,dtype='int64')
		print shape(indexop),indexop[0]
		sources=sources[indexop-1,:]
		print shape(sources)
	
	# If using a subset of an original bigcat
	if bigcat_subset=='y':
		#index=array(sources[:,0],dtype='int64')
		#print index[0], index[-1], shape(index) # for hicat nonconf+nonex
		if opcat_subset=='y':
			index=indexop
			y=y[:,index-1] 
		else: # This is where the index column is defined!!!!!!!!!
			index=array(sources[:,indexcol],dtype='int64') # for normal hipass
			y=y[:,index-1]  
	
	# Remove 'spectra' from memory to free up some space
	del spectra
	
	# *******************
	# Define the columns
	# *******************
	if survey=='2df':
		ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
	if survey=='g':
		ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)
		#def_col_gama(sources)

	# *******************
	# Define the original z and zQ (important for shufflez files)
	# *******************
	
	if genshz=='y': #and binning=='y' and binby='L':
		z_orig=sources[:,z_orig_col]
		zQ_orig=sources[:,q_orig_col]
	else:
		z_orig=z
		zQ_orig=zQ

	# ****************************
	# Calculate the raw luminosity
	# ****************************
	if survey=='2df':
		DL=z_to_DL(z_orig)
		absmag=app_to_absmag(DL,appmag)
		absmag_uncorr=absmag
		
		lum_b=absmag_to_lum(absmag,SolarMag_b)
		Lum=lum_b
		lum_raw=Lum
		L_err=lum_error(z, zerr, lum_raw, merr)
		
		# Find the absolute magnitude in B-band (for use in TF relation below)
		Bappmag = bj_to_B(sb_bj, sb_r)
		Babsmag=app_to_absmag(DL,Bappmag)
		#figure()
		#hist(Babsmag,bins=20)
		#xlabel(r'$M_{B}$')
		
	if survey=='g':
		DL=z_to_DL(z_orig)
		#absmag_ugriz=ma.masked_less(absmag_ugriz,-98.)
		absmag_ugriz = ma.masked_invalid(absmag_ugriz)
		#absmag_ugriz = ma.masked_greater(absmag_ugriz,-15)
		#absmag_ugriz = ma.masked_greater(absmag_ugriz,-6)
		absmag_uncorr=absmag_ugriz
		
		lum_ugriz=absmag_to_lum(absmag_ugriz,SolarMag_ugriz)
		Lum=lum_ugriz
		lum_raw = Lum
		
		#lum_raw=absmag_ugriz #temporary only!
		#L_err=lum_error(z, zerr, lum_raw, merr) # needs changing!
		L_err = Merr_to_Lerr(absmagerr_ugriz, lum_raw) # needs changing!
		
		# Will need additional DL corrections here (eg. Virgo etc.)
		
	# **************************************************************
	# Calculate estimated line width using the Tully-Fisher relation
	# **************************************************************
		
	print 'Prior to selection: lum mean, median = ', mean(lum_raw), median(lum_raw)

	if confuse=='y':

		print 'Prior to source selection:'
		if survey=='2df':
			linewidth = ohitf.tf_W(semimajor, semiminor, tf_offset, tf_slope, Babsmag) #absmag_uncorr is the bj band magnitude 	
		if survey=='g':
			linewidth = ohitf.tf_W(semimajor, semiminor, tf_offset, tf_slope, absmag_uncorr)
		print 'mean(linewidth) = ', mean(linewidth)
		print 'shape(linewidth) = ',shape(linewidth)
	
	
	# ********************************************		
	# Define the stellar mass
	# ********************************************
	if survey=='g':
		fluxscale = sources[:,42] # Note: This is also defined in def_col_gama
		fluxscale = restrict_fluxscale(fluxscale)	
		figure()
		hist(fluxscale,bins=30,histtype='step',lw=3,ec='c')
		xlabel(texstr('Flux scale'))
		
		stellar_mass = sources[:,47]
		stellar_mass_error = sources[:,48]
		figure()
		hist(stellar_mass,bins=30,histtype='step',lw=3,ec='g',label='Original')
		
		stellar_mass = gama_stellarmass_fluxscale(stellar_mass, fluxscale)
		hist(stellar_mass,bins=30,histtype='step',lw=3,ec='y',label='Flux scaled')
		legend(loc=0)
		xlabel(texstr('log(Stellar\ Mass)'))
		
		stellar_mass_raw = stellar_mass.copy()

		
	
		tmp = ones((len(stellar_mass),5))
		stellar_mass_raw_2d = stellar_mass_raw
		stellar_mass_error_2d = stellar_mass_error
		stellar_mass_raw_2d.shape = (len(stellar_mass),1)
		stellar_mass_error_2d.shape = (len(stellar_mass_error),1)
		
		stellar_mass_raw_2d = stellar_mass_raw_2d * tmp
		stellar_mass_error_2d = stellar_mass_error * tmp
		
	# ********************************************
	# Adjust luminosity (and stellar mass) for confusion if required
	# ********************************************

	# Calculate the completeness boost factors

	if do_boostf=='y':
		if survey=='2df':
			boost_fact = boostf[surveynm] # Defined in ohi_vars.py
		if survey=='g':
			boost_fact = utils_gama.completenessF_gama_r(z, appmag_lim, -24, SolarMag_ugriz[mbd[magband]])
			print 'mean(boost_fact) = ', mean(boost_fact)
			print 'Warning: Completeness factor being calculated for r-band ONLY!!!'
			#boostfact.shape=(len(boost_fact),1)
			#%figure()
			#%hist(boost_fact,bins=20,fc='y')
			#%xlabel(texstr('L_r\ completeness\ scaling'))
			
			SMboost_fact, SMmin = utils_gama.SMcompletenessF_gama(z, 10**14)
			sm_msk = where(10**stellar_mass>SMmin)
			print 'Len(SM_msk) = ', shape(sm_msk[0])
			
			print 'mean(SMboost_fact) = ', mean(SMboost_fact)
			
			#%figure()
			#%hist(SMboost_fact,bins=20,fc='m')
			#%xlabel(texstr('\mathcal{M}\ completeness\ scaling'))
			
	else:		
		boost_fact=1.0
		SMboost_fact = 1.0


	# Do the confusion correction
	
	if confuse=='y':
				
		if survey=='2df':
			Lum, L_err, beamsum = account_for_confusion(cat_num, ID, z_orig, ra, dec, lum_raw, L_err, dv, linewidth, boost_fact, doSM='n')
		
		if survey=='g':
			Lum, L_err, beamsum, stellar_mass, dummy_smerr = account_for_confusion(cat_num, ID, z_orig, ra, dec, lum_raw, L_err, dv, linewidth, boost_fact, doSM='y', stellar_mass = 10**stellar_mass_raw_2d, sm_err = 10**stellar_mass_error_2d, SMboostf=SMboost_fact)
			#Lum=ma.masked_less(Lum,10**5)
			
			print 'shape(beamsum)', shape(beamsum)
			
			stellar_mass, dummy_smerr = log10(stellar_mass[:,0]), log10(dummy_smerr[:,0])
			#stellar_mass=ma.masked_less(Lum[:,2],10**5)
	else:
		beamsum=stellar_mass*0.0+1.0
	
####################################################	
	
	# Put any source selection conditions below here
	
	
	if binning=='y' and binby=='z':
		z_orig=z 
		zQ_orig=zQ
		# ie. don't want to use the original redshifts for source selection from this point down, want to use the shuffled z.
	
	# Define the properties to replace in source selection
	
	pdict = {  }
	for i in params_lst:
		pdict[i]=eval(i)	

	
	# # *******************************************
	# # Exclude bad gama abs_mag sources
	# # *******************************************
		
	if absmag_msk=='y':		
		
		#Mmsk=where(logical_and(sources[:,92]<-15,sources[:,92]>-99))
		Mmsk=where(logical_and(sources[:,92]<-8,sources[:,92]>-99))
		
		print "Bad absmag:",len(Mmsk[0]), "out of", len(sources), "spectra passed"
		
		pdict,y = source_select(Mmsk[0],pdict,y)
		for key in pdict:
			exec ("%s = pdict[key]" % key)			
		#sources,Lum,DL,y=source_select(colourmsk[0],sources,Lum,DL,y)
		
		# Re-define the columns
		if survey=='2df':
			ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
		if survey=='g':
			ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)
		#print shape(z)
		print 'final', shape(z)
		
	
	#************************************************************
	# Exclude sources with negative and/or bad quality redshifts
	#************************************************************
	
	#*******************
	# 1) Negative redshifts and stars in hipass sample
	#*******************
	if survey=='2df' and pksdat=='h' and hicatq!='y':
	
		zgood = where( z_orig>0.0025 ) #749km/s (was 0.001 = 300km/s)
		
		print "Stars:",len(zgood[0]), "out of", len(sources), "spectra passed"
		
		pdict,y = source_select(zgood[0],pdict,y)
		for key in pdict:
			exec ("%s = pdict[key]" % key)
		
# 		sources,Lum,DL,y = source_select(zgood[0],sources,Lum,DL,y)
# 		z_orig, zQ_orig = z_orig[zgood[0]], zQ_orig[zgood[0]]
		
		# *******************
		# Re-Define the columns
		# *******************
		if survey=='2df':
			ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
		if survey=='g':
			ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)

	#*******************
	# Deselect edge of GAMA field (to minimise confusion correction edge effects)
	#*******************
	
	if edgecut=='y':
		if survey=='g' and field=='combo' or field=='hipass':
		
# 		g9chop = where( logical_and( logical_and( dec<=(dmax9-edgew), dec>=(dmin9+edgew) ),\
# logical_and( ra<=(ramax9-edgew), ra>=(ramin9+edgew) )))
# 		g15chop = where( logical_and( logical_and( dec<=(dmax15-edgew), dec>=(dmin15+edgew) ),\
# logical_and( ra<=(ramax15-edgew), ra>=(ramin15+edgew) )))

			edgemsk = where( logical_or(\
logical_and( logical_and( dec<=(dmax9-edgew), dec>=(dmin9+edgew) ),\
logical_and( ra<=(ramax9-edgew), ra>=(ramin9+edgew) )),\
logical_and( logical_and( dec<=(dmax15-edgew), dec>=(dmin15+edgew) ),\
logical_and( ra<=(ramax15-edgew), ra>=(ramin15+edgew) ))))
			
	
		# Print number of sources excluded
		print "Edge limit:",len(edgemsk[0]), "out of", len(dec), "spectra passed"
	
		pdict,y=source_select(edgemsk[0],pdict,y)	
		for key in pdict:
			exec ("%s = pdict[key]" % key)
	
		
		# Re-Define the columns
		if survey=='2df':
			ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
		if survey=='g':
			ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)

	#*******************
	# Select good gama region (not extra GAMA II dec strip)
	#*******************
	
	if deccut=='y':
		if survey=='g' and field=='9':
			gooddec=where(dec>-1.)
		if survey=='g' and field=='15':
			gooddec=where(dec<2.)
		if survey=='g' and field=='combo':
			gooddec=where(logical_and(dec<2.,dec>-1))
			print "Are you sure you want to use these dec ranges??"
	
		# Print number of sources excluded
		print "Bad dec:",len(gooddec[0]), "out of", len(dec), "spectra passed"
	
		pdict,y=source_select(gooddec[0],pdict,y)	
		for key in pdict:
			exec ("%s = pdict[key]" % key)
	
		
		# Re-Define the columns
		if survey=='2df':
			ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
		if survey=='g':
			ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)

	#*********************************
	# Use correct GAMA DMU
	#*********************************
			
	if survey=='g':
		survey_class=sources[:,24]
		if dmu=='I':
			class_cut = where(survey_class>=6)
		if dmu=='II':
			class_cut = where(survey_class>=3)
		print "Survey class cut: ",len(class_cut[0]), "out of", len(survey_class), "spectra passed"
		
		pdict,y=source_select(class_cut[0],pdict,y)	
		for key in pdict:
			exec ("%s = pdict[key]" % key)
		
		# Re-Define the columns
		if survey=='2df':
			ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
		if survey=='g':
			ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)


	#*******************
	# Select only GAMA sources below the r<19.8 completeness limit
	#*******************
			
	if completeness_cut=='y':
		rpetro=sources[:,5]
		if dmu=='I':
			goodr = where(rpetro<=19.4)#19.4)
		if dmu=='II':
			goodr = where(rpetro<=19.8)
		print "Completeness cut: ",len(goodr[0]), "out of", len(rpetro), "spectra passed"
		
		pdict,y=source_select(goodr[0],pdict,y)	
		for key in pdict:
			exec ("%s = pdict[key]" % key)
		
		# Re-Define the columns
		if survey=='2df':
			ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
		if survey=='g':
			ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)

		
	#*******************
	# 2) Low quality redshifts
	#*******************
	if survey=='g':
		goodQ=where(zQ_orig>0.)
	if survey=='2df':
		goodQ=where(zQ_orig>2.)
	
	# Print number of sources excluded
	print "Bad Q(z):",len(goodQ[0]), "out of", len(zQ_orig), "spectra passed"
	
	pdict,y=source_select(goodQ[0],pdict,y)	
	for key in pdict:
		exec ("%s = pdict[key]" % key)
		#vars()[key]=pdict[key]

# 	sources,Lum,DL,y=source_select(goodQ[0],sources,Lum,DL,y)
# 	z_orig, zQ_orig = z_orig[goodQ[0]], zQ_orig[goodQ[0]]
	
	# *******************
	# Re-Define the columns
	# *******************
	if survey=='2df':
		ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
	if survey=='g':
		ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)
		
		
	# **************************
	# 3) Exclude continuum sources
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
		
		angsep_msk=ma.masked_less(angsep,beamwidth/sepfac) #Note the /2 (sepfac)!!!
		msk=ma.any(angsep_msk.mask,axis=1)
		non_cont=where(msk==False) 
	
		# Print number of sources excluded
		print "Contin:",len(non_cont[0]), "out of", len(sources), "spectra passed"
	
		[pdict,y] = source_select(non_cont[0],pdict,y)
		for key in pdict:
			exec ("%s = pdict[key]" % key)
		
# 		sources,Lum,DL,y=source_select(non_cont[0],sources,Lum,DL,y)
# 		z_orig, zQ_orig = z_orig[non_cont[0]], zQ_orig[non_cont[0]]
	
	
		# *******************
		# Re-Define the columns
		# *******************
		if survey=='2df':
			ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
		if survey=='g':
			ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)

	# # *******************************************
	# # 4) Exclude full spectrum if source in RFI zone
	# # *******************************************
		
	if rfi_mask=='y':
		for i in bad_i:
			#print shape(z_orig), '=shape z'
			badf=eval('bad'+str(i))			
			badz=freq_to_z(array(badf),freq0)		
			rfimsk=where(logical_or(z_orig>badz[0],z_orig<badz[1]))
			
			#print shape(badz), 'passed badz'
			
			pdict,y = source_select(rfimsk[0],pdict,y)
			for key in pdict:
				exec ("%s = pdict[key]" % key)
# 			sources,Lum,DL,y=source_select(rfimsk[0],sources,Lum,DL,y)
# 			z_orig, zQ_orig = z_orig[rfimsk[0]], zQ_orig[rfimsk[0]]
			
			# Re-define the columns
			if survey=='2df':
				ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
			if survey=='g':
				ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)
			
		print "RFI zone:",len(sources), "spectra remain"
		
	# # **************************
	# # Select one single spectrum to play with
	# # **************************
	
	if singspec=='y':
		ind=where(ID==76496.0)
		
		pdict,y = source_select(ind[0],pdict,y)
		for key in pdict:
			exec ("%s = pdict[key]" % key)
		#sources,Lum,DL,y=source_select(ind[0],sources,Lum,DL,y)
	
	# *******************
	# Re-Define the columns
	# *******************
	if survey=='2df':
		ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
	if survey=='g':
		ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)
		
	# # *******************************************
	# # Exclude sources with bad colours
	# # *******************************************
		
	if colour_msk=='y':		
		
		colourmsk=where(logical_and(colour>-1,colour<2))
		
		print "Bad colours:",len(colourmsk[0]), "out of", len(sources), "spectra passed"
		
		pdict,y = source_select(colourmsk[0],pdict,y)
		for key in pdict:
			exec ("%s = pdict[key]" % key)			
		#sources,Lum,DL,y=source_select(colourmsk[0],sources,Lum,DL,y)
		
		# Re-define the columns
		if survey=='2df':
			ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
		if survey=='g':
			ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)
		#print shape(z)
		print 'final', shape(z)
		


# Choose only emission line sources
# 		
# 	abemma=sources[:,27]
# 	emi_msk=where(abemma==1)[0]
# 	sources,Lum,DL,y=source_select(emi_msk,sources,Lum,DL,y)
# 	
# 	# *******************
# 	# Re-Define the columns
# 	# *******************
# 	if survey=='2df':
# 		ID,z,zQ,appmag,ra,dec=def_col_2df(sources)
# 	if survey=='g':
# 		ID,z,zQ,absmag_ugriz,ra,dec,absmagerr_ugriz=def_col_gama(sources)
# 	


	# ************************************************************************************
	# Plot the comparison of the raw and adjusted luminosity distributions if required
	# ************************************************************************************
	
	print 'mean Lum = ',mean(Lum)
	if confuse=='y':
		if survey=='2df':
			ohiplots.lum_compare(Lum,lum_raw)
		if survey=='g':
			print 'shape(Lum),shape(lum_raw) = ',shape(Lum[:,2]),shape(lum_raw[:,2])
			#%ohiplots.lum_compare(Lum[:,2],lum_raw[:,2])
			#%ohiplots.sm_compare(10**stellar_mass, 10**stellar_mass_raw)

		


	# ***********************************************
	# Return all parameters needed in parent code
	# ***********************************************
	
	if survey=='2df':
		return sources,ID,z,zQ,appmag,ra,dec,x,y,ids,DL,Lum,goodQ,colour,L_err,semimajor,semiminor,absmag_uncorr,Babsmag, lum_raw #absmag,lum_b,lum_raw,contin
	elif survey=='g':
		return sources,ID,z,zQ,absmag_ugriz,ra,dec,absmagerr_ugriz,x,y,ids,DL,Lum,goodQ,colour,L_err,semimajor,semiminor,absmag_uncorr,stellar_mass, lum_raw, stellar_mass_raw[:,0], beamsum #absmag_ugriz,lum_ugriz,lum_raw,,contin
	
