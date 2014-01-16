#!/usr/bin/env python

# Calculates the average confused MHI and L(Bj) to get OmegaHI for the SGP and GAMA data. Also plots the stacked spectrum, the weight function and the confusion histogram

# Differs from OmegaHI_bigcat4 by splitting bits off into functions in different modules

# Plots that are usually output but which I have temporarily commented out to make the program run faster can be found by searching for '#%'

#****************************
# Load packages and modules
#****************************
from numpy import *
from scipy import *
from scipy import integrate
from scipy.interpolate import interp1d
from pylab import *
ion() # turns interactive mode on
import pickle
import time
import os, sys

rcParams['text.usetex']='True'
rcParams['font.sans-serif']= "Helvetica"
rcParams['axes.labelsize']=17
rcParams['font.size']=14
rcParams['font.weight']=600
	
#from constants import *
from OmegaHI_vars import *
from OmegaHI_control import *
from OmegaHI_sourceselect import *
from OmegaHI_plots import *
from utils import *
from OmegaHI_SN import *
from OmegaHI_indivSN import *
from OmegaHI_nvn import *
from OmegaHI_fitting import *
from OmegaHI_binplot import *
#from OmegaHI_stackoverlay import *
from OmegaHI_TF import *
import OmegaHI_Ffactor as ohi_f
import OmegaHI_other as ohi_other

startt=time.time()	


if genshz=='y':
	cat_num_loop=range(1,11)
	#cat_num_loop=range(8,11)
else:
	cat_num_loop=[""]
for cat_num in cat_num_loop: #use this one if no loop needed
	
	cat_num=str(cat_num)
	print 'Catalogue number:',cat_num
		
	OHI_ary_conf=[]
	OHI_err_ary_conf=[]
	dv_ary=[]
	
	#for dv in range(50,750,50):
	if confuse=='y':
		#conf_loop=[122, 152, 187, 231, 297, 353, 464]
		#conf_loop=range(50,750,50)
		conf_loop=[600]#[2000]#[750] # For 2dFGRS paper decided to settle on dv=600
		#conf_loop=range(90,990,50)
		#conf_loop=range(890,990,50)
	else:
		#conf_loop=[""] # use this if you want to use the siglims defined in OHI_vars
		conf_loop=[600]
	
	for dv in conf_loop:	
	
		print 'Delta V = ',dv
		dv_ary+=[dv]
		
		###
		# Comment out if want to use original siglims for HI integration
		if OHIvDV=='y':
			siguplim,siglolim = choose_siglims(dv_to_df(dv))
			#siglolim=1419.7
		
		print 'Loading sources= ',choose_opcat(cat_num)
		sources=genfromtxt(choose_opcat(cat_num))
		# CAUTION: THIS IS TEMPORARY ONLY! I MISSED A SOURCE IN THE EXTRACTION!
		#sources=sources[1::,:]
		print 'shape sources', shape(sources) 
		
		print 'Loading spectra= ',choose_bigcat(cat_num)
		spectra=genfromtxt(choose_bigcat(cat_num))
		print 'shape spectra', shape(spectra)
			
		# Execute the sourceselect code
		# Note that I have taken logN out of the sourceselect input parameters
		if survey=='2df':
			sources,ID,z,zQ,appmag,ra,dec,x,y,ids,DL,Lum,goodQ,colour,L_err,semimajor,semiminor,absmag_uncorr,Babsmag, lum_raw=sourceselect(sources,spectra,cat_num,dv)
		elif survey=='g':
			sources,ID,z,zQ,absmag_ugriz,ra,dec,absmagerr_ugriz,x,y,ids,DL,Lum,goodQ,colour,L_err,semimajor,semiminor,absmag_uncorr,stellar_mass, lum_raw, stellar_mass_raw, beamsum=sourceselect(sources,spectra,cat_num,dv)
		
		# From OmegaHI_sourceselect.py
 		if survey=='2df':
			params_lst = [ 'sources', 'Lum', 'DL', 'L_err', 'z', 'zQ' , 'semimajor', 'semiminor', 'absmag_uncorr','Babsmag', 'lum_raw']#note:different to as defined in OHI_sourselect as is missing lum_raw. does this matter?
		else:
			params_lst = [ 'sources', 'Lum', 'DL', 'L_err', 'z', 'zQ' , 'semimajor', 'semiminor', 'absmag_uncorr', 'stellar_mass', 'lum_raw','stellar_mass_raw', 'beamsum']
		pdict = {  }
		for i in params_lst:
			pdict[i]=eval(i)
		pdict_orig=pdict.copy()
		print "shape Orig = ", shape(pdict_orig['sources'])
# 		rnd=[]
# 		for rni in range(len(ID)):
# 			rnd+=random.gauss(0,0.03)
# 		rnd=array(rnd)
# 		rnd_big=ones(shape(y))*rnd
# 		y=y+rnd_big
		
		
# 		if survey=='2df' and pksdat=='p' and fluxcalib=='y': # This has moved to below
# 			zcut=freq_to_z(1310,freq0)
# 			msk1285=where(z>zcut)[0]
# 			msk1335=where(z<=zcut)[0]
# 			y[:,msk1285]=y[:,msk1285]*1.15
# 			y[:,msk1335]=y[:,msk1335]*1.02
			
		
		
		# Save the full selection for the logN loop to sample from
		sourcesA,LumA,DLA,yA,pdictA=sources,Lum,DL,y,pdict
		
		NvN_mean=[]
		NvN_std=[]
		
		if runnvn=='y':
			#logNloop=range(200,210,iterat)
			logNloop=range(100,maxN,iterat)
		else:
			logNloop=[""]
			
		for logN in logNloop:
		
			print 'logN: ',logN, 'shzn: ',cat_num
			
			# **************************************************
			# Select N sources to stack (movie or NvN analysis)
			# **************************************************
			# Doing this here so that it won't re-load bigcat every time
			
			if runnvn=='y':
				
				logNflt=logN/100.
				Num=round(10**logNflt)
			
				from random import sample
				Nselection=array(sample(range(0,len(sourcesA[:,0])),int(Num)))	
					
				#pdict,y = source_select(Nselection,pdictA,yA) #or Nselection[0]??
# 				for key in pdict:
# 					exec ("%s = pdict[key]" % key)

# 				sources,Lum,DL,y=source_select(Nselection,sourcesA,LumA,DLA,yA)

# NOTE: This is a massive hack. For some reason the source_select line above isn't working. If you add anything extra into pdict, include it here too.

				# Resetting to full arrays (pre-Nselection)
				if survey=='2df':
					(sources, Lum, DL, L_err, z, zQ , semimajor, semiminor, absmag_uncorr,Babsmag, lum_raw) = (pdictA['sources'], pdictA['Lum'], pdictA['DL'], pdictA['L_err'], pdictA['z'], pdictA['zQ'] , pdictA['semimajor'], pdictA['semiminor'], pdictA['absmag_uncorr'],pdictA['Babsmag'],pdictA['lum_raw'])
				if survey=='g':
					(sources, Lum, DL, L_err, z, zQ , semimajor, semiminor, absmag_uncorr,stellar_mass, lum_raw, stellar_mass_raw, beamsum) = (pdictA['sources'], pdictA['Lum'], pdictA['DL'], pdictA['L_err'], pdictA['z'], pdictA['zQ'] , pdictA['semimajor'], pdictA['semiminor'], pdictA['absmag_uncorr'], pdictA['stellar_mass'],pdictA['lum_raw'],pdictA['stellar_mass_raw'], pdictA['beamsum'])
				
				# Doing this Nselection
				y=yA[:,Nselection]
				
				Lum=Lum[Nselection]
				DL=DL[Nselection]
				semimajor=semimajor[Nselection]
				absmag_uncorr=absmag_uncorr[Nselection]
				zQ=zQ[Nselection]
				sources=sources[Nselection]
				z=z[Nselection]
				semiminor=semiminor[Nselection]
				L_err=L_err[Nselection]		
				lum_raw = lum_raw[Nselection]
				beamsum = beamsum[Nselection]		
				if survey=='2df':
					Babsmag=Babsmag[Nselection]
				if survey=='g':
					stellar_mass = stellar_mass[Nselection]
					stellar_mass_raw = stellar_mass_raw[Nselection]
					


				# Re-Define the columns
				if survey=='2df':
					ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
				if survey=='g':
					ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)
			
			#***********************************************************************
			# 								Set up binning						   *
			#																	   *
			#***********************************************************************
			
			sources_all=sources
			y_all=y
			z_all=z
			Lum_all=Lum
			DL_all=DL
			lum_raw_all = lum_raw
			beamsum_all = beamsum
			if survey=='g':
				stellar_mass_all = stellar_mass
				stellar_mass_raw_all = stellar_mass_raw
			
			flag=0
			ary_OHI=[]
			ary_MHI=[]
			ary_ML=[]
			midbin_ary=[]
			meanbin=[]
			rmsbin=[]
			OHI_err_ary=[]
			ML_err_ary=[]
			colour_ary=[]
			integ_noise_ary=[]
			bin_len_ary=[]
			peaksn_ary=[]
			gastostar_ary=[]
			gastostar_err_ary=[]
			lum_raw_ary = []
			Lum_ary = []
			
			# Define new parameters that don't need to be used previously (not used for source selection)
			if survey=='g':
				
# 				if genshz=='y':
# 					print 'WARNING: Using temporary columns!!!'
# # 					stellar_mass = sources_all[:,20]
# # 					stellar_mass_all = stellar_mass				
# 					uminusr = sources_all[:,6]
# 					uminusr_all = uminusr			
# 					gminusi = sources_all[:,7]
# 					gminusi_all = gminusi
# 					
# 				else:
# 					fluxscale = sources_all[:,42]
# 					#fluxscale = ma.masked_greater(fluxscale,2)
# 					fluxscale = restrict_fluxscale(fluxscale)
# 				
# 					stellar_mass = sources_all[:,47]
# 					stellar_mass = gama_stellarmass_fluxscale(stellar_mass, fluxscale)
# 					stellar_mass_all = stellar_mass
				
				uminusr = sources_all[:,63]
				uminusr_all = uminusr
			
				gminusi = sources_all[:,61]
				gminusi_all = gminusi
				
				sfr = sources_all[:,98]
				sfr = ma.masked_less(sfr,0)
				sfr = log10(sfr)
				sfr_all = sfr
				
# 				sdens = sources_all[:,102]
# 				sdens = ma.masked_less(sdens,0)
# 				sdens = log10(sdens)
# 				sdens_all = sdens
			
			
			
			binprop=eval(bindic[binby][2])    # Which variable to bin by
			binstring=bindic[binby][1]        # xlabel for plots
			binstring_simple=bindic[binby][0] # simple string for other uses
			
			#binprop=colour #log10(Lum_all) #z    
			#binstring=r'$b_{j}-r$' #r"Log(L) (L$_{\odot}$)" # r'$\rm{z}$'    
			#binstring_simple='colour'#'log(L)' #'log(L)' #'z'
			
			# Define the number of bins
			if binning=='y':
				#nbins=int(nbins)		# Note: nbins is defined in raw input
				nbins=int(bindic[binby][4])
			else:
				nbins=1
				
			# Needed for even bin ranges
			if binchoice=='r':
				bin_inc=(max(binprop)-min(binprop))/float(nbins)
			
			# Needed for even number of sources per bin
# 			if binchoice=='N':
# 				sorted=unique(binprop)
# 				bin_inc2=len(sorted)/nbins
					
			ary_right=[min(binprop)]
			
			#********************************
			# Choose elements in required bin
			#********************************
			for binnum in range(0,nbins,1): # remember range doesn't include the end value
			
				sources=sources_all
				#newL=newL_all
				Lum=Lum_all
				z=z_all
				y=y_all
				DL=DL_all
				lum_raw = lum_raw_all
				beamsum = beamsum_all
				if survey=='g':
					stellar_mass = stellar_mass_all
					stellar_mass_raw = stellar_mass_raw_all
					uminusr = uminusr_all
					gminusi = gminusi_all
					sfr = sfr_all
					#sdens = sdens_all
				
				flag+=1
			
				# If using even bin ranges
				if binchoice=='r':
					binA=min(binprop)+binnum*bin_inc
					binB=min(binprop)+(binnum+1)*bin_inc
			
				
				# If using equal numbers of sources per bin
				if binchoice=='N':
					binA = percentile(binprop,binnum*100/nbins)
					binB = percentile(binprop,(binnum+1)*100/nbins)
						
					#binA=sorted[bin_inc2*binnum]
					#binB=sorted[bin_inc2*(binnum+1)-1]
				
				print binstring, 'between',binA,'and',binB
				
				midbin=(binA+binB)/2
				midbin_ary+=[midbin]
				ary_right+=[binB]
				
				if binnum==0:
					bin_msk=where(logical_and(binprop>=binA,binprop<=binB))
				else:
					if binchoice=='r':
						bin_msk=where(logical_and(binprop>=binA,binprop<=binB))
					if binchoice=='N':
						bin_msk=where(logical_and(binprop>binA,binprop<=binB))
					# NOTE: For old binning system, take out the = for binA here!
					
			
				# *********************************	
				# Calculate mean and rms of binned property
				# *********************************		
				binprop_mean=mean(binprop[bin_msk[0]])
				meanbin+=[binprop_mean]
				binprop_rms=std(binprop[bin_msk[0]])
				rmsbin+=[binprop_rms]
				#lum_raw_ary+=[mean(lum_raw[bin_msk[0]])]
				Lum_ary+=[mean(Lum[bin_msk[0]])]
				
				print 'Bin mean and rms =',binprop_mean, binprop_rms
				
				# Also save the mean colour in the bin
				colour_ary+=[mean(colour)]
				
# 				params_lst = [ 'sources', 'Lum', 'DL', 'L_err', 'z', 'zQ' ]
# 				pdict = {  }
# 				for i in params_lst:
# 					pdict[i]=eval(i)	
				
				#print "shape Orig = ", shape(pdict_orig['sources'])
				#print "shape current =", shape(pdict['sources'])
				pdict,y = source_select(bin_msk[0],pdict_orig,y)				
				#print "shape new = ", shape(pdict['sources'])
				for key in pdict:
					exec ("%s = pdict[key]" % key)
				#sources,Lum,DL,y=source_select(bin_msk[0],sources_all,Lum_all,DL_all,y_all)
				
				# gmi etc colour don't need to have column selection here because they are the binprop and are not used anywhere else. stellar_mass can be a binprop, but is also used for gastostar calculations.
				#stellar_mass = stellar_mass[bin_msk[0]]
				
				# Temp replacement. Figure out what is going wrong.
			# 	sources=sources[bin_msk[0]]
			# 	Lum=Lum[bin_msk[0]]
			# 	DL=DL[bin_msk[0]]
			# 	y=y[:,bin_msk[0]]
			
				
				# *** NOTE: Make sure source selection columns ABOVE this point. ***
				#******************************
				# Redefine columns
				#******************************
				if survey=='2df':
					ID,z,zQ,appmag,ra,dec,colour,semimajor,semiminor,sb_bj,sb_r=def_col_2df(sources)
				if survey=='g':
					ID,z,zQ,absmag_ugriz,ra,dec,colour,absmagerr_ugriz, semimajor, semiminor=def_col_gama(sources)
				
				# *********************************
				# Calculate error on luminosity
				# *********************************
				#if survey=='g':
					#lumerr_ugriz=0.921034*absmagerr_ugriz*Lum
				
				
				bin_length=len(ID)
				bin_len_ary+=[bin_length]
				print 'No. in bin = ', bin_length
				print binstring, 'actually between', min(binprop),'and',max(binprop)
				
				# **************************************************************************
				# 							WORK ON SPECTRA								   *
				# 																		   *
				# **************************************************************************
				
				#z_for_Hz = mean(z) # use =0 if H0 wanted in OHI calc, otherwise use mean(z) if H(z) required
				z_for_Hz = 0.0
				
				# ****************************
				# mask default flux values
				# ****************************
				#y_ma=ma.masked_equal(y,-999.)
				y_ma=ma.masked_less(y,-100)
				y_ma=ma.masked_greater(y_ma,100)
				#y_ma=ma.masked_equal(y_ma,0) #for atca
				
				# Scale up the flux due to the Parkes calibration
				if survey=='2df' and pksdat=='p' and fluxcalib=='y': # Should also do this for gama
					zcut=freq_to_z(1310,freq0)
					msk1285=where(z>zcut)[0]
					msk1335=where(z<=zcut)[0]
					y_ma[:,msk1285]=y_ma[:,msk1285]*1.15
					y_ma[:,msk1335]=y_ma[:,msk1335]*1.02
				
				
				# ********************
				# Exclude bad spectra
				# ********************
				if excl=='y':
				
					# Bounds that I have chosen to determine bad spectra
					cliphi=0.03
					cliplo=-0.03
					
					# Find min and max value of each spectrum
					ymins=ma.min(y_ma,axis=0)
					ymaxs=ma.max(y_ma,axis=0)
					
					# Find which columns have 'bad' spectra (ie. outside bounds)
					bad_spect=array(ma.where(logical_or(ymins<=cliplo,ymaxs>=cliphi)))[0]
					good_spect=array(ma.where(logical_and(ymins>cliplo,ymaxs<cliphi)))[0]
					
					# Make a masking template and apply it
					mask_templ=zeros(shape(y_ma),dtype=int8) # make mask template
					mask_templ[:,bad_spect]=1 # Set columns containing bad spectra to flag=1
					y_ma=ma.masked_array(y_ma,mask=mask_templ)
					
					# Make an annotation file of the bad spectra
					ann_ary=make_ann_file(ra[bad_spect],dec[bad_spect])
					#savetxt(savetxt('../improve_spect/badspect_g9.ann',ann_ary,fmt='%4s')
				
				
				# *************************************************
				# Fit and subtract polynomial to each spectra
				# *************************************************
			# 	if fit=='y':
			# 		y_0=ma.filled(y_ma,fill_value=0.0) # Set blank values to 0. Need to do this?
			# 		p=ma.polyfit(x[:,0],y_0,4)
			# 		newy=y-polyval(p,x)
			# 		#newy=cubicfit(x[:,0],y_0)
			# 		y_ma=ma.MaskedArray(newy,y_ma.mask)
				
				# ****************************
				# Blank bad frequency channels
				# ****************************
				
			# For alternative blanking code see older versions of this code
				
				x_ma=ones(shape(y))*x
				
				
				if blank=='y':
					
					if noshift=='y':
						shift=ones(shape(z)) 	#use this for noshift case
					else:
						shift=1+z
						
					shift=ma.reshape(shift,(len(shift),1))
					
					print 'Flagged freqs:'
					for iter_b in bad_i:
						bad_rng=eval('bad'+str(iter_b))
						bad_shift=bad_rng*shift
						print bad_rng
						x_ma=masked_between(x_ma,bad_shift[:,0],bad_shift[:,1])
					
					# This makes sure the x-values are only defined where y!=-999
					#x_ma_yrng=ma.MaskedArray(x_ma,y_ma.mask)
					
					y_ma=ma.MaskedArray(y_ma,x_ma.mask)
					
					# Mask the 1285 band only for the G15 sources in the combo case
					if survey=='g' and field=='combo':
						bad_rng_big = [1252,1310.0]
						g15_msk = ma.masked_less(ra, 211).mask
						shift_ma = ma.MaskedArray(shift, g15_msk)
						shift_ma=ma.filled(shift_ma,0) # Making sure not everything is masked by accident
						bad_shift_big = bad_rng_big * shift_ma
						x_ma_big=masked_between(x_ma,bad_shift_big[:,0],bad_shift_big[:,1])
						y_ma=ma.MaskedArray(y_ma,x_ma_big.mask)
					
				# ***************************************
				# Correct for shifting factor (cosmology)
				# ***************************************
				y_ma = y_ma/(1+z)
					
				
				# *****************************************
				# Plot some individual shifted flux spectra
				# *****************************************
# 				figure()
# 				plot(x,y_ma[:,0:500]*1000)
# 				#plot(x,y_ma*1000)
# 				xlabel(texstr('Frequency\ (MHz)'))
# 				ylabel(texstr('Flux\ density\ (mJy)'))
# 				siglims_plot(siglolim,siguplim,Colour='k')
			# # 	#savefig('../images/sgp_noshift_overlay.gif')
				
				#%figure()
				#%plot_indiv(x,y_ma,'f',siglolim,siguplim)
				#%title(texstr('Unweighted'))	
	
				# *****************************************
				# Calculate the S/N of each individual spectrum
				# *****************************************		
				# Note: currently using the sigreg to define each profile width
				
				#indiv_S,indiv_N,indiv_SN = SN_calc(x,y_ma,siglolim,siguplim,xlolim,xuplim,spec_res,indiv='y',nosig_msk=nosig_msk)
	#From here down
	
				if indivsn=='y':
					print 'After source selection:'
					if survey=='2df':
						linewidth = tf_W(semimajor, semiminor, tf_offset, tf_slope, Babsmag)#absmag_uncorr) # CAUTION: Should be using B-band not bj here. See _sourceselect
					if survey=='g':
						linewidth = ohitf.tf_W(semimajor, semiminor, tf_offset, tf_slope, absmag_uncorr)
					print 'mean(linewidth) = ', mean(linewidth)
					print 'shape(linewidth) = ',shape(linewidth)
					
					# Note: There are some undefined values here which may end up giving the wrong widths. Not even sure if this calculation is correct!
					lwmhz = linewidth[:,0] * freq0 / c
					indiv_siglolim = freq0 - lwmhz/2.0
					indiv_siguplim = freq0 + lwmhz/2.0
					
					#indiv_siglolim = freq0 / (1 + (linewidth/2)/c) 
					#indiv_siguplim = freq0 / (1 - (linewidth/2)/c)
					#lwmhz = indiv_siguplim - indiv_siglolim
					
					figure()
					tmp1=ma.masked_invalid(lwmhz)
					msktf=where(tmp1.mask==False)
					hist(lwmhz[msktf],bins=20,histtype='step',lw=3,ec='k')
					xlabel(texstr('Line\ width\ (MHz)'))

				
				# NOTE: This is only configured for f_or_m=='f'
				
				if indivsn=='y':
					yfit,nosig_msk=basefit(x,y_ma,xlolim,xuplim,siglolim_orig,siguplim_orig,excl_sig='y',multid_y='y',plt='n')
					
					# Note: I am NOT using the fitted spectrum to calc the S/N as it artificially boosts the signal 
					#indiv_S,indiv_N,indiv_SN = SN_calc(x,y_ma,siglolim,siguplim,xlolim,xuplim,spec_res,f_or_m,indiv='y',nosig_msk=nosig_msk)
							
					indiv_S,indiv_N,indiv_SN = calc_indivsn(x_ma, y_ma, indiv_siglolim, indiv_siguplim)
					
					# Run function from OmegaHI_indivSN
					indiv_SN_investig(x,y_ma,yfit,ID,indiv_S,indiv_N,indiv_SN,z,indiv_siglolim, indiv_siguplim, cut=3.7)
				
				
				# ***********************
				# Weighting and scaling
				# ***********************
				
				# Finding the noise level of each spectra
				x_ma2=masked_between(x_ma,siglolim,siguplim)
				x_ma2=masked_between(x_ma2,min(x[:,0]),xlolim)
				x_ma2=masked_between(x_ma2,xuplim,max(x[:,0]))
				y_ma2=ma.MaskedArray(y_ma,x_ma2.mask)
				sigma=ma.std(y_ma2,axis=0)
				
				figure()
				hist(sigma*1000,bins=20)
				xlabel(texstr('\sigma\ (mJy)'))
				ylabel(texstr('Count'))
				
				# Finding the number of sources contributing to stack at each freq
				num_stack=ma.count(y_ma,axis=1)
				num_stack.shape=(len(num_stack),1)
				# y_ma=y_ma/sigma         # Optional direct scaling by sigma
				
				
				# ***********************
				# Choose weight factor!!!!
				# ***********************
				
				#weight=sigma**(-2)/(1+beamsum)
				#weight=array(1.0+0.0*sigma)
				
				#weight=sigma**(-2) # This one is the best!!!
				#weight=sigma**(-2)*Lum[:,2]**2
				#weight=DL**(-4.0)
				
				if survey=='g' and csample=='L':
					wfac=Lum[:,2]
				if survey=='g' and csample=='SM':
					wfac=10*stellar_mass
	
				weight=sigma**(-2.0)*DL**(-1.0)
				
				#weight=sigma**(-2.0)*DL**(-4.0)
				#weight=sigma**(-2.0)*DL**(-4.0)*wfac**2
				
				#weight=(num_stack)**0.5
				
				# Make it 2 dimensional (for use with GAMA)
				bigweight = zeros((len(weight),5))
				for i in range(5):
					bigweight[:,i] = bigweight[:,i] + weight
				
				# ***********************
				# Convert to HI mass axis
				# ***********************
				
				# massy=y_ma*(2.356*10**5)*DL**2 		#This is for y in Jy*km/s
				
				if weighting=='y':
					y_ma_save=y_ma
					y_ma=y_ma*weight
				
				mass_array=y_ma*(4.98*10**7)*DL**2 # Is DL being multiplied correctly???
				
				# if weighting=='y':
				# 	mass_array=mass_array*weight
				
				# Do a baseline fit to each spectrum individually
				if noshift!='y':
					mass_fit=basefit(x,mass_array,xlolim,xuplim,siglolim_orig,siguplim_orig,excl_sig='y',multid_y='y',plt='n')
					mass_fit=mass_fit[0] # Note sure what the extra dimension output is, but only the first one is the array we need
				else:
					mass_fit=mass_array
				
				#mass_fit=mass_array
				
				#%figure()
				#%plot_indiv(x,mass_array,'m',siglolim,siguplim)
				#%title(texstr('Weighted'))
				
				#%figure()
				#%plot_indiv(x,mass_array/weight,'m',siglolim,siguplim)
				#%title(texstr('Unweighted'))
				
				# ******************************
				# Calculate MHI for each spectra
				# ******************************
				
				sigreg=where(logical_and(x<siguplim,x>siglolim))
				MHI=ma.sum(mass_array[sigreg[0],:],axis=0)*spec_res # Array of indiv HI masses.
				# The same, but using the fitted mass spectra
				MHI_fit=ma.sum(mass_fit[sigreg[0],:],axis=0)*spec_res
				
				# NOTE: the DL could be multiplied here instead!!! Check whether [0] needs to be here!!! 
				
				
				#***************************************************************************
				# 							CALCULATING FINAL VALUES					   *
				# 																		   *
				#***************************************************************************
				
				
				#************************
				# Find average luminosity
				#************************
				
				if survey=='2df':
					#Lav=mean(Lum)
					if weighting=='y':
						Lav=ma.sum(Lum*weight)/ma.sum(weight)
					else:
						Lav=mean(Lum) #should be ma.mean?
						
				if survey=='g':
					#Lav=ma.mean(Lum,axis=0)
					if weighting=='y':
						Lav=ma.sum(Lum*bigweight,axis=0)/ma.sum(weight)
						Lav_noweight = mean(Lum,axis=0)
					else:
						Lav=mean(Lum,axis=0) #should be ma.mean?

					#Lav_err=((ma.sum(lumerr_ugriz**2,axis=0))**0.5)/len(lumerr_ugriz[:,2]) #do this in quadrature!
				
				print 'mean Luminosity =',print_sci(Lav,survey) #'%e' %Lav
				#print 'mean Luminosity =',print_sci(Lav,'2df') #'%e' %Lav
				#print 'mean Luminosity error =',print_sci(Lav_err,survey) #'%e' %Lav
				print 'mean(Lum) (no weighting) = ', print_sci(mean(Lum,axis=0),survey)
				#print 'mean(Lum) (no weighting) = ', print_sci(mean(Lum),'2df')
				
				Lav_err = error_on_av(L_err)
				#print "<L> % error = ", '%e' %(100*Lav_err/Lav)
				print "<L> % error = ", print_sci(100*Lav_err/Lav,survey)
				
				#******************
				# Find average mass
				#******************
				
				# Have you applied the same mask to the weights as to the y-vals????
				
				if weighting=='y':
					MHIav=ma.sum(MHI)/sum(weight)
					MHIav_fit=ma.sum(MHI_fit)/sum(weight)
				else:
					MHIav=ma.mean(MHI)
					MHIav_fit=ma.mean(MHI_fit)
				print 'mean MHI[_fit] =','%e' %MHIav,', ','%e' %MHIav_fit
				#print 'median MHI =','%e' %ma.median(MHI)
				
				#**************
				# Find OmegaHI
				#**************
				
				if boost_rhoL=='y':
					rhoL=rhoL*boostrhoL_dct[surveynm] # Defined in ohi_vars
				
				rhoHI=rhoHI_calc(MHIav/Lav,rhoL)
				OmegaHI=rho_to_OmegaHI(rhoHI,z_for_Hz)
				print "OmegaHI =",print_sci(OmegaHI,survey)
				
				rhoHI_fit=rhoHI_calc(MHIav_fit/Lav,rhoL)
				OmegaHI_fit=rho_to_OmegaHI(rhoHI_fit,z_for_Hz)
				print "OmegaHI with fit =",print_sci(OmegaHI_fit,survey)
				
				
				
				print "MASS TO LIGHT METHOD:"
				# Find mass-to-light ratio
				#if survey=='2df':
				if survey=='g':
					MHI.shape=(len(MHI),1)
					MHI_fit.shape=(len(MHI_fit),1)
					
				MassToLight=MHI/Lum
				MassToLight_fit=MHI_fit/Lum
				
				if survey=='g':
					MassToLight=MHI/Lum
					MassToLight_fit=MHI_fit/Lum
				
				
				if weighting=='y':
					Av_MassToLight=ma.sum(MassToLight)/ma.sum(weight) # Note that this only works without sigma weighting!
					Av_MassToLight_fit=ma.sum(MassToLight_fit)/ma.sum(weight)				
				else:
					Av_MassToLight=ma.mean(MassToLight)
					Av_MassToLight_fit=ma.mean(MassToLight_fit)
				#print "Masstolight [without / with fit] =",print_sci(Av_MassToLight,survey),print_sci(Av_MassToLight_fit,survey)
				print "Masstolight [without / with fit] =",'%.3e' %Av_MassToLight, '%.3e' %Av_MassToLight_fit
				
				rhoHI2=rhoHI_calc(Av_MassToLight,rhoL)
				rhoHI2_fit=rhoHI_calc(Av_MassToLight_fit,rhoL)
				
				OmegaHI2=rho_to_OmegaHI(rhoHI2,z_for_Hz)
				OmegaHI2_fit=rho_to_OmegaHI(rhoHI2_fit,z_for_Hz)
				print "OmegaHI from <M/L> [without / with fit] =",print_sci(OmegaHI2,survey),print_sci(OmegaHI2_fit,survey)
				
				
				#***************************************************************************
				# 								STACKING ARRAYS							   *
				#	   																	   *
				#***************************************************************************
			
				# If you want to use the fit to the individual mass spectrum, rather than the fit to the final stack, use indivfit=y and fit=n
				if indivfit=='y':
					mass_array=mass_fit
				
				weightfn_mass=ma.count(mass_array,axis=1)
				weightfn_flux=ma.count(y_ma,axis=1)
				
				if weighting=='y':
					mn_flux=ma.sum(y_ma,axis=1)/sum(weight)
					md_flux=ma.median(y_ma,axis=1)/sum(weight)
					
					mean_stack=ma.sum(mass_array,axis=1)/sum(weight) # #ie. no division by n
					med_stack=ma.median(mass_array,axis=1)/sum(weight)
											
									
					if survey=='g':
					
					#	stellar_mass=ma.masked_less(stellar_mass,7) # Condition should be removed when using proper catalogue
						
						# Masking outliers
						#print 'New mask shape =', shape(where((MHI[:,0]/(10**stellar_mass)/weight)>100)[0])
						#stellar_mass=ma.masked_where((MHI[:,0]/(10**stellar_mass)/weight)>100,stellar_mass)
						
						gastostar_stack=ma.sum(mass_array/10**stellar_mass,axis=1)/sum(weight) 
						Mstar_av_noweight = ma.mean(10**stellar_mass)
						print "<M*> (unweighted) =",'%.3e' %Mstar_av_noweight
						Mstar_av=ma.sum((10**stellar_mass)*weight)/ma.sum(weight)
						
						gastostar_indiv = ( mass_array/10**stellar_mass ) #/weight ?
					
						mass_array.shape=(ma.size(mass_array,axis=0),ma.size(mass_array,axis=1),1)
						#Lum=ma.masked_less(Lum,10**7) # Condition should be removed when using proper catalogue
						Lum.shape=(1,ma.size(Lum,axis=0),ma.size(Lum,axis=1))
					
					ml_stack=ma.sum(mass_array/Lum,axis=1)/sum(weight) # I think this correctly weights (M/L)
					ml_indiv = ( mass_array/Lum ) #/weight?
				
				else:		
					mn_flux=ma.mean(y_ma,axis=1)
					md_flux=ma.median(y_ma,axis=1)	
								
					mean_stack=ma.mean(mass_array,axis=1)
					med_stack=ma.median(mass_array,axis=1)
										
					if survey=='g':
						# Probably for binning purposes?
						#stellar_mass=ma.masked_less(stellar_mass,7) # Condition should be removed when using proper catalogue
						#mass_array = ma.masked_where((mass_array/10**stellar_mass)>50.0, mass_array)
						gastostar_stack=ma.mean(mass_array/10**stellar_mass,axis=1)
						Mstar_av = ma.mean(10**stellar_mass)
						
						gastostar_indiv = mass_array/10**stellar_mass
						
						mass_array.shape=(ma.size(mass_array,axis=0),ma.size(mass_array,axis=1),1)
						#Lum=ma.masked_less(Lum,10**7) # Condition should be removed when using proper catalogue
						Lum.shape=(1,ma.size(Lum,axis=0),ma.size(Lum,axis=1))
					
					ml_stack = ma.mean( mass_array/Lum, axis=1 )
					ml_indiv = mass_array/Lum

				#%figure()
				#%plot_indiv(x,ml_indiv[:,:,mbd[magband]]/weight,'ml',siglolim,siguplim)
				#%title(texstr('Unweighted'))
				
				#%figure()
				#%plot_indiv(x,gastostar_indiv/weight,'mmstar',siglolim,siguplim)
				#%title(texstr('Unweighted'))
					
				DL_ary=(DL**2*ones(shape(y_ma)))
				DL_ma=ma.MaskedArray(DL_ary,y_ma.mask)
				DL_mn=ma.mean(DL_ma,axis=1)
				DL_md=ma.median(DL_ma,axis=1)
				
				if fit=='y':
				
					if genshz=='y':
						excl_sig_a='n'
					else:
						excl_sig_a='y'
					if survey=='g':
						multid_y_a='y'
					else:
						multid_y_a='n'
			
					mean_stack,mskt1=basefit(x,mean_stack,xlolim,xuplim,siglolim_orig,siguplim_orig,excl_sig=excl_sig_a,plt='n')
					med_stack,mskt2=basefit(x,med_stack,xlolim,xuplim,siglolim_orig,siguplim_orig,excl_sig=excl_sig_a)
					mn_flux,mskt1=basefit(x,mn_flux,xlolim,xuplim,siglolim_orig,siguplim_orig,excl_sig=excl_sig_a)
					md_flux,mskt2=basefit(x,md_flux,xlolim,xuplim,siglolim_orig,siguplim_orig,excl_sig=excl_sig_a)
					
					ml_stack,mskt1=basefit(x,ml_stack,xlolim,xuplim,siglolim_orig,siguplim_orig,excl_sig=excl_sig_a, multid_y=multid_y_a)
					
					if survey=='g':
						gastostar_stack,mskt1=basefit(x,gastostar_stack,xlolim,xuplim,siglolim_orig,siguplim_orig,excl_sig=excl_sig_a)
			
			###
			
				#***************************************************************************
				# 									PLOTTING							   *
				#	   																	   *
				#***************************************************************************
			
				if nocontin=='y' and showplots=='y':
					figure()
					plot(ra,dec,'ko')
					a=genfromtxt(contincat)
					plot(a[:,0],a[:,1],'ro')
					xlabel('RA')
					ylabel('Dec')
			
				#*******************************
				# Find appropriate signal range
				#*******************************
				
	# 			if survey=='2df' and pksdat=='h':
	# 				newlim=find_sigrange(x[:,0],mn_flux,20)
	# 				siglolim,siguplim=newlim['xlo'],newlim['xhi']
	# 				print 'New limits',siglolim,siguplim,\
	# 				'at',"%.2f" %newlim['ylo'],"%.2f" %newlim['yhi']
	# 				sigreg=where(logical_and(x<siguplim,x>siglolim))
					
				#***************************
				# To output images for NvN movie
				#***************************			

				if movie=='y':
					figure()
					plot_movie_stack(x,mn_flux,'f',siglolim,siguplim)
					ntitle = str(int(10**(logN/100.)))
					title('N='+ntitle)
					ith_image = str( int((logN-logNloop[0])/float(iterat)) )
					if len(ith_image)==1:
						ith_image='0'+ith_image				
					savefig('../movie/'+ith_image+'.png')

				#***************************
				# Plot flux and mass stacks
				#***************************
				
				if showplots=='y':
					figure()
					plot_stack(x,mn_flux,'f',siglolim,siguplim)
					#plot_stack_2xax(x,mn_flux,'f',siglolim,siguplim)
					
					figure()
					plot_stack(x,mean_stack,'m',siglolim,siguplim)
					#subplot(322)
					#subplot(122)

					figure()
					plot_indiv(x,ml_stack,'ml',siglolim,siguplim)
					figure()
					plot_stack(x,ml_stack[:,2],'ml',siglolim,siguplim)
					
					if survey=='g':
						figure()
						plot_stack(x, gastostar_stack, 'mmstar', siglolim, siguplim)
		
				#***********************************
				# Plot overlay with shuffled cat
				#***********************************
				
				if plotoverlay=='y':
					if avtyp=='mn':
						if survey=='2df':
							arrdict = {'f':mn_flux, 'm':mean_stack, 'ml':ml_stack}
						if survey=='g':
							arrdict = {'f':mn_flux, 'm':mean_stack, 'mmstar':gastostar_stack, 'ml':ml_stack[:,mbd[magband]]}
						
						inarr = arrdict[plotpar] 
	
					figure()
					stack_overlay(binnum,x,inarr,choose_rndary(binnum,plotpar),plotpar)
					if binning=='y':
						title(r'$'+'%.3f' %binA+' < '+binstring_simple+' < '+'%.3f' %binB+'$')
					ylim(-3,5)
						
					if saveplot=='y':
						#eval(fmdic[f_or_m][1])
						savefig(saveplot_name(binnum))
	# 			
	# 			#**********************
	# 			# Plot noise components
	# 			#**********************
	# 			# This one is the N per channel normal plot
				#%figure() # Plot in new figure window
				#subplot(221)
				#%plot_Nweight(x,weightfn_flux,'N')
	# 			
	# 			#figure()
	# 			subplot(222)
	# 			noisefac_mn=DL_mn/(weightfn_flux)**0.5
	# 			noisefac_md=DL_md/(weightfn_flux)**0.5
	# 			plot_noise2(x,noisefac_mn,noisefac_md,r'<DL$^{2}$>/N$^{0.5}$','r','k',r'mean(DL$^2$)',r'median(DL$^2$)')
	# 			
	# 			#figure()
	# 			subplot(223)
	# 			noisefac=1./weightfn_flux**0.5
	# 			plot_noise(x,noisefac,r'1/N$^{0.5}$','g')
	# 			
	# 			#figure()
	# 			subplot(224)
	# 			noisefac1=DL_mn
	# 			noisefac2=DL_md
	# 			plot_noise2(x,noisefac1,noisefac2,r'<DL$^{2}$>','b','k',r'mean(DL$^2$)',r'median(DL$^2$)')
	# 			
				#**************************************************
				# Plot stacked flux and mass scaled by noise factor
				#**************************************************
				
			# 	figure()
			# 	subplot(121)
			# 	flux_stack(x,mn_flux/noisefac,md_flux/noisefac)
			# 	title('Flux scaled by noise factor')
			# 	
			# 	#figure()
			# 	subplot(122)
			# 	mass_stack(x,mean_stack/noisefac_mn,med_stack/noisefac_md)
			# 	title('Mass scaled by noise factor')
			
				
				#***************************************************************************
				# 							AVERAGE THEN INTEGRATE METHOD				   *
				#	   																	   *
				#***************************************************************************
				
				
				# Set the signal boundaries based on hicat
				if confuse=='y' and cscale=='y': 
					siglolimC,siguplimC=vdict[dv]['f_min'],vdict[dv]['f_max']

				# Do any requested profile fitting
				if fit_sin=='y':
					inp={'mn':mean_stack,'md':med_stack}
					#vars()[inp[avtyp]]=remove_sinusoid(x,inp[avtyp],'m')
					mean_stack=remove_sinusoid(x[:,0],mean_stack,'m')
				if fit_gaussian=='y': # Note this does NOT replace the spectrum
					gauss_mass,[ga,gm,gs]=fit_gauss(x[:,0],mean_stack,'m')
					#fit_gauss(x[:,0],mn_flux,'f')
		
					# Correct for profile offset from freq0
					if confuse=='y' and cscale=='y': 
						centre_offset=freq0-gm 
						siglolimC-=centre_offset
						siguplimC-=centre_offset
					
						siglims_plot(siglolimC,siguplimC)
			
				# Find the integrated flux within the signal reg
				Sint_mn=integrate_spectrum(mn_flux,sigreg[0],spec_res)
				Sint_md=integrate_spectrum(md_flux,sigreg[0],spec_res)
				
				print 'Integrated (mean) flux = ', Sint_mn, 'Jy MHz'
				
				# Define the signal region over which to find the integrated mass
				if confuse=='y' and cscale=='y':		
					s_bounds=x_mask(x,siglolimC,siguplimC,None,None,excl_sig='n')
					
				else:
					s_bounds=sigreg[0]
				
				MHI_stack_av=integrate_spectrum(mean_stack,s_bounds,spec_res) 
				MHI_stack_md=integrate_spectrum(med_stack,s_bounds,spec_res)
				
				ml_stack_av=integrate_spectrum(ml_stack,s_bounds,spec_res) 
				
				if survey=='g':
					gastostar_stack_av=integrate_spectrum(gastostar_stack,s_bounds,spec_res)
					print  "Stacked <MHI/M*> =",'%.3e' %gastostar_stack_av
				
				# Scale up the mass by the 'confusion fraction' if required
				if confuse=='y' and cscale=='y':
					MHI_stack_av=MHI_stack_av/(vdict[dv]['frac']/100.)
					MHI_stack_md=MHI_stack_md/(vdict[dv]['frac']/100.)
				
				ary_MHI+=[MHI_stack_av]
				#if weighting=='y':
					#MHI_stack_av=MHI_stack_av/sum(weight)
				print "Stacked mean mass =",'%.3e' %MHI_stack_av
				#print "Stacked median mass =",'%.3e' %MHI_stack_md
				print "Stacked <M/L> =", print_sci(ml_stack_av,survey)
				
				

				
				if f_factor=='y':
					if confuse=='y':
						Ff = Ffd[surveynm] # Read from OmegaHI_vars
						weight_avMavL = Hfd[surveynm]
					else:
						#Ff = ohi_f.Ffactor(Lum)
						Ff = ohi_f.Ffactor_weight(weight,Lum[0,:,2], alpha_ml, phi_rhoL, alpha_rhoL, Lstar_rhoL, typ='avMX') # Calculate it from OHI_Ffactor.py. Also known as weight_avML
						weight_avMavL = ohi_f.Ffactor_weight(weight,Lum[0,:,2], alpha_ml, phi_rhoL, alpha_rhoL, Lstar_rhoL,typ='avMavX') # The _ml and _rhoL values come from OHI_vars
						#weight_avMavL = ohi_f.Hfactor(Lum[0,:,2], alpha_ml, phi_rhoL, alpha_rhoL, Lstar_rhoL) 
						
					print 'Ff = ',Ff
					print 'Hf = ', weight_avMavL
	
	
				# OmegaHI from median
				rhoHI_stack_md=rhoHI_calc(MHI_stack_md/Lav,rhoL)
				OmegaHI_stack_md=rho_to_OmegaHI(rhoHI_stack_md,z_for_Hz)
				#OmegaHI_stack_md=OmegaHI_calc(MHI_stack_md/Lav,rhoL)
				#print "*median* Stacked OmegaHI =",print_sci(OmegaHI_stack_md,survey)				

				# Omega HI from <M>/<L>			
				rhoHI_stack=rhoHI_calc(MHI_stack_av/Lav,rhoL)
				OmegaHI_stack=rho_to_OmegaHI(rhoHI_stack,z_for_Hz)				
				print '<M>/<L> = ', print_sci(MHI_stack_av/Lav,survey)
				print "Stacked OmegaHI =",print_sci(OmegaHI_stack,survey)
				if f_factor=='y':	
					OmegaHI_stack, rhoHI_stack = OmegaHI_stack*weight_avMavL, rhoHI_stack*weight_avMavL			
					print "Stacked OmegaHI from <M>/<L> *H =",print_sci(OmegaHI_stack,survey)

				# Omega HI from <M>/<L>	(no <L> weight)		
				rhoHI_stack2=rhoHI_calc(MHI_stack_av/Lav_noweight,rhoL)
				OmegaHI_stack2=rho_to_OmegaHI(rhoHI_stack2,z_for_Hz)				
				print '<M>/<L> (no weight)= ', print_sci(MHI_stack_av/Lav_noweight,survey)
				print "Stacked OmegaHI (no weight) =",print_sci(OmegaHI_stack2,survey)
				if f_factor=='y':	
					OmegaHI_stack2, rhoHI_stack2 = OmegaHI_stack2*weight_avMavL, rhoHI_stack2*weight_avMavL			
					print "Stacked OmegaHI from <M>/<L> *H  (no <L> weight)=",print_sci(OmegaHI_stack2,survey)


				# Calculate OHI from the <M/L> method
				rhoHI_ml_stack=rhoHI_calc(ml_stack_av,rhoL)
				OmegaHI_ml_stack=rho_to_OmegaHI(rhoHI_ml_stack,z_for_Hz)				
				#OmegaHI_stack=OmegaHI_calc(MHI_stack_av/Lav,rhoL)
				print "Stacked OmegaHI from <M/L> =",print_sci(OmegaHI_ml_stack,survey)
				if f_factor=='y':
					OmegaHI_ml_stack=OmegaHI_ml_stack*Ff
					print "Stacked OmegaHI from <M/L>*Ff =",print_sci(OmegaHI_ml_stack,survey)

				### Omega HI and F-factors with stellar masses ###

				if survey=='g':
				
					if f_factor=='y':
						if confuse=='y':
							weight_avMMstar = Ffmstard[surveynm] # Read from OmegaHI_vars
							weight_avMavMstar = Hfmstard[surveynm]
						else:
							#Ff = ohi_f.Ffactor(Lum)
							weight_avMMstar = ohi_f.Ffactor_sm_weight(weight, 10**stellar_mass, alpha_mmstar, phi1_rhoMstar, phi2_rhoMstar, alpha1_rhoMstar, alpha2_rhoMstar, xstar_rhoMstar, typ='avMX') # Calculate it from OHI_Ffactor.py. Also known as weight_avML
							weight_avMavMstar = ohi_f.Ffactor_sm_weight(weight, 10**stellar_mass, alpha_mmstar, phi1_rhoMstar, phi2_rhoMstar, alpha1_rhoMstar, alpha2_rhoMstar, xstar_rhoMstar, typ='avMavX')
							#weight_avMavMstar = ohi_f.Hfactor(10**stellar_mass, alpha_mmstar, phi_rhoMstar, alpha_rhoMstar, Mstarstar_rhoMstar) # The _mmstar and _rhoMstar values come from OHI_vars
						
						print 'Ff_mstar = ',weight_avMMstar
						print 'Hf_mstar = ', weight_avMavMstar


					# Calculate OHI from the <M>/<M*> method
					rhoHI_mmstar_stack1 = rhoHI_calc(MHI_stack_av/Mstar_av,rhoMstar)
					OmegaHI_mmstar_stack1 = rho_to_OmegaHI(rhoHI_mmstar_stack1,z_for_Hz)
					print "<M>/<M*> =",'%.3e' %(MHI_stack_av/Mstar_av)
					print "Stacked OmegaHI from <M>/<M*> =",'%.3e' %OmegaHI_mmstar_stack1
					if f_factor=='y':
						OmegaHI_mmstar_stack1 = OmegaHI_mmstar_stack1 * weight_avMavMstar
						print "Stacked OmegaHI from <M>/<M*>*Hf =",'%.3e' %(OmegaHI_mmstar_stack1)

					# Calculate OHI from the <M>/<M*> method (no <M*> weight)
					rhoHI_mmstar_stack1nw = rhoHI_calc(MHI_stack_av/Mstar_av_noweight,rhoMstar)
					OmegaHI_mmstar_stack1nw = rho_to_OmegaHI(rhoHI_mmstar_stack1nw,z_for_Hz)
					print "<M>/<M*> (no weight) =",'%.3e' %(MHI_stack_av/Mstar_av_noweight)
					print "Stacked OmegaHI from <M>/<M*> (no weight) =",'%.3e' %OmegaHI_mmstar_stack1nw
					if f_factor=='y':
						OmegaHI_mmstar_stack1nw = OmegaHI_mmstar_stack1nw * weight_avMavMstar
						print "Stacked OmegaHI from <M>/<M*>*Hf (no <M*> weight) =",'%.3e' %(OmegaHI_mmstar_stack1nw)

					# Calculate OHI from the <M/M*> method
					rhoHI_mmstar_stack = rhoHI_calc(gastostar_stack_av,rhoMstar)
					OmegaHI_mmstar_stack = rho_to_OmegaHI(rhoHI_mmstar_stack,z_for_Hz)	
					print "Stacked OmegaHI from <M/M*> =",'%.3e' %OmegaHI_mmstar_stack
					if f_factor=='y':
						OmegaHI_mmstar_stack = OmegaHI_mmstar_stack * weight_avMMstar
						print "Stacked OmegaHI from <M/M*>*Ff =",'%.3e' %(OmegaHI_mmstar_stack)								
				

				
				print 'Max no. stacked = ', max(weightfn_mass)
				print 'No. stacked at f0 = ', weightfn_mass[find_nearest(x,freq0)[0]]
	
				
				if avtyp=='mn':
				
					if survey=='2df':
						#if f_factor=='y':
						#	ary_OHI=ary_OHI+[OmegaHI_ml_stack*Ff]
						#else:
						ary_OHI=ary_OHI+[OmegaHI_ml_stack]
				
					if survey=='g':
						ary_OHI=ary_OHI+[OmegaHI_stack[mbd[magband]]] #choosing just r-band FIX THIS!
						
				if avtyp=='md':
					
					if survey=='2df':
						ary_OHI=ary_OHI+[OmegaHI_stack_md]
				
					if survey=='g':
						ary_OHI=ary_OHI+[OmegaHI_stack_md[mbd[magband]]] #choosing just r-band FIX THIS!
						

				
				# Calculate the rms in the stacked spectrum excluding the signal region
				#sigmsk=mask_sigreg(x,siglolim,siguplim)
				if fit=='y':
					[xlo,xhi]=[xlolim,xuplim]
				else:
					[xlo,xhi]=[min(x)[0],max(x)[0]]
				sigmsk=x_mask(x,xlolim,xuplim,siglolim,siguplim,excl_sig='y')
				rms_all=std(mean_stack[sigmsk])
				print 'rms excl signal region = ','%.2e' %rms_all
				
				print ' '
				print 'mean z = ',mean(z)
				print ' '
				
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
				
				#***********************************
				# Save a spectrum array (from utils)
				#***********************************
				if savearr=='y':
				
				# 	save_spect(x,md_flux,'../ra_offset/md_flux_blank.txt')
				# 	save_spect(x,mn_flux,'../ra_offset/mn_flux_blank.txt')
					outspect_flux, outspect_mass, outspect_ml, outspect_mmstar = choose_outspect(cat_num,binnum,logN)
	
# Uncomment this if you want to output *both* the mean and md spect for *either* flux or mass
# 					if f_or_m=='m':
# 						save_spect(x,mean_stack,outspect_mn)
# 						#save_spect(x,med_stack,outspect_md) #Not interested in median at the moment
# 					if f_or_m=='f':
# 						save_spect(x,mn_flux,outspect_mn)
# 						#save_spect(x,md_flux,outspect_md)

					save_spect(x,mn_flux,outspect_flux)
					save_spect(x,mean_stack,outspect_mass)
					if use_ml=='y':
						if survey=='2df':
							save_spect(x,ml_stack,outspect_ml)
						if survey=='g':
							save_spect(x,ml_stack[:,mbd[magband]],outspect_ml)
							save_spect(x,gastostar_stack,outspect_mmstar)
							
				#***********************************
				# Find the S/N
				#***********************************
					
				if calcsn=='y':
		
					# For mass only!
					if f_or_m=='m':
						if avtyp=='mn':
							y_ax=mean_stack
						if avtyp=='md':
							y_ax=med_stack
					
					# For flux only!
					if f_or_m=='f':
						if avtyp=='mn':
							y_ax=mn_flux
						if avtyp=='md':
							y_ax=md_flux
					
					# This finds the integrated and peak mass and flux errors
					print 'For mass/flux:'
					integ_sig,integ_noise,SN_integ,peak_flux,std_noise,std_noise_err,SN_peak = SN_calc(x,y_ax,siglolim,siguplim,xlolim,xuplim,spec_res,f_or_m,shuf_dir=choose_shufdir(binnum,logN,f_or_m,ml='n'))
					
					# This finds the integrated <M/L> error
					if use_ml=='y':
						print 'For M/L:'
						ML_err = SN_calc(x,ml_stack[:,2],siglolim,siguplim,xlolim,xuplim,spec_res,'ml',shuf_dir=choose_shufdir(binnum,logN,f_or_m,ml='y',spectyp='ml'))[1]
						if survey=='g':
							print 'For M/M*:'
							gastostar_err = SN_calc(x,gastostar_stack,siglolim,siguplim,xlolim,xuplim,spec_res,'mmstar',shuf_dir=choose_shufdir(binnum,logN,f_or_m,ml='n',spectyp='mmstar'))[1]
					
					NvN_mean+=[std_noise]
					NvN_std+=[std_noise_err]
					
					integ_noise_ary+=[integ_noise]
					peaksn_ary+=[SN_peak]
					
				#***********************************
				# Find the error on OHI and rhoHI and the mass-to-light ratio
				#***********************************
				
				#mean_Mstar = mean(10**stellar_mass) #Mstar_av already defined above with proper weighting
				print '<M*> = ', '%.3e' %mean(10**stellar_mass)
				
				if use_ml=='y':				
					ary_ML+=[ml_stack_av] # This is using <M/L>
					print "<M/L> = ", print_sci(ary_ML[0],survey)					
				else:
					# This is the original
					ary_ML+=[MHI_stack_av/Lav]
					print "<M>/<L> = ", print_sci(ary_ML[0],survey)	
					
				if survey=='g':
					gastostar_ary+=[gastostar_stack_av]	
				
				if calcsn=='y':
		
# 					if survey=='2df'
# 						# This is for the <M/L> method
# 						if use_ml=='y':
# 						
# 							#ML_err = integ_noise # this is now done above
# 							print "<M/L> error = ", print_sci(ML_err,'2df')
# 						
# 							OHI_err, rhoHI_err = OHI_error_v2(ml_stack_av, ML_err, rhoL, rhoL_err, z=z_for_Hz)
# 							print "OHI error = ", print_sci(OHI_err,survey)
# 					
# 						# This is for the <M>/<L> method
# 						else:
# 							OHI_err1,rhoHI_err1 = OHI_error(OmegaHI_stack,MHI_stack_av,integ_noise,rhoL,rhoL_err,Lav,z_for_Hz,Lav_err=0)
# 							print "OHI error = ", print_sci(OHI_err1,survey)
# 						
# 							OHI_err,rhoHI_err = OHI_error(OmegaHI_stack,MHI_stack_av,integ_noise,rhoL,rhoL_err,Lav,z_for_Hz,Lav_err)
# 							print "OHI error incl Lav_err = ", print_sci(OHI_err,survey)
# 						
# 							print "% diff = ", 100*(OHI_err-OHI_err1)/OHI_err
# 						
# 							dummy, ML_err = OHI_error(1,MHI_stack_av,integ_noise,1,0,Lav,z_for_Hz,Lav_err)
# 							print "<M>/<L> error = ", print_sci(ML_err,survey)

					# OHI error from <M/L>
					print "<M/L> error = ", print_sci(ML_err,'2df')					
					OHI_errL, rhoHI_errL = OHI_error_v2(ml_stack_av, ML_err, rhoL, rhoL_err, z=z_for_Hz)
					print "OHI error = ", print_sci(OHI_errL,survey)
					
					# OHI error from <M/SM>
					print '<M/SM> = ', print_sci(gastostar_stack_av,'2df')	
					print "<M/SM> error = ", print_sci(gastostar_err,'2df')		
					OHI_errSM, rhoHI_errSM = OHI_error_v2(gastostar_stack_av, gastostar_err, rhoMstar, rhoMstar_err, z=z_for_Hz)
					print "OHI error = ", print_sci(OHI_errSM,'2df')
					
					# OHI error from <M>/<L>
					print '<M>/<L> = ',MHI_stack_av/Lav
					print '<M>/<L> (unweighted) = ',MHI_stack_av/Lav_noweight
					dummy, ML_err1 = OHI_error(1,MHI_stack_av,integ_noise,1,0,Lav,z_for_Hz,Lav_err)
					print "<M>/<L> error = ", print_sci(ML_err1,survey)

					OHI_err1L,rhoHI_err1L = OHI_error(OmegaHI_stack,MHI_stack_av,integ_noise,rhoL,rhoL_err,Lav,z_for_Hz,Lav_err=0)
					print "OHI error = ", print_sci(OHI_err1L,survey)				
					OHI_err1Lnw,rhoHI_err1Lnw = OHI_error(OmegaHI_stack,MHI_stack_av,integ_noise,rhoL,rhoL_err,Lav_noweight,z_for_Hz,Lav_err=0)
					print "OHI error (no weight) = ", print_sci(OHI_err1Lnw,survey)
				
# 					OHI_err,rhoHI_err = OHI_error(OmegaHI_stack,MHI_stack_av,integ_noise,rhoL,rhoL_err,Lav,z_for_Hz,Lav_err)
# 					print "OHI error incl Lav_err = ", print_sci(OHI_err,survey) 				
# 					print "% diff = ", 100*(OHI_err-OHI_err1)/OHI_err

					# OHI error from <M>/<SM>
					print '<M>/<SM> = ',MHI_stack_av/Mstar_av
					print '<M>/<SM> (unweighted) = ',MHI_stack_av/Mstar_av_noweight
					dummy, Mmstar_err1 = OHI_error(1,MHI_stack_av,integ_noise,1,0,Mstar_av,z_for_Hz,Lav_err=0)
					print "<M>/<SM> error = ", print_sci(Mmstar_err1,'2df')

					OHI_err1SM,rhoHI_err1SM = OHI_error(1,MHI_stack_av,integ_noise,rhoMstar,rhoMstar_err,Mstar_av,z_for_Hz,Lav_err=0)
					print "OHI error = ", print_sci(OHI_err1SM,'2df')
					# Without Mstar weighting
					OHI_err1SMnw,rhoHI_err1SMnw = OHI_error(1,MHI_stack_av,integ_noise,rhoMstar,rhoMstar_err,Mstar_av_noweight,z_for_Hz,Lav_err=0)
					print "OHI error (no weight) = ", print_sci(OHI_err1SMnw,'2df')
					
					OHI_err,rhoHI_err = OHI_errL,rhoHI_errL # temporary only!!!
						
					if f_factor=='y':
						OHI_err,rhoHI_err = OHI_err*Ff ,rhoHI_err*Ff
						print "OHI error incl F = ", print_sci(OHI_err,survey)
				
					
					OHI_err_ary+=[OHI_err]
					OHI_err_ary_conf+=[OHI_err]
					ML_err_ary+=[ML_err]
					if survey=='g':
						gastostar_err_ary+=[gastostar_err]
				
				#Space between output in different bins	
				print " "

				OHI_ary_conf+=[ary_OHI]
				
				# Various L_err plots (functions defined in OmegaHI_plots.py)
# 				figure()
# 				subplot(224)
# 				plot_z_v_zerr(z, zerr)
# 				subplot(222)
# 				plot_L_v_Ler(L_err, Lum)
# 				subplot(221)
# 				Ler_hist(L_err, Lum)
# 				subplot(223)
# 				plot_z_v_Ler(L_err, Lum, z)
			
			#**********************************	
			# Plot histogram of binned sources
			#**********************************
			figure()
			#subplot(224)
			#grid()
			for i in range(0,nbins):
				# Should change this to read in the array of binned values, not re-calculate them here.
				#bin_count=binprop[where(logical_and(binprop>ary_right[i],binprop<=(ary_right[i]+bin_inc)))]
				bin_count=binprop[where(logical_and(binprop>=ary_right[i],binprop<=(ary_right[i+1])))]
				# bin_count now defined at start of binning loop
				#hist(bin_count,range=[min(binprop),max(binprop)],bins=20,label=str(len(bin_count)))#,fc='k',ec='w')
				hist(bin_count,range=[min(binprop),max(binprop)],bins=20,label=str(bin_len_ary[i]))#,fc='k',ec='w')
								
			#if binning=='y':
				#vlines(boundaries,min(ylim()),max(ylim()),linestyles='dotted',lw='1.3',color='k')
				
			legend(loc=0,numpoints=1,frameon=False)
			xlabel(binstring)
			ylabel(r'$\rm{Count}$')
			
			# Plot the absolute magnitude distribution
			if binby=='L':
				figure()
				for i in range(0,nbins):
					bin_count=binprop[where(logical_and(binprop>=ary_right[i],binprop<=(ary_right[i+1])))]
					if survey=='g':
						Mag_count=lum_to_absmag(10**bin_count,SolarMag_ugriz[mbd[magband]])
					if survey=='2df':
						Mag_count=lum_to_absmag(10**bin_count,SolarMag_b)
					hist(Mag_count,range=[-26,-8],bins=20,label=str(bin_len_ary[i]))
				xlabel(r'$M_r$')
				ylabel(texstr('Count'))
				legend(loc=0,numpoints=1,frameon=False)

			#**************************************
			# Plot OmegaHI/MHI variation over the bins
			#**************************************


			if binning=='y':	

				if savebin=='y':
					# Save bin properties and bin ranges
					
					if survey=='2df':
						bininfoary=array((meanbin,ary_MHI,integ_noise_ary,ary_OHI,OHI_err_ary, ary_ML, ML_err_ary,rmsbin,bin_len_ary,peaksn_ary,Lum_ary))
					if survey=='g':
						bininfoary=array((meanbin,ary_MHI,integ_noise_ary,ary_OHI,array(OHI_err_ary)[:,mbd[magband]], array(ary_ML)[:,mbd[magband]], ML_err_ary,rmsbin,bin_len_ary,peaksn_ary,gastostar_ary,gastostar_err_ary,Lum_ary))
					#print '*** Saving temporary binning array!!! ***'
					#bininfoary=array((meanbin,ary_MHI,ary_OHI,rmsbin))
					
					# binsavenm function defined in OmegaHI_control.py
					bininfonm=binsavenm('bininfo_',binby,survey)
					binboundsnm=binsavenm('binbounds_',binby,survey)
					
					bincolournm=binsavenm('bincolour_',binby,survey)		
					
					savetxt(bininfonm,bininfoary.T)
					savetxt(binboundsnm,ary_right)
					savetxt(bincolournm,colour_ary)
				
				
				# Plot the MHI and OHI vs binprop

				ary_OHI=array(ary_OHI)
 				ary_MHI=array(ary_MHI)	
 				ary_ML=array(ary_ML)
 				
 				if survey=='g':
					gastostar_ary=array(gastostar_ary)

				if calcsn!='y':
					figure()
					binplotrun(meanbin,ary_MHI,None,None,0,0,binstring,r'$M_{\rm{HI}}$',ary_right)

					if survey=='2df':
						figure()
						binplotrun(meanbin,ary_ML,None,None,0,0,binstring,r'$M_{\rm{HI}}/L$',ary_right)
					
					if survey=='g':				
						figure()
						binplotrun(meanbin,log10(ary_ML[:,2]),None,None,0,0,binstring,r'$M_{\rm{HI}}/L$',ary_right)				
						# Overplotting the Toribio line
						xt=linspace(9.2,11)
						alphat = -0.4
						logbetat = 3.734
						logmlt=logbetat+(alphat)*xt
						plot(xt,logmlt,'g:',lw=3,label=texstr('Toribio'))
						
						figure()
						binplotrun(meanbin,log10(gastostar_ary),None,None,0,0,binstring,r'$M_{\rm{HI}}/M*$',ary_right)
						# Overplotting the Huang (2012) line
						xh=linspace(7.5,12)
						alphah = -0.724
						logbetah = 7.266
						logmlh=logbetah+(alphah)*xh
						alphah2 = -0.288
						logbetah2 = 3.206
						logmlh2=logbetah2+(alphah2)*xh

						plot(xh,logmlh,'b:',lw=3,label=texstr('Huang'))
						plot(xh,logmlh2,'m:',lw=3,label=texstr('Huang'))

 				
 				if calcsn=='y':
					ML_err_ary=array(ML_err_ary)
					OHI_err_ary=array(OHI_err_ary)			
				
# 					figure()
# 					binplotrun(meanbin,ary_MHI,None,integ_noise_ary,nobinM,nobinM_er,binstring,r'$M_{\rm{HI}}$',ary_right)
# 					
# 					figure()
# 					binplotrun(meanbin,ary_OHI*1000,None,OHI_err_ary*1000,nobinO*1000,nobinO_er*1000,binstring,r'$\Omega_{HI} \times 10^{-3}$',ary_right)
		
		
		# Plot the redshift distribution and the full 2dfgrs Nz (scaled to the same field size). Can also be done by running Nz_plot.py
# 		if survey=='2df' and pksdat=='p':
# 			ohi_other.nz_n_full2df(z)
		
	#*************************************************************
	# Plot variation of OHI with Delta V (confused velocity width)
	#*************************************************************
	# (might not currently work without calcsn=='y')
	# Commented out when using only one dV value.
# 	if confuse=='y':
# 		figure()
# 		OHI_ary_conf=ravel(array(OHI_ary_conf))
# 		OHI_err_ary_conf=ravel(array(OHI_err_ary_conf))
# 
# 		#bin_plot2(dv_ary,OHI_ary_conf*1000,None,OHI_err_ary_conf,0,'OmegaHI (x10**-3)','Confused velocity')
# 		
# 		yval=OHI_ary_conf*1000
# 		#if f_factor=='y': #now taken care of above where OHI_ary is defined
# 		#	yval=yval*Ff
# 		if calcsn=='y':
# 			yerrval=OHI_err_ary_conf
# 		else:
# 			yerrval=None
# 		
# 		
# 		bin_plot2(dv_ary,yval,None,yerrval,0,'OmegaHI (x10**-3)','Confused velocity')
# 		
# 		# Overlay non-corrected OHI estimate
# 		# 0.61
# 		hlines(1.18,0,900,linestyle='dashed',color='red',label='No confusion correction')
# 		hlines(1.18+0.05,0,900,linestyle='dotted',color='r',label='Error margin')
# 		hlines(1.18-0.05,0,900,linestyle='dotted',color='r')
# 		
# 		# Overlay HIPASS OHI estimate
# 		hlines(0.26,0,900,linestyle='dashed',color='g',label='Zwaan et al. 2005')
# 		hlines(0.26+0.03,0,900,linestyle='dotted',color='g')
# 		hlines(0.26-0.03,0,900,linestyle='dotted',color='g')
# 		
# 		#legend(loc=5)
# 		xlabel(r'$\Delta V\ (\rm{kms}^{-1})$',size='large')
# 		ylabel(r'$\Omega_{HI}\times10^{3}$',size='large')
# 		#xlim(0,710)
# 		
# 		# Save the values in the plot for later use
# 		outary_c=array([dv_ary,OHI_ary_conf,OHI_err_ary_conf]).T
		#savetxt('../catalogues/sgp_ohi_vs_v_normal.txt',outary_c)
		
		# To load this in and make the appropriate plot, use:
		#dv_ary,OHI_ary_conf,OHI_err_ary_conf=genfromtxt('../catalogues/sgp_ohi_vs_v_test1.txt').T
		# Then copy and paste the commands above, from bin_plot2...
			
			
	#show()	

# The linewidth distributions and plots AFTER source selection. Use this to save the lw arrays and use lw_compare.py to plot them against eachother
print 'After source selection:'
#linewidth = tf_W(semimajor, semiminor, tf_offset, tf_slope, Babsmag)#absmag_uncorr) # CAUTION: Should be using B-band not bj here. See _sourceselect
linewidth = tf_W(semimajor, semiminor, tf_offset, tf_slope, absmag_uncorr)[:,2]
print 'mean(linewidth) = ', mean(linewidth)
print 'shape(linewidth) = ',shape(linewidth)
linewidth_mhz = dv_to_df(linewidth)
savetxt('../ghipass/lw_ghipass2.txt', linewidth)
savetxt('../ghipass/lw_mhz_ghipass2.txt', linewidth_mhz)

# Examining the completeness of the sample


zplot = linspace(min(z),max(z))
zplot_msk=where(zplot>0.014)[0]

Llim = lum_completeness_lim(z, appmag_lim, SolarMag_ugriz[mbd[magband]])
Lr = Lum[0,:,2]
Lr_raw = lum_raw[:,2]
Lmsk = where(Lr>Llim)[0]
Lplot = lum_completeness_lim(zplot, appmag_lim, SolarMag_ugriz[mbd[magband]])
#savetxt('../ghipass/Lcomplete_IDs.txt',ID[Lmsk])

# This was the 'by eye' method that I used
SMlim_me = sm_completeness_lim(DL, SM_flux_lim)
# Using Aaron's SM completeness cut that accounts for colour bias
SMlim = utils_gama.SMcompleteness_AR(z)
SMplot = utils_gama.SMcompleteness_AR(zplot)

sm_msk = where(10**stellar_mass>SMlim)[0]
LnSMmsk = where(logical_and(Lr>Llim, 10**stellar_mass>SMlim))[0]
#savetxt('../ghipass/L+SMar_complete_IDs.txt',ID[LnSMmsk])


#%fig=figure()
#%ax1=subplot(111)
#plot(z[Lmsk], Lr[Lmsk],'k.')
#%sc = scatter(z,Lr,c=gminusi,s=20,hold=True,marker='o',edgecolor='None',vmin=0.2,vmax=1.2)#,cmap='autumn')
#%plot(zplot, Lplot, 'k--',lw=3)
#%xlabel(r'$z$')
#%ylabel(r'$L_r$')
#%yscale('log')
#%cax=colorbar(sc)
#%cax.ax.set_ylabel(r'$g-i$')


#%fig=figure()
#%ax1=subplot(111)
#plot(z[sm_msk], stellar_mass[sm_msk],'k.')
#%sc = scatter(z,stellar_mass,c=gminusi,s=15,hold=True,marker='o',edgecolor='None',vmin=0.2,vmax=1.2)
#plot(z, log10(SMlim_me), 'k.')
#%plot(zplot[zplot_msk], log10(SMplot)[zplot_msk], 'k--',lw=3)
#%if field=='hipass':
#%	hlines(8.025-2*log10(100/70.),0.003,0.0140,color='red',linestyle='--',lw=3)

#%xlabel(r'$z$')
#%ylabel(r'$\log(\mathcal{M})$')
#%cax=colorbar(sc)
#%cax.ax.set_ylabel(r'$g-i$')

boost_fact = utils_gama.completenessF_gama_r(z, appmag_lim, -24, SolarMag_ugriz[mbd[magband]])
SMboost_fact, tmp = utils_gama.SMcompletenessF_gama(z, 10**14)
figure()
boostf_plot(boost_fact, SMboost_fact)

outar=array([ID,z,Lr,stellar_mass,gminusi,boost_fact,SMboost_fact,Lr_raw,stellar_mass_raw]).T
#savetxt('../gcombo/output_pars_gcombo_Lrcomplete.txt',outar)

if movie=='y':
	os.chdir('../movie/')
	movcmd = 'ffmpeg -r 12 -qscale 1 -i %02d.png timelapse.mp4'
	os.system(movcmd)

nw=sum(weight)/max(weight)
nwz=(sum(weight*DL))/sum(weight)
print 'sum(weight)/max(weight)', nw
print 'nwz', nwz
print 'Effective volume = ', sky_vol_weight(2*4,2*12,z,weight)
print 'Effective volume (- local) = ', sky_vol_weight(2*4,2*12,z,weight)-sky_vol(2*4,2*12,min(z))
print 'Full Volume = ', sky_vol(2*4,2*12,max(z))-sky_vol(2*4,2*12,min(z))
endt=time.time()
print "Program took %1.2f seconds" % (endt - startt)
		
	# 
	# #*******************************************************************************
	# # # 								END PROGRAM!								 *
	# # #*****************************************************************************


