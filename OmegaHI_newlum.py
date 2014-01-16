#!/usr/bin/env python

# Need to make this multidimensional for processing gama

from numpy import *
from pylab import *
from utils import *
from constants import *
from OmegaHI_vars import *
from OmegaHI_control import *
import OmegaHI_other as ohiother
import OmegaHI_plots as ohiplots

def account_for_confusion(cat_num, ID, z, ra, dec, lum_raw, L_err, dv, linewidth, boostf, doSM='n', stellar_mass=None, sm_err=None, SMboostf=1):
# Note: boostf is the factor we are using to include the contribution of confused sources that are below the magnitude limit of the optical catalogue. Set to f=1 if no correction needed.
	
	#********************
	# Read in catalogues
	#********************
	confused_cat=open(choose_confcat(dv,cat_num)) # ConfusedID defined in OmegaHI_control
	confused_match=confused_cat.readlines()

	#**********************
	# Loop over each refgal
	#**********************
	gauss_bigary=[]
	counter=0
	countconfused=[]
	countconfused2=[]
	isep=[]
	lum_fact_ary=[]
	zdiffs=[]
	newL=zeros(shape(lum_raw))
	newL_err=zeros(shape(lum_raw))
	
	if doSM=='y':
		newsm=zeros(shape(stellar_mass))
		newsm_err=zeros(shape(stellar_mass))
	
	for cindex in confused_match:
		
		# Skip header lines
		if cindex[0][0]=='#':
			continue
		
		counter+=1
		clist=cindex.split()
		
		# Get the parameters of the refgal
		refgal_id=float(clist[0]) # ID of refgal
		refgal_msk=where(ID==refgal_id)[0]
		
		refgal_ra=ra[refgal_msk]*pi/180.0 # in radians
		refgal_dec=dec[refgal_msk]*pi/180.0
		refgal_vel = z_to_vel(z[refgal_msk])
		refgal_lum=lum_raw[refgal_msk]
		refgal_lum_err=L_err[refgal_msk]
		if doSM=='y':
			refgal_sm=stellar_mass[refgal_msk]
			refgal_sm_err=sm_err[refgal_msk]
		
		# number of galaxies the refgal is confused with
		numgals_conf=size(clist)-1 
		countconfused=countconfused+[numgals_conf]
		
		#print counter, refgal_id, numgals_conf
		
		#****************************************
		# Loop over galaxies confused with refgal
		#****************************************
		if numgals_conf>0:
		
			lum_wgt=[]
			lum_err_wgt=[]										
			gauss_weight_ary=[]
			
			for ncindex in clist[1::]:
				
				nid=float(ncindex)
				idmsk=where(ID==nid)[0]
				
				nra=ra[idmsk]*pi/180.0
				ndec=dec[idmsk]*pi/180.0
				nlum=lum_raw[idmsk]
				nlum_err=L_err[idmsk]
				nlw=linewidth[idmsk]
				nv=z_to_vel(z[idmsk])
				
				zdiffs+=[z[idmsk] - z[refgal_msk]]
				
				if survey=='g' and do_boostf=='y':
					nboostf = boostf[idmsk]
					if doSM=='y':
						nSMboostf = SMboostf[idmsk]
				else:
					nboostf = boostf
					nSMboostf = SMboostf
				
				if usetf=='y':
					
					lum_fact = ohiother.spectral_confusion( nlw, nv, refgal_vel, dv )
					lum_fact_ary+=[lum_fact]
					#print lum_fact
					
					nlum = nlum*lum_fact
					nlum_err = nlum_err * lum_fact
					
					#if counter == 100:
						#print 'refgal_vel, nv, nlw, lum_fact = ', refgal_vel, nv, nlw, lum_fact
				
				#print 'shape(lum_fact)', shape(lum_fact)
				
				# Calculate angular sep with refgal
				nsep=angular_sep(refgal_ra,nra,refgal_dec,ndec)
				isep+=[nsep*60.]
				
				# Find the weighted luminosity (+ error) of the confused galaxy
				gauss_weight=convolve_gaussian(nsep,beamwidth_big)
				
				gauss_weight_ary+=[gauss_weight[0]]
				
				if survey=='g':
					gauss_weight.shape = (len(gauss_weight),1)
				
				nlum_wgt=nlum*gauss_weight
				lum_wgt+=[nlum_wgt]
				
				#print 'nlum_wgt', nlum_wgt
				
				nlum_err_wgt=nlum_err*gauss_weight
				lum_err_wgt+=[nlum_err_wgt]				
				
			#print isep
			lum_wgt=array(lum_wgt)
			lum_err_wgt=array(lum_err_wgt)
			
			#print 'gauss_weight_ary = ', gauss_weight_ary
			gauss_bigary += [sum(gauss_weight_ary)]
			#print 'gauss_bigary', gauss_bigary
			
			#print 'lum_fact_ary=',shape(lum_fact_ary[:,0,:]),lum_fact_ary[:,0,:]
			
			# Count how many confused after spectral axis accounted for
			if usetf=='y':
 				#print 'shape(lum_wgt)',shape(lum_wgt)
# 				print 'lum_wgt', lum_wgt
# 				print 'any(lum_wgt,axis=0 !=0)', any(lum_wgt,axis=0 !=0)
				#print 'any(lum_wgt,axis=2 !=0)',any(lum_wgt,axis=1 !=0)
				#print any(array(lum_wgt)[:,0,:]!=0,axis=1)
				if survey=='g':
					numgals_conf2=len(where(any(array(lum_wgt)[:,0,:]!=0,axis=1)==True)[0])
					#print numgals_conf2
				if survey=='2df':
					numgals_conf2=len(lum_wgt[:,0])
				
				countconfused2+=[numgals_conf2]
			
			finalL=nboostf*ma.sum(lum_wgt,axis=0)+refgal_lum
			finalL_err=(ma.sum(lum_err_wgt**2,axis=0)+refgal_lum_err**2)**0.5
			
			if doSM=='y':
				finalsm=nSMboostf*ma.sum(lum_wgt,axis=0)+refgal_sm
				finalsm_err=(ma.sum(lum_err_wgt**2,axis=0)+refgal_sm_err**2)**0.5
			
			#print 'shape(finalL)',shape(finalL)
		
		else:
			gauss_bigary += [0.0]
			finalL=refgal_lum
			finalL_err=refgal_lum_err
			if doSM=='y':
				finalsm=refgal_sm
				finalsm_err=refgal_sm_err
		
		newL[refgal_msk[0]]=finalL
		newL_err[refgal_msk[0]]=finalL_err
		if doSM=='y':
			newsm[refgal_msk[0]]=finalsm
			newsm_err[refgal_msk[0]]=finalsm_err
		
		#print 'shape(newL)', shape(newL)
			
	gauss_bigary = array(gauss_bigary)
	
	#**********************
	# Plots
	#**********************
	
	#Uncomment if want to see all the plots

	figure()
	hist(gauss_bigary,bins=20,fc='w',ec='k')
	xlabel(r'$\sum b_i$')
	ylabel(r'$\rm{Count}$')
	
	print 'shape(gauss_bigary) = ', shape(gauss_bigary)
	print 'min(gauss_bigary) = ', min(gauss_bigary)
	
	# Plot histogram of angular separations between confused galaxies
# 	figure()
# 	subplot(211)
# 	isep=array(isep)
# 	hist(isep,bins=20)
# 	xlabel('Angular separation (arcmin)')
# 	ylabel('Count')
# 	print 'min angsep=',min(isep)
# 	print 'max angsep=',max(isep)
# 	print'mean angsep=',mean(isep)
# 	
	# Plot histogram of number of confused galaxies
	#subplot(212)
# 	figure()
# 	countconfused=array(countconfused)
# 	hist(countconfused,bins=20)
# 	ylabel('Count')
# 	xlabel('Number of confused sources')
# 	
# 	figure()
# 	plot(z,countconfused,'ko')
# 	xlabel('z')
# 	ylabel('Number of confused sources')
# 	
# 	figure()
# 	zdiffs=array(zdiffs)
# 	hist(zdiffs,bins=20)
# 	xlabel(r'$z_{\rm{conf}}-z_{\rm{central}}$')
# 	print 'Mean z difference = ', mean(zdiffs)
# 	
# 	# Plot histogram of number of confused galaxies (including spectral axis)
# 	if usetf=='y':
# 		#figure()
# 		countconfused2=array(countconfused2)
# 		hist(countconfused2,bins=20)
# 		#ylabel('Count')
# 		#xlabel('Number of confused sources (2)')
# 		
# 	# Plot the luminosity factors
# 	if usetf=='y':
# 		figure()
# 		lum_fact_ary=array(lum_fact_ary)
# 		hist(lum_fact_ary,bins=15)
# 		xlabel(r'$\rm{Luminosity\ factors}$')
# 		ylabel(r'$\rm{Count}$')
# 	
# 	
# 	# Plot original and adjusted luminosity and lum err histograms
# 
# 	ohiplots.lum_compare(newL,lum_raw)
# 	ohiplots.lum_err_compare(newL_err,L_err)
	
	#**********************
	# Count and print stuff
	#**********************
	
	# Count stuff
	if usetf=='y':
		countconfused=countconfused2
		
	#nonc_count=len(where(countconfused==0)[0])
	#nonc_frac=nonc_count/float(len(countconfused))*100
	#conf_count=len(where(countconfused>0)[0])
	#conf_frac=conf_count/float(len(countconfused))*100
	
	# Print stuff
	#print "Num not confused={0}={1}%".format(nonc_count,round(nonc_frac,1))
	#print "Num confused={0}={1}%".format(conf_count,round(conf_frac,1))
	print 'mean num=',ma.mean(countconfused)
	print 'max num=',ma.max(countconfused)
	print 'min num=',ma.min(countconfused)
	
	#savetxt(outname,newL)
	

	out = [newL, newL_err, gauss_bigary]
	if doSM=='y':
		out = [newL, newL_err, gauss_bigary, newsm, newsm_err]
	
	return out

		
