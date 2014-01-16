#!/usr/bin/env python

print " "
print "*** Have you checked the iterators cat_num and logN??? ***"

from OmegaHI_vars import *

# Use parameter file or manual inputs?
ask='y'  #ask=raw_input('Use param file? (y or n): ')
if ask=='y':
	from OmegaHI_params import *
# else:

# 	#***************************
# 	# Ask which processes to run
# 	#***************************


# *************************
# Optical catalogue to load
# *************************
def choose_opcat(cat_num):

	if survey=='2df':
		
		#For SGP
		if pksdat=='p':
			if genshz=='y':
				#opcat='../catalogues/p669+2df_zcat.txt'
				opcat='../shufflez_sgp2/p669+2df_shufflez'+cat_num+'+i.txt'
				#opcat='../sgp_block/sgpya_shz/p669+2df_shufflez'+cat_num+'+i2.txt'
			else:
				opcat='../catalogues/p669+2df_zcat+i.txt'
				#opcat='../sgp_block/sgpya_zcat.txt'
				#opcat='../asym_test/p669+2df_tester3.txt'
				
				#opcat='../gaussnoise_sgp/p669+2df_zcat_gauss1_sel.txt'
				#opcat='../atca_stack/atca+2df_zcat_hiz_match.txt'
				#opcat='../lh/lh_s0_combo.txt' # Laura's stuff only!!!
		
		# For HICAT
		if pksdat=='h' and hicatq=='y':
			if genshz=='y':
				#opcat='../shufflez_hicat/hicat_shufflez'+cat_num+'.txt'
				#opcat='../hicatB/hicatB_shufflez'+cat_num+'.txt'
				opcat='../gaussnoise_hicat/shufflez/hicat_orig_subset_shufflez'+cat_num+'.txt'
			else:
				#opcat='../catalogues/hipass_sourcecat_sel_nonconf+nonex_orig.txt'
				#opcat='../gaussnoise_hicat/hicat_gauss_0030_sel.txt'
				opcat='../gaussnoise_hicat/hicat_gauss_0005.txt'
				#opcat='../catalogues/hipass_sourcecat_sel_nonconf+nonex_subset.txt'
				#opcat='../hicatB/hicat_near_2df.txt'
		
		# For HIPASS		
		if pksdat=='h' and hicatq=='n':
			if genshz=='y':
				#opcat='../shufflez_hipass2/hipass+2df_shufflez'+cat_num+'.txt'
				opcat='../shufflez_hipass3/hipass_shufflez'+cat_num+'.txt'
				#opcat='../hipass/hya_shz/hya_shufflez'+cat_num+'.txt'

			else:
				opcat='../hipass/hipass+2df_vGT300_zcat+i_v2.txt'
				#opcat='../hipass/hya_zcat_sel.txt'
				#opcat='../hicatB/2df_notnear_hicat.txt'
				#opcat='../shufflez_hicat/hicat_shufflez1.txt'
			#opcat='../hipass/non-hipass_sources.txt'
			#opcat='../hipass/bestmatch_h+2.txt'
			#opcat=opcat[:,4::] # Only for bestmatch!!!		
			#opcat='../catalogues/2dfgrs+hipass_zcat.txt'
			#opcat='../shufflez_hipass/hipass+2df_shufflez'+cat_num+'.txt'
		
	# For GAMA
	if survey=='g':
		if genshz=='y':
		
			if csample=='SM':
				ctg='_SM'
			else:
				ctg=''
				
			#opcat='../g'+field+'/shufflez/g'+field+'_tc16+SM+i+sfr_1310_shufflez'+cat_num+'_sel.txt'
			opcat='../gcombo/shufflez_new/gcombo'+ctg+'_shufflez'+cat_num+'.txt'
			
			if field=='hipass':
				#opcat='../ghipass/shufflez/ghipass_shufflez_test'+cat_num+'+SM+i+sfr.txt'
				opcat='../ghipass/shufflez_new/ghipass'+ctg+'_shufflez'+cat_num+'.txt'
		else:
			#opcat='../g9/p669+g9_az+SM+i+sfr_volim1.txt'
			#opcat='../gcombo/p669+gcombo_TC16_volim_z0607.txt'
			#opcat= '../ghipass/ghipass_TC16_volim_z0306.txt'
			
			if csample=='L':
				ctag='_Lrcomplete_all.txt'
			elif csample=='SM':
				ctag='_Lr+SMar_complete.txt'
			else:
				ctag='.txt'
			
			opd={'15':'../g15/p669+g15_TC16+SM+i+sfr_1310.txt', \
'9':'../g9/p669+g9_TC16+SM+i.txt',\
'combo':'../gcombo/p669+gcombo_TC16+SM+i+sfr_1310_nobadMag'+ctag,\
'hipass':'../ghipass/ghipass_TC16+SM+i+sfr'+ctag} #_Lr+SMar_complete.txt #_Lrcomplete_all.txt
			opcat=opd[field]
			
			#opcat='../photz/g9_zLT0.1_spec+phot.txt'
				
	return opcat


# ****************************
# Radio spectra bigcat to load
# ****************************
# Default 'no' for selecting subset of bigcat

# Has been moved to OHI_params.py
# if survey=='2df' and pksdat=='h':
# 	bigcat_subset=raw_input('Use subset of original bigcat? (y or n):')
# 	print 'Have you checked the ID column number??? (Currently set to 42 in _sourceselect)'
# else:
# 	bigcat_subset='n'

def choose_bigcat(cat_num):
	
	if survey=='2df':
		if pksdat=='p':
			if genshz=='y':
				bigcat='../shufflez_sgp/sgp_blk_bigcat_shz'+cat_num+'.txt'
				#bigcat='../sgp_block/sgpya_shz/sgpya_bigcat_shz'+cat_num+'.txt'
			elif noshift=='y':				
				bigcat='../sgp_block/sgp_blockcat_bigcat_noshift.txt'
				#bigcat='../asym_test/test3+n_bigcat_noshift.txt'
				#bigcat='../sgp_block/sgpya2_bigcat_noshift.txt'
			else:
				bigcat='../sgp_block/sgp_blockcat_bigcat.txt'
				#bigcat='../asym_test/test3+n_bigcat.txt'
				#bigcat='../sgp_block/sgpya2_bigcat.txt'
				
				#bigcat='../gaussnoise_sgp/sgp_bigcat_gauss1.txt'
				#bigcat='../flagging/sgp_flag10ng_bigcat.txt'	
				#bigcat='../catalogues/sgp_noshift_bigcat.txt'
				#bigcat='../ra_offset/sgp_bigcat_offset.txt'
				#bigcat='../catalogues/bigcat_sgp.txt'
				#bigcat='../atca_stack/atca_sgp_bigcat_1280.txt'
				#bigcat='../continuum/sgp_p3_combo_bigcat.txt'
			
			#bigcat='../lh/lh_s0_combo_bigcat.txt'  # Laura's stuff only!!!
			
		# For HICAT
		if pksdat=='h' and hicatq=='y':
			if genshz=='y':
				#bigcat='../shufflez_hicat/hicat_bigcat_shz'+cat_num+'.txt'
				bigcat='../gaussnoise_hicat/shufflez/hicat_orig_subset_bigcat_shz'+cat_num+'.txt'
				#bigcat='../hicatB/hicatB_bigcat_shz'+cat_num+'.txt'
			elif noshift=='y':
				#bigcat='../hipass/hicat_bigcat_noshift.txt'
				bigcat='../hicatB/hicatB_bigcat_noshift.txt'
			else:
				#bigcat='../hipass/hicat_bigcat.txt'
				bigcat='../gaussnoise_hicat/hicat_gauss_0005_bigcat.txt'
				#bigcat='../gaussnoise_hicat/hicat_bigcat_1000sel.txt'
				#bigcat='../hicatB/hicatB_bigcat.txt'
		
		# For HIPASS		
		if pksdat=='h' and hicatq=='n':
			if genshz=='y':
				#bigcat='../shufflez_hipass/hipass_shufflez'+cat_num+'_bigcat.txt'
				bigcat='../shufflez_hipass3/hipass_bigcat_shz'+cat_num+'.txt'
				#bigcat='../hipass/hya_shz/hya_bigcat_shz'+cat_num+'.txt'
			elif noshift=='y':
				bigcat='../catalogues/hipass_bigcat_noshift_v2.txt'
			else:
				bigcat='../catalogues/hipass_bigcat_v2.txt'
				#bigcat='../hipass/hya_bigcat.txt'
				#bigcat='../hicatB/2df_hicatB_bigcat.txt'
				
			#bigcat='../hipass/hipass_dets_bigcat.txt' #for bestmatch
		
	if survey=='g':
		if genshz=='y':
		
			if csample=='SM':
				ctg='_SM'
			else:
				ctg=''
		
			#bigcat='../g'+field+'/shufflez/g'+field+'_bigcat_test_shz'+cat_num+'+SM+i.txt'
			#bigcat='../g'+field+'/shufflez/g'+field+'_bigcat_test_shz'+cat_num+'+SM+i+env.txt'
			#bigcat='../g'+field+'/shufflez/g'+field+'_tc16+SM+i+sfr_1310_bigcat_shz'+cat_num+'.txt'
			bigcat='../gcombo/shufflez_new/gcombo_bigcat'+ctg+'_shz'+cat_num+'.txt'
			
			if field=='hipass':
				#bigcat='../ghipass/shufflez/ghipass_bigcat_shz'+cat_num+'+SM+i+sfr.txt'
				bigcat='../ghipass/shufflez_new/ghipass_bigcat'+ctg+'_shz'+cat_num+'.txt'
			
		elif noshift=='y':
			#bigcat='../g'+field+'/g'+field+'_az_bigcat_gD_1290_noshift.txt'
			bigd={'15':'../g15/g15_az_bigcat_gD_noshift.txt', \
'9':'../g9/g9_tc16+SM_bigcat_noshift.txt', \
'combo':'../gcombo/gcombo_az_bigcat_gD_1310_noshift.txt',\
'hipass':'../gama+hipass/hipass+gama_bigcat_noshift.txt'}
			bigcat=bigd[field]
			
		else:
			#bigcat='../g9/g9+oct_az+SM_volim1_bigcat.txt'
			#bigcat='../g9/g9+oct_az_bigcat_tc16sel.txt'
			#bigcat='../gcombo/gcombo_tc16_volim_z0607_bigcat.txt'
			#bigcat= '../ghipass/ghipass_tc16_volim_z0306_bigcat.txt'
			
			if csample=='L':
				ctag='_Lrcomplete_all'
			elif csample=='SM':
				ctag='_Lr+SMar_complete'
			else:
				ctag=''
			
			bigd={'15':'../g15/g15_tc16+SM_1310_bigcat.txt', \
'9':'../g9/g9_tc16+SM_bigcat.txt', \
'combo':'../gcombo/gcombo_tc16_1310_nobadMag'+ctag+'_bigcat.txt',\
'hipass':'../ghipass/ghipass_tc16+SM+i'+ctag+'_bigcat.txt'} #_Lr+SMar_complete # _Lrcomplete_all
			bigcat=bigd[field]
			#bigcat='../photz/g9_zLT0.1_phot_bigcat.txt'
			
		
	return bigcat

# ****************************
# Catalogue of source subset IDs
# ****************************
if opcat_subset=='y':
	subsetcat='../hipass/hicat_index_sigma008.txt'

# ****************************
# Catalogue of confused IDs to load
# ****************************
if confuse=='y':

	def choose_confcat(dv,cat_num):
	
		if survey=='2df' and pksdat=='p':
			#confusedID='../catalogues/p669+2df_confused'+str(dv)+'.txt'
			confusedID='../catalogues/p669+2df_confused1000_bigbeam.txt'
			#confusedID='../catalogues/p669+2df_confused1000.txt'
			#confusedID='../sgp_block/sgpya_shz/p669+2df_confused1000_-1.txt'
		if survey=='2df' and pksdat=='h':
			#confusedID='../hipass/hipass+2df_vGT300_confusedID.txt'
			#confusedID='../hipass/hipass+2df_vGT300_confused'+str(dv)+'.txt'
			confusedID='../hipass/hipass+2df_vGT300_confused2000_bigbeam.txt'
			#confusedID='../hipass/hipass+2df_vGT300_confused2000.txt'
			#confusedID='../hipass/hya_confused2000.txt'

			
		if survey=='g':
		
			if csample=='L':
				ctag='_Lrcomplete_all.txt'
				ctag2='.txt'
			elif csample=='SM':
				ctag='_Lr+SMar_complete.txt'
				ctag2='_SM.txt'
			else:
				ctag='.txt'
				ctag2='.txt'
				
			if field=='hipass':
				confusedID = '../g'+field+'/g'+field+'_confused2000_tc16_bigbeam'+ctag
			elif field=='combo' and genshz=='y':
				#confusedID = '../g'+field+'/shufflez/g'+field+'_confused2000_bigbeam_shz'+cat_num+'.txt'
				confusedID = '../g'+field+'/p669+gcombo_confused2000_shufflez'+ctag2
			else:
				confusedID = '../g'+field+'/p669+g'+field+'_confused2000_tc16_bigbeam_nobadMag'+ctag #_Lr+SMar_complete.txt' #_Lrcomplete_all.txt
			
				
		return confusedID

# ****************************
# Function to choose signal limits if OHI measurement restricted to DV.
# ****************************		
if OHIvDV=='y':

	def choose_siglims(dv):
		siguplim = freq0 + dv/2.
		siglolim = freq0 - dv/2.
		return siguplim,siglolim
	
# **************************************************
# Continuum sources (note delimiter!)
# **************************************************


if survey=='2df':

	if pksdat=='p':
		contincat='../continuum/sgp_nvss_gt200.txt'
		cdelim=' '
	
	if pksdat=='h':
		contincat='../continuum/nvss_hicat_gt200.txt'
		cdelim=' '
	
if survey=='g':
	#contin=genfromtxt('../Omega_HI/g9_continuum_pos.csv',delimiter=',')
	contincat='../continuum/g'+field+'_nvss_gt200.txt'
	cdelim=' '
	

# **************************************************
# Name tag of output saved array
# **************************************************	

# Now outputs the name for both the mass and flux (mean) spectra
if savearr=='y':

	def choose_outspect(cat_num,binnum,logN):
		
		if survey=='2df' and pksdat=='p':
			outdir='../shufflez_sgp3/'+bindic[binby][3]+'bins_3/' #zbins/ #/nvn/'
			#outdir='../shufflez_sgp3/nvn/'
			#outdir='../shufflez_sgp3/nobin/'
			
			#outspect='sgp_shz'+cat_num+'_wsigma_fluxcalib+dv600+boostL3+bigbeam2' # without binning or NvN
			#outspect='sgp_shz'+cat_num+'_N'+str(logN)+'_wsigma' # for NvN
			#outspect='sgp_shz'+cat_num+'_N'+str(logN)+'_wsigma_fluxcalib_bigbeam' # for NvN
			
			#outspect='sgp_shz'+cat_num+'_fluxcalib+dv600+boostL3_bigbeam2_bin'+str(binnum)+'_'+bindic[binby][3] #For binning
			#outspect='sgp_shz'+cat_num+'_fluxcalib+dv600+boostL3_bigbeam2_bin'+str(binnum)+'_'+bindic[binby][3] #For binning
			outspect='sgp_shz'+cat_num+'_bigbeam2_evenL_nocc_bin'+str(binnum)+'_'+bindic[binby][3] #For binning
			
		### For HIPASS  ###
		if survey=='2df' and pksdat=='h' and hicatq=='n':
			#outdir='../shufflez_hipass3/nobin/'
			#outdir='../shufflez_hipass3/nvn/'
			outdir='../shufflez_hipass3/'+bindic[binby][3]+'bins_8/' #zbins/ #/nvn/'
						
			#outspect='hipass_shz'+cat_num+'_wsigma+dv600+boostL+bigbeam' # For normal
			#outspect='hipass_shz'+cat_num+'_N'+str(logN)+'_wsigma_nofit' # For NvN
			
			#outspect='hipass_shz'+cat_num+'_fluxcalib+dv600+boostL3_bigbeam2_bin'+str(binnum)+'_'+bindic[binby][3] #For binning
			outspect='hipass_shz'+cat_num+'_bigbeam2_evenL_nocc_bin'+str(binnum)+'_'+bindic[binby][3] #For binning
		
		# FOR HICAT#	
		if survey=='2df' and pksdat=='h' and hicatq=='y':
			outdir='../gaussnoise_hicat/'
			#outdir='../gaussnoise_hicat/shufflez/'
			
			#outspect='hicatB_shz'+cat_num+'_orig'
			outspect='hicat_gauss_0005'
			#outspect='hicat_orig_subset_shz'+cat_num
			
		if survey=='g':
			if genshz=='y':		
			
				if csample=='SM':
					ctg='_SMcomplete'
				else:
					ctg=''
			
				if confuse=='y':
					cotg='_cc'
				else:
					cotg=''

				if runnvn=='y':					
					outdir='../g'+field+'/shufflez_new/NvN/'					
				if binning=='y':
					outdir='../g'+field+'/shufflez_new/'+bindic[binby][3]+'bins/'
				else:
					outdir='../g'+field+'/shufflez_new/nobin2/'
					
				#outdir='../g'+field+'/shufflez/nobin/'
				#outdir='../g'+field+'/shufflez/nvn/'

				if runnvn=='y':
					outspect='g'+field+ctg+cotg+'_shz'+cat_num+'_N'+str(logN) #for NvN
					#outspect='g'+field+'_shz'+cat_num+'_N'+str(logN) # for NvN

				if binning=='y':
					outspect='g'+field+ctg+cotg+'_shz'+cat_num+'_'+str(binnum)+'_'+bindic[binby][3] #For binning

				else:
					#outspect='g'+field+'_test_shz'+cat_num
					#outspect='g+h_test_v2_shz'+cat_num
					#outspect='g'+field+'_shz'+cat_num+'_tc16_test' #For no binning
					#outspect='g'+field+'_SMcomplete_cc_shz'+cat_num #For no binning
					outspect='g'+field+ctg+cotg+'_shz'+cat_num #For no binning


				
			#else:
				#outdir='../gama+hipass/'
				#outspect='g+h_test_stack'
				#outspect='g'+field+'_test_stack'

		
# Uncomment if want to output both mean and median spectra
# 		outspect_mn=outdir+'mn_'+outspect
# 		outspect_md=outdir+'md_'+outspect

		#return outspect_mn, outspect_md
		
# Uncomment if want to output both flux and mass spectra
		outspect_flux=outdir+'mn_'+outspect+'_flux.txt'
		outspect_mass=outdir+'mn_'+outspect+'_mass.txt'
		outspect_ml=outdir+'mn_'+outspect+'_ml.txt'
		outspect_mmstar=outdir+'mn_'+outspect+'_mmstar.txt'
		
		return outspect_flux, outspect_mass, outspect_ml, outspect_mmstar

# **************************************************
# Files to use for S/N calc
# **************************************************		

def choose_shufdir(binnum,logN,f_or_m,ml='n',spectyp=''):
	if spectyp=='':
		spectyp=f_or_m
	dfm={'m':'mass', 'f':'flux', 'ml':'ml', 'mmstar':'mmstar'}

#####################################
# SGP
	if survey=='2df' and pksdat=='p':
	
		#det_stack='../NvN/md_p669+2df_full_stack.txt' #This can be defined here or as the output of choose_outspect
		
		if binning!='y':
			if ml=='n': # Note that this is different to use_ml!!! (Defined in OmegaHI_parent)
				#shuf_dir=lambda x:'../shufflez_sgp2/mn_sgp_shz'+x+'_corr150_flux.txt' # for non-binning case with confusion correction
				#shuf_dir=lambda x:'../shufflez_sgp3/nobin/mn_sgp_shz'+x+'_bigbeam2_nocalib_flux.txt' # for non-binning case with confusion correction
				shuf_dir=lambda x:'../shufflez_sgp3/nobin/mn_sgp_shz'+x+'_wsigma_fluxcalib+dv600+boostL3+bigbeam2_'+dfm[f_or_m]+'.txt'
				#shuf_dir=lambda x:'../sgp_block/sgpya_shz/mn_sgpya_shz'+x+'bigcontin_'+dfm[f_or_m]+'.txt' # for non-binning case. Use this one!
			if ml=='y':
				shuf_dir=lambda x:'../shufflez_sgp3/nobin/mn_sgp_shz'+x+'_wsigma_fluxcalib+dv600+boostL3+bigbeam2_'+'ml'+'.txt' # for the <ML> error
				#shuf_dir=lambda x:'../sgp_block/sgpya_shz/mn_sgpya_shz'+x+'bigcontin_'+'ml'+'.txt' # for the <ML> error
			
			#shuf_dir=lambda x:'../shufflez_sgp3/nvn/mn_sgp_shz'+x+'_N'+str(logN)+'_wsigma_'+dfm[f_or_m]+'.txt' #for NvN case
		
		if binning=='y':
			#if ml=='n':
				#shuf_dir=lambda x:'../shufflez_sgp2/cbins/mn_sgp_shz'+x+'_bin'+str(binnum)+'_c_mass.txt' #for binning case
				#shuf_dir=lambda x:'../shufflez_sgp3/'+bindic[binby][3]+'bins_3/mn_sgp_shz'+x+'_fluxcalib+dv600+boostL3_bigbeam2_bin'+str(binnum)+'_'+bindic[binby][3]+'_mass.txt' #for binning case. Use this!
			shuf_dir=lambda x:'../shufflez_sgp3/'+bindic[binby][3]+'bins_3/mn_sgp_shz'+x+'_bigbeam2_evenL_nocc_bin'+str(binnum)+'_'+bindic[binby][3]+'_'+dfm[spectyp]+'.txt' #for binning case. Use this!
			#if ml=='y':
			#	shuf_dir=lambda x:'../shufflez_sgp3/'+bindic[binby][3]+'bins_3/mn_sgp_shz'+x+'_fluxcalib+dv600+boostL3_bigbeam2_bin'+str(binnum)+'_'+bindic[binby][3]+'_ml.txt' #This is for <ml> error
		
		
	
#####################################
# HIPASS	
	if survey=='2df' and pksdat=='h' and hicatq=='n':
	
		#det_stack='../shufflez_hipass/md_hipass+2df_full_stack.txt' #This can be defined here or as the output of choose_outspect
		
		
		if binning!='y':
			shuf_dir=lambda x:'../shufflez_hipass3/nobin/mn_hipass_shz'+x+'_wsigma+dv600+boostL+bigbeam2_'+dfm[f_or_m]+'.txt'
			#shuf_dir=lambda x:'../hipass/hya_shz/mn_hya_shz'+x+'_'+dfm[f_or_m]+'.txt'
			if ml=='y':
				shuf_dir=lambda x:'../shufflez_hipass3/nobin/mn_hipass_shz'+x+'_wsigma+dv600+boostL+bigbeam2_'+'ml'+'.txt'
				#shuf_dir=lambda x:'../hipass/hya_shz/mn_hya_shz'+x+'_ml.txt'
			
			#shuf_dir=lambda x:'../shufflez_hipass3/nvn/mn_hipass_shz'+x+'_N'+str(logN)+'_wsigma_'+dfm[f_or_m]+'.txt' #for NvN case
		
		
		if binning=='y':
			#shuf_dir=lambda x:'../shufflez_hipass2/cbins/mn_hipass_shz'+x+'_bin'+str(binnum)+'_c_mass.txt' #for binning case
			#shuf_dir=lambda x:'../shufflez_hipass3/'+bindic[binby][3]+'bins_8/mn_hipass_shz'+x+'_bin'+str(binnum)+'_'+bindic[binby][3]+'_mass.txt'
			#shuf_dir=lambda x:'../shufflez_hipass3/'+bindic[binby][3]+'bins_8/mn_hipass_shz'+x+'_fluxcalib+dv600+boostL3_bigbeam2_bin'+str(binnum)+'_'+bindic[binby][3]+'_mass.txt' #for binning case
			shuf_dir=lambda x:'../shufflez_hipass3/'+bindic[binby][3]+'bins_8/mn_hipass_shz'+x+'_bigbeam2_evenL_nocc_bin'+str(binnum)+'_'+bindic[binby][3]+'_'+dfm[spectyp]+'.txt' #for binning case
			
# 			if ml=='y':
# 				#shuf_dir=lambda x:'../shufflez_hipass3/'+bindic[binby][3]+'bins_8/mn_hipass_shz'+x+'_wsigma+dv600+boostL_bin'+str(binnum)+'_'+bindic[binby][3]+'_ml.txt'
# 				shuf_dir=lambda x:'../shufflez_hipass3/'+bindic[binby][3]+'bins_8/mn_hipass_shz'+x+'_fluxcalib+dv600+boostL3_bigbeam2_bin'+str(binnum)+'_'+bindic[binby][3]+'_ml.txt' #for binning case
				
		
		#shuf_dir=lambda x:'../shufflez_hipass/NvN/'+avtyp+'_hipass_shz'+x+'_N'+str(logN)+'.txt' #for NvN case

#####################################
# HICAT		
	if survey=='2df' and pksdat=='h' and hicatq=='y':
	
		#det_stack='../shufflez_hipass/md_hipass+2df_full_stack.txt' #This can be defined here or as the output of choose_outspect
		
		
		if binning=='n':
			#shuf_dir=lambda x:'../hicatB/mn_hicatB_shz'+x+'_wDL_mass.txt'
			#shuf_dir=lambda x:'../shufflez_hipass2/mn_hipass_shufflez'+x+'_wsig2_mass.txt'
			shuf_dir=lambda x:'../gaussnoise_hicat/shufflez/mn_hicat_orig_subset_shz'+x+'_'+dfm[spectyp]+'.txt'
		
		
		if binning=='y':
			#shuf_dir=lambda x:'../shufflez_hipass2/cbins/mn_hipass_shz'+x+'_bin'+str(binnum)+'_c_mass.txt' #for binning case
		
			shuf_dir=lambda x:'../shufflez_hipass2/'+bindic[binby][3]+'bins_8/mn_hipass_shz'+x+'_bin'+str(binnum)+'_'+bindic[binby][3]+'_mass.txt' #for binning case
		
		#shuf_dir=lambda x:'../shufflez_hipass/NvN/'+avtyp+'_hipass_shz'+x+'_N'+str(logN)+'.txt' #for NvN case

#####################################	
# GAMA
	
	if survey=='g':
	
		if csample=='SM':
			ctg='_SMcomplete'
		else:
			ctg=''
			
		if confuse=='y':
			cotg='_cc'
		else:
			cotg=''
	
		if binning=='n':
# 			if field=='hipass':
# 				#shuf_dir=lambda x:'../ghipass/shufflez/mn_ghipass_test_shz'+x+'_'+dfm[f_or_m]+'.txt' # for non-binning case
# 				#shuf_dir=lambda x:'../ghipass/shufflez_new/nobin/mn_ghipass_SMcomplete_cc_shz'+x+'_'+dfm[spectyp]+'.txt' # for non-binning case
# 				shuf_dir=lambda x:'../ghipass/shufflez_new/nobin/mn_ghipass'+ctg+cotg+'_shz'+x+'_'+dfm[spectyp]+'.txt' # for non-binning case
# 				
# 			else:
# 				#shuf_dir=lambda x:'../g'+field+'/shufflez/nobin/mn_g'+field+'_shz'+x+'_tc16_test_'+dfm[spectyp]+'.txt' # for non-binning case
			shuf_dir=lambda x:'../g'+field+'/shufflez_new/nobin2/mn_g'+field+ctg+cotg+'_shz'+x+'_'+dfm[spectyp]+'.txt' # for non-binning case
		
		if binning=='y':
			#shuf_dir=lambda x:'../g'+field+'/shufflez/'+bindic[binby][3]+'bins/mn_g'+field+'_shz'+x+'_test_'+str(binnum)+'_'+bindic[binby][3]+'_'+dfm[f_or_m]+'.txt' # for non-binning case
#			shuf_dir=lambda x:'../g'+field+'/shufflez/'+bindic[binby][3]+'bins/mn_g'+field+'_shz'+x+'lumraw_evenL_'+str(binnum)+'_'+bindic[binby][3]+'_'+dfm[spectyp]+'.txt' # for binning case
			shuf_dir=lambda x:'../g'+field+'/shufflez_new/'+bindic[binby][3]+'bins/mn_g'+field+ctg+cotg+'_shz'+x+'_'+str(binnum)+'_'+bindic[binby][3]+'_'+dfm[spectyp]+'.txt' # for binning case
			
		
	print shuf_dir
	return shuf_dir
	

# ****************************
# Random array to overlay
# ****************************

if plotoverlay=='y':
	
	def choose_rndary(binnum, spectyp=''):
		if spectyp=='':
			spectyp=f_or_m
		dfm={'m':'mass', 'f':'flux', 'ml':'ml', 'mmstar':'mmstar'}
		
	#####################################	
		if survey=='2df' and pksdat=='p':
			#rnd_ary='../shufflez_sgp2/zbins/mn_sgp_shz1_bin'+str(binnum)+'_z_mass.txt'
			#rnd_ary='../shufflez_sgp2/mn_sgp_shz2_N_final_mass.txt' # Use this for paper plot
			#rnd_ary='../gaussnoise_sgp/mn_gauss3_stack_flux.txt'
		
			if binning=='y':
				rnd_ary='../shufflez_sgp3/zbins_3/mn_sgp_shz1_bin'+str(binnum)+'_'+bindic[binby][3]+'_'+fmdic[f_or_m][0]+'.txt'
			if binning=='n':
				rnd_ary='../shufflez_sgp3/nobin/mn_sgp_shz2_wsigma_fluxcalib+dv600+boostL3+bigbeam2_mass.txt'
				#rnd_ary='../sgp_block/mn_sgpya_full_mass.txt'

	#####################################			
		if survey=='2df' and pksdat=='h' and hicatq=='n':
			if binning=='y':
				rnd_ary='../shufflez_hipass3/Lbins_8/mn_hipass_shz1_fluxcalib+dv600+boostL3_bigbeam2_bin'+str(binnum)+'_'+bindic[binby][3]+'_'+fmdic[f_or_m][0]+'.txt'
			if binning=='n':
				#rnd_ary='../hicatB/mn_hipass_notnear_hicat_mass.txt'
				rnd_ary='../shufflez_hipass3/nobin/mn_hipass_shz2_wsigma_mass.txt'
				#rnd_ary='../hipass/mn_hya_full_mass.txt'

	#####################################			
		if survey=='2df' and pksdat=='h' and hicatq=='y':
			rnd_ary='../hicatB/mn_hicatB_shz1_orig_mass.txt'

	#####################################
		if survey=='g':
			#rnd_ary='../photz/mn_g9_photz_blnk_fit.txt' #photoz
			#rnd_ary='../shufflez_gama/md_g15_shufflez_blnk_fit21.txt'
			#rnd_ary='../g9_flag/mn_g9_blnk_fit_nocont.txt'
			#rnd_ary='../gama+hipass/shufflez/mn_g+h_test_shz1_mass.txt'
			#rnd_ary = '../g9/shufflez/SMbins/mn_g9_shz1_test_cc_0_SM_mass.txt'
			#rnd_ary='../gcombo/shufflez/nobin/mn_gcombo_shz1_tc16_test_'+dfm[spectyp]+'.txt'
			
			rnd_ary='../gcombo/shufflez_new/nobin2/mn_gcombo_SMcomplete_cc_shz1_'+dfm[spectyp]+'.txt' #_SMcomplete_cc_ (for SM) _cc_ (for ML)
			#rnd_ary='../ghipass/shufflez_new/nobin2/mn_ghipass_shz1_'+dfm[spectyp]+'.txt'
			
			if binning=='y':
				rnd_ary='../gcombo/shufflez/Lbins/mn_gcombo_shz1lumraw_evenL_'+str(binnum)+'_'+bindic[binby][3]+'_'+dfm[spectyp]+'.txt'

		
		return rnd_ary

# *******************************************
# Name of saved overlay stacked spectrum plot
# *******************************************
		
def saveplot_name(binnum):
	#name= '../binning2/images/bigbeam_fit2/binfigs_'+surveynm+'_'+bindic[binby][3]+'_'+str(binnum)+'_'+f_or_m+'.pdf'
	name= '../gcombo/images/binfigs_evenL_'+surveynm+'_'+bindic[binby][3]+'_'+str(binnum)+'_'+f_or_m+'.pdf'
	#name= '../binning/images/shz_test1.pdf'
	
	return name


# ****************************
# Name of saved binned info
# ****************************


if binning=='y':

	def binsavenm(typ,binby,survey,tag='_8'): #'_evenL_test3' '_test_cc' '_evenL_bigbeam2_nocc2'
	
	
		outnm='../ghipass/bininfo_new/'+typ+binby+'_'+surveynm+tag+'.txt'
		#outnm='../binning2/bininfo/'+typ+binby+'_'+survey+'_'+surveynm+tag+'.txt'
	
		return outnm

