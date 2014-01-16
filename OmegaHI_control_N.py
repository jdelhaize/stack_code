#!/usr/bin/env python

# *********************
# Which surveys to use
# *********************

survey=raw_input('Use gama or sgp? (g or 2df): ')

if survey=='2df': 
	pksdat=raw_input('--> p669 or hipass? (p or h):')

#***************************
# Ask which processes to run
#***************************

confuse=raw_input('Account for confusion? (y or n): ')	
binning=raw_input('z binning?: ')
if binning=='y':
	nbins=raw_input('--> How many bins?: ')
	nobin=raw_input('--> Value for nobin line? (0 if none): ')
else:
	nbins=1
nocontin=raw_input('Exclude spectra near continuum sources? (y or n): ')
#nocontin='n'
fit=raw_input('Fit baseline to stack? (y or n): ')

blank=raw_input('Blank bad freq channels? (y or n): ')
weighting=raw_input('Use weighting? (y or n): ')
excl=raw_input('Exclude particularly bad spectra? (y or n): ')
noshift=raw_input('Stack without shifting? (y or n): ')
#noshift='n'
savearr=raw_input('Save output array?')

cat_num=3
cat_num=str(cat_num)

# *************************
# Optical catalogue to load
# *************************
def choose_opcat(cat_num):

	if survey=='2df':
		if pksdat=='p':
			#opcat='../catalogues/p669+2df_zcat.txt'
			opcat='../shufflez/p669+2df_shufflez'+cat_num+'.txt'
			#opcat='../atca_stack/atca+2df_zcat_hiz_match.txt'
		
		if pksdat=='h':
			#opcat='../hipass/non-hipass_sources.txt'
			#opcat='../hipass/bestmatch_h+2.txt'
			#opcat=opcat[:,4::] # Only for bestmatch!!!		
			#opcat='../catalogues/2dfgrs+hipass_zcat.txt'
			opcat='../hipass/hipass+2df_vGT300_zcat.txt'
		
	if survey=='g':
		opcat='../catalogues/gama9v6.txt.combined'
		#opcat='../Evenf_trial/g9_evenf_cat.txt'
		
	return opcat


# ****************************
# Radio spectra bigcat to load
# ****************************
# Default 'no' for selecting subset of bigcat
bigcat_subset='n'

if survey=='2df':
	if pksdat=='p':
		#bigcat='../catalogues/bigcat_sgp.txt'
		#bigcat='../catalogues/sgp_noshift_bigcat.txt'
		#bigcat='../ra_offset/sgp_bigcat_offset.txt'
		bigcat='../shufflez/sgp_shufflez'+cat_num+'_bigcat.txt'
		#bigcat='../atca_stack/atca_sgp_bigcat_1280.txt'
		
	if pksdat=='h':
		#bigcat='../hipass/hipass_dets_bigcat.txt' #for bestmatch
		bigcat='../catalogues/hipass_bigcat.txt'
		#bigcat='../catalogues/hipass_bigcat_noshift.txt'
		
		bigcat_subset=raw_input('Use subset of original bigcat? (y or n):')
	
if survey=='g':
	bigcat='../catalogues/bigcat_g9.txt'
	#bigcat='../catalogues/g9_noshift_bigcat.txt'

# ****************************
# Original luminosities to load(?)
# ****************************


# ****************************
# Catalogue of confused IDs to load
# ****************************
if confuse=='y':
	if pksdat=='p':
		confusedID='../catalogues/p669+2df_confusedID2.txt'
	if pksdat=='h':
		confusedID='../hipass/hipass+2df_vGT300_confusedID.txt'

# **************************************************
# Continuum sources (note delimiter!)
# **************************************************	
if survey=='2df':

	if pksdat=='p':
		contincat='../Omega_HI/sgp_continuum_pos.csv'
		cdelim=','
	
	#if pksdat=='h':
		#
	
if survey=='g':
	#contin=genfromtxt('../Omega_HI/g9_continuum_pos.csv',delimiter=',')
	contincat='../Omega_HI/g9_NVSS_contin_200.csv'
	cdelim=','

# **************************************************
# Name tag of output array
# **************************************************	
if savearr=='y':
	outspect='shufflez'+cat_num+'_stack'
