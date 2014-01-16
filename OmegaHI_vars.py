#!/usr/bin/env python

# This module defines values for use with OmegaHI_parent.py

from numpy import *
from constants import *
#from OmegaHI_control import *
from OmegaHI_params import *
import utils

#*************************************
# Solar magnitudes in different bands
#*************************************

SolarMag_b=5.29 #Vega. 5.21 AB (according to http://mips.as.arizona.edu/~cnaw/sun.html)
 #5.27 # This is a VEGA Absmag according to http://www.ucolick.org/~cnaw/sun.html. AB is 5.18
SolarMag_ugriz=array([6.38, 5.15, 4.71, 4.56, 4.54]) # These are AB solar magnitudes as in Table 3 Driver (2012). he got them from Hill (2010)  # [6.39,5.07,4.62,4.52,4.48])

# More updated version? http://mips.as.arizona.edu/~cnaw/sun.html


#****************************************
# Luminosity density from various sources
#****************************************
if survey=='2df':
	rhoL=1.82*10**8*h #Bj luminosity density from Norberg 2002
	rhoL_err=0.17*10**8*h # Previously I was using 0.12 for some reason! Why??
if survey=='g':
	#rhoL= array([2.055, 1.97, 2.28, 2.545, 3.035])*10**8  # currently using the average of methods 1 and 2 in Driver (2012) Table 4.
	#rhoL_err = array([0.03, 0.03, 0.03, 0.04, 0.05])*10**8 # Average of upper and lower, then error propogation through the mean between the two methods

	rhoL= array([2.08, 1.98, 2.29, 2.58, 3.07])*10**8  # Summed (Method 2) version in Driver (2012) Table 4.
	rhoL_err = array([0.055,0.055,0.06,0.07,0.085])*10**8 # Average of upper and lower, then error propogation through the mean between the two methods

# Stellar mass density	
	rhoMstar = 3.2E8 *h # From Baldry 2012, transformed from h=0.7 to h=1
	rhoMstar_err = 0.5E8*h
	


#****************
# Redshift and apparent magnitude error
#****************

if survey=='2df': 
	zerr=85.0/c # From Colless
	merr=0.164 # From Norberg 2002. Previously was using 0.15 (wrong section in Norberg)
if survey=='g':
	print '**Using wrong zerr and merr**'
	zerr=65.0/c
	merr=0.05
	
#******************************
# Radio survey parameters
#******************************

beamwidth=(15.5/60.) #in degrees

if survey=='2df' and pksdat=='p':
	#beamwidth_big=(15.5/60.) # Use this for small-beam confusion correction
	beamwidth = (15.5**2 + 14.4**2)**0.5 / 60. # Use this for big-beam continuum exclusion
	beamwidth_big = (15.5**2 + 14.4**2)**0.5 / 60. # for big-beam confusion correction
	#beamwidth=(15.5/60.) # Use this if want small-beam for continuum exclusion
	
if survey=='2df' and pksdat=='h':
	beamwidth_big = 2**0.5 * 15.5 / 60.
	#beamwidth_big=beamwidth #for yaxis=average (small beam)
	
if survey=='g':
	if field=='hipass':
		beamwidth = 2**0.5 * 15.5 / 60.
		beamwidth_big = 2**0.5 * 15.5 / 60.
	else:
		beamwidth = (15.5**2 + 14.4**2)**0.5 / 60. # Use this for big-beam continuum exclusion
		beamwidth_big = (15.5**2 + 14.4**2)**0.5 / 60. # for big-beam confusion correction
	

spec_res=0.0625 #Spectral channel resolution in MHz

#******************************
# Optical survey parameters
#******************************
if survey=='g':

	# Apparent r-band magnitude limit
	if dmu=='I':
		appmag_lim = 19.4
	if dmu=='II':
		appmag_lim = 19.8
	
	# Stellar mass 'flux' completeness limit
	SM_flux_lim = 3.85 # Where SM_flux_lim = log10(SM/DL**2)


#******************************
# Field coordinate ranges
#******************************
if survey=='g':	
	dmax9 = 3.
	dmin9 = -1.
	ramax9 = 141.
	ramin9 = 129.

	dmax15 = 2
	dmin15 = -2
	ramax15 = 223.5
	ramin15 = 211.5

	edgew = beamwidth_big / 2.

#******************************
# Tully-Fisher relation to use
#******************************
if survey=='2df':
	tf_offset=-1.409604
	tf_slope=-8.586856
if survey=='g':
	tf_offset = -0.565
	tf_slope = -8.15
	
#******************************
# M/L = beta * L**(alpha_ml) (or SM)
#******************************
if survey=='2df':
	alpha_ml = -0.375 # Karachentsev et al 2008
if survey=='g':
	alpha_ml = -0.4 #-0.4 # Toribio et al 2011
	alpha_mmstar = -0.7 #-0.69 #-0.724 #-0.724#-0.724 #-0.288 # Huang 2012 -0.288 for logSM<9, -0.724 for logSM>9

#******************************
# Luminosity (and stellar mass) density function schechter parameters
#******************************	
if survey=='2df':
	(phi_rhoL, alpha_rhoL, Lstar_rhoL) = (1.61e-2, -1.21, 9.54993e9) # From Norberg et al
if survey=='g':
	(phi_rhoL, alpha_rhoL, Lstar_rhoL) = (1.24e-2, -1.12, utils.absmag_to_lum(-20.86, SolarMag_ugriz[2])) # From Driver 2012 Table 4 (for r-band only!)
	
	# Stellar mass double schechter function parameters
	#(phi1_rhoMstar, phi2_rhoMstar, alpha1_rhoMstar, alpha2_rhoMstar, xstar_rhoMstar) = (5.235E7, 2.592E8, -0.4701, 0.6577, 10**10.3518)
	(phi1_rhoMstar, phi2_rhoMstar, alpha1_rhoMstar, alpha2_rhoMstar, xstar_rhoMstar) = (1.15E-2, 2.33E-3, -0.34, -1.47, 10**10.35)
	#(1.068E8, 5.291E8, -0.4701, 0.6577, 10**10.66)

# **************************************
# Define integration ranges for HI stack
# **************************************

if movie=='y': #or survey=='g': # to show the continuum pattern
	xlolim=1400  #1400. # For fitting
	xuplim=1440. # For fitting
# elif survey=='2df' and pksdat=='h' and hicatq=='y':
# 	xlolim=1400  #1400. # For fitting
# 	xuplim=1440. # For fitting
else:
	xlolim=1413.5 #1410  # For fitting
	xuplim=1427.3 #1430. # For fitting
	
	#xlolim=1400.0 #1410  # For fitting
	#xuplim=1440.0 #1430. # For fitting

# Ranges for integrating:
if survey=='2df':

	if pksdat=='p':

		#hipass
#  		siglolim=1418.83
#  		siguplim=1421.77	

		# real
		siglolim_orig= 1418.23 # used for the paper old=1418.77 
		#siglolim_orig=1418.8415 # tmp. Even spacing with siguplim
		siguplim_orig= 1421.97 # old=1421.45 
		
		# Gauss 4
		#siglolim = 1415.0
		#siguplim = 1428.72
		
		if singspec=='y':
			# TEMP! Only for source 73620
# 			siglolim=1420.24 
# 			siguplim=1420.97

			# For source 76496
			siglolim_orig=1419.26 
			siguplim_orig=1420.82
		
	if pksdat=='h' and hicatq=='n':
		siglolim_orig=1418.83
		#siglolim_orig=1419.7 #tmp!
		siguplim_orig=1421.77

	if pksdat=='h' and hicatq=='y':
		# HICAT only
		#siglolim_orig=1418.51 #900 1413.07 #500 1414.72 #1417.88 250# 1417.97 150#1418.51 90 #1419.08 zerr=0 #narrow = 1419.52 
		#siguplim_orig=1422.28 #900 1432.26 #500 1426.13  #1423.46  250# 1422.76 150 #1422.28 90 #1421.73 zerr=0 #orig=1421.77 #narrow=1421.29
		
		siglolim_orig=1417.97 
		siguplim_orig=1422.76

if survey=='g':

	if field=='15':
		siglolim_orig=1419.07 
		siguplim_orig=1421.76
		
	if field=='9' or field=='combo' or field=='hipass' or field=='15':	
		siglolim_orig= 1418.984 #1419.03 # Use sgp ranges for now
		siguplim_orig= 1421.827 #1421.45 

	#siglolim=1416.4 # With g9 extension
	#siguplim=1421.7 

# include or OHIvDL=='y' if want to use original siglims for HI integration
if OHIvDV!='y':
	siglolim = siglolim_orig
	siguplim = siguplim_orig
	
# The _orig are used for baseline fitting but not necessarily for integrating the mass.
	
# ****************************
# Bad frequency channels
# ****************************	
if survey=='2df':
	if pksdat=='p':
		
		bad1=[1252.94,1254.94]  # Low-freq cube edge
		bad14=[1252.94,1260.0] # Low-freq cube edge 7MHz chop
		bad2=[1260.0,1290.52]  # Broad-band RFI
		#bad2=[1260.00,1278.35]     # sept09_flag4 only
		bad3=[1279.65,1280.15]
		bad4=[1312.26,1312.83]
		bad5=[1315.38,1316.55]
		bad6=[1365.07,1367.07]   # High-freq cube edge
		bad15=[1360.0,1367.07]   # High-freq cube edge 7MHz chop
		bad7=[1302.85,1303.33]   # sept09_flag4 only
		bad8=[1300.29,1300.69]   # sept09_flag4 only
		bad9=[1301.27,1301.68]   # sept09_flag4 only
		bad10=[1302.93,1303.17]  # Lower edge of 1335 cube
		bad11=[1299.74,1300.33]  # sgp_ng only
		bad12=[1316.87,1317.09]  # Upper edge of 1285 cube
		bad13=[1335.2,1335.6]
		bad16=[1349.68,1350.18]

		#bad_i=[1,2,3,4,5,6,7,8,9]
		#bad_i=[1,3,4,5,6,10,12]
		bad_i=[3,4,5,11,16,1,6]#,14,15]
		#if indivsn=='y':
		#	bad_i=[3,4,5,11,16]	
		

# FOR ATCA DATA ONLY!!!
# 		bad1=[1258.46,1283.00]   #Low-freq cube edge
# 		bad2=[1298.39,1301.71]# Broad-band RFI
# 		bad3=[1286.83,1292.25]
# 		bad4=[1306.56,1308.88]
# 
# 		bad_i=[1,2,3,4]
	
	if pksdat=='h':
		# HIPASS - MW
		bad1=[1418.4,1422.4] #galactic
		bad2=[1367.7,1368.3]
		bad3=[1398.9,1399.6]
		bad4=[1399.7,1400.5]
		bad5=[1407.7,1408.3]
		bad6=[1424.84,1425.33]

		#bad_i=[1,2,3,4,5,6]
		bad_i=[1]
	
if survey=='g':
	bad1=[1253.0,1253.5]   #Low-freq cube edge Note that min and max different for SGP and G9
	#bad2=[1260.0,1280.52]  # Broad-band RFI
	bad2=[1252.0,1280.52] # G9 broad-band to edge of cube
	#bad2=[1259.0,1292.27] # Broad-band RFI G15
	bad4=[1312.30,1313.0] #G15 tmp!
	#bad4=[1312.0,1313.0]
	#bad5=[1315.79,1316.33]
	#bad5=[1315.8,1316.4] #G15 tmp!
	bad5=[1315.5,1316.5] # This seems to be the touchy one!
	bad6=[1366.3,1368.0]   #High-freq cube edge
	bad7=[1311.5,1317.11]  # Excluding "band overlap" region. Temp only!
	#bad8=[1252.0,1290.0] # Broad-band RFI G15 to edge of cube
	bad8=[1252.0,1310.13] # Broad-band RFI G15 to edge of cube
	bad9=[1283.8,1284.7]
	bad10=[1344.94,1346.28]

	if field=='9':
		bad_i=[1,2,4,5,6]
	if field=='15':
		#bad_i=[4,5,6,7,8]
		bad_i=[1,4,6,8]
	if field=='combo':
		bad_i=[1,2,4,5,6]
		
	if field=='hipass':
		# HIPASS - MW
		bad1=[1418.4,1422.4] #galactic
		bad2=[1367.7,1368.3]
		bad3=[1398.9,1399.6]
		bad4=[1399.7,1400.5]
		bad5=[1407.7,1408.3]
		bad6=[1424.84,1425.33]

		#bad_i=[1,2,3,4,5,6]
		bad_i=[1]

# ****************************	
# Define survey name-tags
# ****************************
if survey=='2df' and pksdat=='p':
	surveynm='sgp'
if survey=='2df' and pksdat=='h' and hicatq=='n':
	surveynm='hipass'
if survey=='2df' and pksdat=='h' and hicatq=='y':
	surveynm='hicat'
if survey=='g':
	surveynm='g'+field
	
# ****************************
# Max NvN value
# ****************************
# Note that this is 100*truncate(Log10(Nall))+iterat
# iterat included so that last number in range is included
#if runnvn=='y' or plotnvn=='y':
iterat=10
if movie=='y':
	iterat=2
if survey=='2df' and pksdat=='p':
	maxN=350+iterat #(Down from 360 due to added source-selection criteria)
if survey=='2df' and pksdat=='h':
	maxN=410+iterat
if survey=='g' and field=='9':
	maxN=350+iterat
if surveynm=='g15':
	maxN=330+iterat
if survey=='g' and field=='hipass':
	maxN=290+iterat
if survey=='g' and field=='combo':
	maxN=360+iterat
	
if confuse=='y' and cscale=='y':
	vdict={122:dict(N=0.67, frac=49.038, f_min=1420.12, f_max=1420.70,Col='b'),\
152:dict(N=0.84, frac=57.934, f_min=1420.05, f_max=1420.77,Col='c'),\
187:dict(N=1.04, frac=72.826, f_min=1419.97, f_max=1420.85,Col='m'),\
231:dict(N=1.28, frac=78.710, f_min=1419.86, f_max=1420.96,Col='g'),\
297:dict(N=1.64, frac=90.805, f_min=1419.71, f_max=1421.11,Col='y'),\
353:dict(N=1.96, frac=95.193, f_min=1419.57, f_max=1421.25,Col='k'),\
464:dict(N=2.58, frac=98.961, f_min=1419.31, f_max=1421.51,Col='r')}



# ****************************
# Non-binning values
# ****************************

# Define values of MHI/OHI without binning with/without confusion correction

# Nobin parameters = pksdat : [ nobinM, nobinM_err, nobinO, nobinO_err, nobinO_confused, nobinO_err_confused, M/L, M/L_err, M/L_confused, M/L_err_confused, zmean, Lmean, cmean ]

nobindic={\
'p' : [ 5.01E9, 4.574725E8, 2.27E-4, 2.554804E-5, 4.94E-4, 5.566232E-5, 3.480423E-1, 3.161323E-2, 7.588488E-1, 6.892744E-2 ], \
'h': [ 1.4E9, 6.51184E7, 4.59E-4, 3.670035E-5, 4.91E-4, 3.927227E-5, 7.012222E-1, 3.177902E-2, 7.499505E-1, 3.398736E-2 ] , \
}

# Choose non-binning values to plot
if survey=='2df':
	nobinM = nobindic[pksdat][0]
	nobinM_er = nobindic[pksdat][1]
	
	if confuse=='y':
		nobinO = nobindic[pksdat][2] 
		nobinO_er = nobindic[pksdat][3]
	else:
		nobinO = nobindic[pksdat][4] 
		nobinO_er = nobindic[pksdat][5]


# *******************************************************	
# Define various properties for each binning parameter
# *******************************************************

# Binning dictionary = binby : [binstring_simple, binstring, binprop, savename, nbins ]

if survey=='2df':
	bindic={\
'colour' : [ 'colour', r'$b_{J}-r_{F}$', 'colour', 'c', 3 ] ,\
'z' : [ 'z', r'$z$', 'z', 'z', 3 ] ,\
'L' : [ 'log(L)', r'$\log(L/(L_{\odot}\,h^{-2}))$', 'log10(Lum_all)', 'L', 8 ],\
'lum_raw' : [ 'log(L)', r'$\log(L/(L_{\odot}\,h^{-2}))$', 'log10(lum_raw_all)', 'L', 8 ]
}

if survey=='g':
	bindic={\
'colour' : [ 'colour', r'$b_{J}-r_{F}$', 'colour', 'c', 3 ] ,\
'z' : [ 'z', r'$z$', 'z', 'z', 3 ] ,\
'L' : [ 'log(L_{r})', r'$\log(L_{r}/(L_{\odot}\,h^{-2}))$', 'log10(Lum_all[:,2])', 'L', 3 ],\
'lum_raw' : [ 'log(L_{r})', r'$\log(L_{r}/(L_{\odot}\,h^{-2}))$', 'log10(lum_raw_all[:,2])', 'L', 4 ],\
'SM' : [ 'log(M*)', r'$\log(M*)$', 'stellar_mass', 'SM', 3 ],\
'SM_raw' : [ 'log(M*)', r'$\log(M*)$', 'stellar_mass_raw', 'SM', 8 ],\
'gmi' : [ 'g-i', r'$g-i$', 'gminusi', 'gmi', 3 ],\
'umr' : [ 'u-r', r'$u-r$', 'uminusr', 'umr', 3 ],\
'sfr' : [ 'log(SFR)', r'$\log(\rm{SFR})$', 'sfr', 'sfr', 3 ],\
'sdens' : [ 'log(surface density)', r'$\log(\rm{Environmental\ density})$', 'sdens', 'sdens', 3 ]\
}		

# ****************************
# Flux-or-mass dictionary
# ****************************		
# f_or_m dictionary = f_or_m: [fullname, plotlim]		
fmdic={\
'f': [ 'flux', 'ylim(-4,8)' ] ,\
'm': [ 'mass', 'ylim(-2,5)' ]\
}

	
# ****************************
# Magnitude band dictionary
# ****************************
mbd = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4}
	
	
# ****************************	
# Define the M/L 'luminosity-weighting' correction factor. Found by using no confusion correction and inputing the luminosities to OmegaHI_Ff.py
# ****************************

# Luminosity

# Ffd = { 'sgp': 1.347, 'hipass': 0.6182, 'g9':1.194, 'gcombo':0.817, 'ghipass':0.422 }#'sgp': 1.3489
# Hfd = { 'sgp': 1.679, 'hipass': 1.189, 'gcombo':1.311, 'ghipass':1.043 }

# Using spectrum weighting
Ffd = { 'gcombo':1.872, 'ghipass':1.235 }#'sgp': 1.3489
Hfd = { 'gcombo':2.290, 'ghipass':1.409 }

#Stellar mass

# Ffmstard = { 'gcombo': 0.946, 'ghipass':0.565} #0.651} # 0.828 for linear binning 0.432 for log
# Hfmstard = { 'gcombo': 1.442, 'ghipass':1.129} #0.992} # 1.473, 1.526 log

# Using spectrum weighting
Ffmstard = { 'gcombo': 0.925, 'ghipass':0.355}
Hfmstard = { 'gcombo': 1.417, 'ghipass':1.032}
	
# ****************************	
# Define the factor to boost the confused luminosities by because of optical catalogue magnitude limits (I think this only depends on the <z> of the sample)
# ****************************
#boostf = {'sgp':1.46,'hipass':1.05}
boostf = {'sgp':1.62,'hipass':1.05}

# ****************************	
# Get rid of evolution correction in rhoL, (but keep k-correction) based on values in Norberg 02 fig 8)
# ****************************
#boostrhoL_dct = {'sgp':1.13,'hipass':1.04}
boostrhoL_dct = {'sgp':0.85,'hipass':0.97}

