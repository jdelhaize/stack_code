#!/usr/bin/env python

# This is a parameter file that can be used instead of manual inputs at the start of the OmegaHI_parent program


survey='2df' #'2df'   #'Use gama or 2dF? (g or 2df): ')

if survey=='2df': 
	pksdat='p' #'--> p669 or hipass? (p or h): ')
	hicatq='n' #'--> Using hicat? (y or n): ')
if survey=='g':
	field='9' #' --> Which gama field? (9 or 15 or combo): ')
	dmu='I'   # Which GAMA DMU version to use? (I or II)

opcat_subset='n' # Select a subset of the opcat sources?

#'Use subset of original bigcat? y for normal, n when using shufflez cats #always 'n' now
if survey=='2df' and pksdat=='h':
	bigcat_subset='n'
	indexcol=42
else:
	bigcat_subset='n'

binning='y'             #'z binning?: ')
binby='L'          # 'z', 'L' or 'colour' or 'SM', 'gmi', 'umr', 'sfr'
binchoice='r'		# Bin by even what? 'N'=even numbers of galaxies in each bin, 'r'=even bin ranges
if binning=='y':
	#nbins=2         #'--> How many bins?: ') NOTE: Now defined in bindic in _vars.py
	savebin='y'     #'--> Save binned info?: ')
# else:
# 	nbins=1
	# Note: further binning plotting parameters in OmegaHI_vars.py
	
confuse='n'  #'Account for confusion? (y or n): ')	
usetf='n'    # Use tully-fisher method to decide spectral confusion?
OHIvDV='n'  # Make OHI vs DV plot? (Changes the siglim used for HI mass integration)


do_boostf='y'  # Include factor in confusion correction to account for sources below the optical catalogue magnitude limit 

f_factor='y' # Apply F-factor correction to OmegaHI?

cscale='n'  # Scale up final mass? (Only for calibration against HICAT)

fluxcalib='y' # Scale up fluxes due to Parkes calibration? (y only for SGP)

nocontin='y'  #'Exclude spectra near continuum sources? (y or n): ') (n for hipass, y for sgp)
if nocontin=='y':
	sepfac=1.0   #'----> 1.0 if excluding full beamwidth, 2.0 if half a beamwidth
	
colour_msk='n'  #'Exclude sources with bad colours? (y or n): ') NOTE: Check ranges in OHI_sourceselect.py (y for HIPASS, n for SGP)

completeness_cut='n' # Exclude galaxies with r>19.8 or r>19.4


fit='y'  #'Fit baseline to stack? (y or n): ')

indivfit='n'  # Fit a baseline to the individual mass spectra and use this for final stack values? (y or n)

blank='y'   #'Blank bad freq channels? (y or n): ')

weighting='y'   #'Use weighting? (y or n): ')

rfi_mask='y'    #'Exclude full spectrum if source in RFI zone? (y or n): ')

deccut='n'

absmag_msk='n' # Mask GAMA sources with bad absmag (and stellar masses)

excl='n'      #'Exclude particularly bad spectra? (y or n): ')

fit_sin='n'  #'Fit and remove sinusoid to stacked mass spectrum? (y or n): ')

fit_gaussian='n'

boost_rhoL='y'  # To change the luminosity density from z=0 to median(z)

showplots='n'

avtyp='mn'        #'Use mean or median stats? (mn or md):')

noshift='n'   #'Stack without shifting? (y or n): ')

savearr='n'   #'Save output array? (y or n): ')

magband='r'

use_ml='y'  # The output spectrum to save and use for S/N is the <M/L> stack?

genshz='n'    #'Generate stacked shuffled spectra? (y or n):
z_orig_col=-2 # Column with original (not shuffled) redshifts
q_orig_col=-1

runnvn='n'         #'Run NvN? (y or n): ')
movie='n'			# Make NvN movie?

f_or_m='m'       #'Use frequency or mass? (f or m): ')

calcsn='y'       #'Calculate S/N? (y or n): ')

plotoverlay='n'		#'Plot stack overlay? (y or n): ')
if plotoverlay=='y':
	plotpar='mmstar'  # m, f, ml, mmstar
	saveplot='n'     # Save the overlaid plot

indivsn='n' 		#'Investigate individual spectra? (y or n): ')

singspec='n'    #'Run process for only a single spectrum? (y or n): ')