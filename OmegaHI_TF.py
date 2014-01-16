#!/usr/bin/env python

from scipy import *
from numpy import *
from pylab import *

rcParams['text.usetex']='True'
rcParams['font.sans-serif']= "Helvetica"
rcParams['axes.labelsize']=17
rcParams['font.size']=14
	
from utils import *
from constants import *
from OmegaHI_vars import *
#from OmegaHI_control import *

def tf_W(semimajor, semiminor, tf_offset, tf_slope, absmag_uncorr):

	# Calculate the inclination of the galaxies from their axes ratios, and plot the distribution of i
	inclin=find_incl(semimajor, semiminor)
	inclin=ma.masked_invalid(inclin)
	
	# Find which sources are nans
	msk=where(inclin.mask==False)[0]
	print 'shape(inclin[msk])', shape(inclin[msk])

	# Calculate the rotational velocities using a Tully-Fisher relation, and plot the relation
	rotvel = TF_vel( absmag_uncorr, tf_offset, tf_slope )
	
	if survey=='g':
		print 'Dividing LW by 2!'
		rotvel = rotvel/2.
	
	if survey=='g':
		inclin.shape=(len(inclin),1)
		inclin=inclin*ones(shape(rotvel))
	print 'shape(inclin)', shape(inclin)

	# Calculate the estimated HI linewidth from the rotational velocities and plot the distbn.
	linewidth = rotvel_to_W(rotvel, inclin)

# 	figure()
# 	hist(inclin[msk]*180/pi,bins=20)
# 	ylabel(r'$\rm{Count}$')
# 	xlabel(r'$i\ (\rm{degrees})$')

# 	figure()
# 	plot(rotvel[msk],absmag_uncorr[msk],'ko')
# 	xscale('log')
# 	ylim(-16,-23)
# 	grid()
# 	xlabel(r'$\rm{Rotational\ Velocity\ (km/s)}$')
# 	ylabel(r'$\rm{Absolute\ Magnitude\ (mag)}$')
	
# 	figure()
# 	hist(linewidth[msk],bins=20)
# 	xlabel(r'$\rm{Line\ width\ (km/s)}$')
# 	ylabel(r'$\rm{Count}$')
	
	print 'min(linewidth[msk])', ma.min(linewidth[msk],axis=0), 'max(linewidth[msk])' , ma.max(linewidth[msk],axis=0)
	
	return linewidth