#!/usr/bin/env python

from OmegaHI_params import *
#from OmegaHI_control import *
#from constants import *
from numpy import *
from pylab import *


#*****************************************
# Plot OmegaHI/MHI variation over the bins
#*****************************************

def plotpoints(x,y,xer,yer,mkrcol='k',lab=None,fmt='o'):

	errorbar(x,y,yerr=yer,xerr=xer,ecolor=mkrcol,color=mkrcol,fmt=fmt,elinewidth=1.0,label=lab, ms=9,mew=1.0)


def binplotrun(x,y,xer,yer,nobin,nobiner,xlab,ylab,boundaries,mkrcol='k',lcol='r',lab=None,fmt='o'):			

	if calcsn!='y':
		yer=None #No y-error bars if S/N not being calculated
	
	# Plot MHI vs binprop
	
	#figure()
	#print 'Shape y =', shape(y)
	#print 'Shape yer =', shape(yer)
	plotpoints(x,y,xer,yer,mkrcol,lab,fmt)
	
	#xlim(min(boundaries)-abs(min(boundaries)/10.),max(boundaries)+abs(max(boundaries))/10.)
	
	# Will need to change x[0] and xer[0] to the real values! Need to make the edge colour different (not error colour)
	#print x
	#print nobin
	#nobin_x = ( max(boundaries) + min(boundaries) ) / 2.
	#errorbar(nobin_x,nobin,yerr=[nobiner],xerr=[[min(boundaries)],[max(boundaries)]],ecolor=mkrcol,color='w',fmt='s',elinewidth=0.8,ms=7,els='--')
	
# 	hlines(nobin,min(boundaries),max(boundaries),lcol,linestyle='dashed')
# 	hlines(nobin-nobiner,min(boundaries),max(boundaries),lcol,linestyle='dotted')
# 	hlines(nobin+nobiner,min(boundaries),max(boundaries),lcol,linestyle='dotted')
	
	# Draw vertical lines at the bin boundaries
	#vlines(boundaries,min(ylim()),max(ylim()),linestyles='dotted',lw='1.3',color='k')
	
	xlabel(xlab)
	ylabel(ylab)


