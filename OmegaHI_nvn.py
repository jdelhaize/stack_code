#!/usr/bin/env python

from numpy import *
#from matplotlib.pyplot import *
from pylab import *
import utils

def nvn_plot(maxN,iterat,NvN_mean,NvN_std,mkr,avtyp,fmt,lab,ec='k',linecol='k--',plotline='y'):
	x=arange(100,maxN,iterat)/100.
# 	y=log10(NvN_mean)
# 	error=log10(NvN_std)

	beamwidth=15.5/60.0
	pksrad=beamwidth/2.0 #deg
	fieldsize=42.0 #deg^2
	obstime=75.0 #hours
	n=10.0**(x)
	nbeams=13.0
	beamarea=pksrad**2*pi*nbeams
	coverage=fieldsize/beamarea
	t_hrs=(obstime/coverage)*n
	t=t_hrs/24.0 #days
	
 	time=t #change later!
	y=log10(NvN_mean)
	#error=NvN_std # This is the normal way!
	error=NvN_std / (array(NvN_mean) * log(10)) 
	#errorbar(x,y,yerr=error,ecolor='k',fmt=None)
	errorbar(x,y,yerr=error,ecolor=ec,fmt=fmt,label=lab) # If no colour bar wanted
	#scatter(x,y,c=log10(time),s=30,hold=True,marker=mkr)#,cmap='autumn')
	
	# Fit a straight line
	#m,c=polyfit(x,y,1)
	
	c, m, covar = utils.fit_line(x, y, error, [1,-0.5])
	
	#Exclude the last few points from the fit
	#c, m, covar = utils.fit_line(x[0:15], y[0:15], error[0:15], [1,-0.5])
	
	merr = sqrt( covar[0][0] )
	cerr = sqrt( covar[1][1] ) # I think!!
	
	print 'm=',m,'+/-',merr,' c=',c,'+/-',cerr
	
	if plotline=='y':
		# Plot the line of best fit, or the Gaussian gradient (-0.5)
		plot(x,m*x+c,linecol)#,label=r'$\rm{fit:}\ $$m=%.3f$' %m)
		#plot(x,-0.5*x+c,'r--')#,label=r'$m=-0.5$')
		
	# 	m,c=polyfit(x[0:20],y[0:20],1)
	# 	print m,c
	# 	plot(x,-0.5*x+c,'k--')
	
	print "min x,y,t =",10**min(x),10**min(y),min(time)
	print "max x,y,t =",10**max(x),10**max(y),max(time)

	#colorbar()
	xlabel(r'$\log(N)$',size='large')
	ylabel(r'$\log(\sigma/\rm{mJy})$',size='large')
	#ylabel('Log(Integration Time - Days)')
	grid()
	
	
	ary=[x,y,t]
	ary=array(ary).T
	print shape(ary)
	return ary

#[x,y,t]=plotting(gama10c,'o')
#ylim(-1.2,0.3)

def nvn_plot2(maxN,iterat,NvN_mean,NvN_std,mkr,avtyp,fmt,lab,ec='k',linecol='k--',plotline='y'):
	x=arange(100,maxN,iterat)/100.
# 	y=log10(NvN_mean)
# 	error=log10(NvN_std)

	beamwidth=15.5/60.0
	pksrad=beamwidth/2.0 #deg
	fieldsize=42.0 #deg^2
	obstime=75.0 #hours
	n=10.0**(x)
	nbeams=13.0
	beamarea=pksrad**2*pi*nbeams
	coverage=fieldsize/beamarea
	t_hrs=(obstime/coverage)*n
	t=t_hrs/24.0 #days
	
 	time=t #change later!
	y=log10(NvN_mean)
	#error=NvN_std # This is the normal way!
	error=NvN_std / (array(NvN_mean) * log(10)) 
	#errorbar(x,y,yerr=error,ecolor='k',fmt=None)
	errorbar(x,y,yerr=error,ecolor=ec,fmt=fmt,label=lab) # If no colour bar wanted
	#scatter(x,y,c=log10(time),s=30,hold=True,marker=mkr)#,cmap='autumn')
	
	# Fit a straight line
	#m,c=polyfit(x,y,1)
	
	#c, m, covar = utils.fit_line(x, y, error, [1,-0.5])
	
	#Exclude the last few points from the fit
	c, m, covar = utils.fit_line(x[0:15], y[0:15], error[0:15], [1,-0.5])
	
	merr = sqrt( covar[0][0] )
	cerr = sqrt( covar[1][1] ) # I think!!
	
	print 'm=',m,'+/-',merr,' c=',c,'+/-',cerr
	
	if plotline=='y':
		# Plot the line of best fit, or the Gaussian gradient (-0.5)
		plot(x,m*x+c,linecol)#,label=r'$\rm{fit:}\ $$m=%.3f$' %m)
		#plot(x,-0.5*x+c,'k--')#,label=r'$m=-0.5$')
		
	# 	m,c=polyfit(x[0:20],y[0:20],1)
	# 	print m,c
	# 	plot(x,-0.5*x+c,'k--')
	
	print "min x,y,t =",10**min(x),10**min(y),min(time)
	print "max x,y,t =",10**max(x),10**max(y),max(time)

	#colorbar()
	xlabel(r'$\log(N)$',size='large')
	ylabel(r'$\log(\sigma/\rm{mJy})$',size='large')
	#ylabel('Log(Integration Time - Days)')
	grid()
	
	
	ary=[x,y,t]
	ary=array(ary).T
	print shape(ary)
	return ary

#[x,y,t]=plotting(gama10c,'o')
#ylim(-1.2,0.3)