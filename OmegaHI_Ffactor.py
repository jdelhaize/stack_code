from numpy import *
from pylab import *
from scipy import integrate
import OmegaHI_vars as ohi_vars
#ion()

# Load luminosities
#lums=genfromtxt('NLcc_SGP.txt')

########################################
# def Ffactor(lums, a, phi, alpha, Lstar, typ):
# 
# 	# Plot the histogram to get the count in each bin
# 	figure()
# 	out=hist(lums,bins=200)
# 	
# 	loglums=log10(lums)
# 	#out=hist(loglums,bins=200)
# 	
# 	count=out[0]
# 	binedg=out[1] # The lower edge of the bin 
# 	binw=binedg[1]-binedg[0]
# 	binc=binedg+binw/2
# 	binc=binc[0:-1] # The bin centres
# 	#binc10=10**binc
# 	binc10=binc
# 	
# 	#a = -0.4364 # From my HIPASS fit.
# # 	if survey=='2df':
# # 		a=-0.375 # Karachentsev value is -0.375
# # 	if survey=='g':
# # 		a=-0.60	# Toribio et al 2011
# 		
# 	#a = ohi_vars.alpha_ml
# 		
# 	#(phi, alpha, Lstar) = (1.61e-2, -1.21, 9.54993e9)
# 	#(phi, alpha, Lstar) = (ohi_vars.phi_rhoL, ohi_vars.alpha_rhoL, ohi_vars.Lstar_rhoL)
# 	
# 	if typ=='avMX':
# 		pow1 = a
# 		pow2 = 0.0
# 	if typ=='avMavX':
# 		pow1 = a+1
# 		pow2 = 1.0
# 		
# 	sgpbot=sum(count*binc10**pow1)
# 	Ns=sum(count * binc10**pow2) #=sum(count)
# 
# # 	if typ=='avMX':
# # 		sgpbot=sum(count*binc10**a)
# # 		Ns=len(lums) #=sum(count)
# # 
# # 	if typ=='avMavX':	
# # 		sgpbot=sum(count*binc10**(a+1))
# # 		Ns=sum(count * binc10)
# 	
# 	# According to MZ:
# # 	sgpbot=sum(count*binc10**(a+1))
# # 	Ns=sum(count*binc10)
# 	
# 	def ldens(L):
# 		out=L*(phi/Lstar)*(L/Lstar)**alpha*e**(-L/Lstar)
# 		return out
# 	
# 	rho=integrate.quad(ldens,10**5,10**12)[0]
# 	print 'rho from Ff code =','%.3e' %rho
# 	#rho_gamma = phi*Lstar*()
# 	norb=1.8*10**8
# 	#rho=norb
# 	
# 		
# 	def int2(L):
# 		out=L*(phi/Lstar)*(L/Lstar)**alpha*e**(-L/Lstar)*L**a
# 		return out
# 	 
# 	top=integrate.quad(int2,10**5,10**12)[0] 
# 	#top=16760 if using HIPASS fit. =58851 if using kara
# 	
# 	F = (top/rho)*(Ns/sgpbot)
# 	
# 	print typ, 'factor ML = ', '%.3f' %F
# 
# 	return F
	
########################################
# Using W(L) instead of N(L)
def Ffactor_weight(weight, lums, a, phi, alpha, Lstar, typ):

	# Plot the histogram to get the count in each bin
	figure()
	out=hist(lums,bins=200)
	
	loglums=log10(lums)
	#out=hist(loglums,bins=200)
	
	count=out[0]
	binedg=out[1] # The lower edge of the bin 
	binw=binedg[1]-binedg[0]
	binc=binedg+binw/2
	binc=binc[0:-1] # The bin centres
	#binc10=10**binc
	binc10=binc
	
	wsum=[]
	for i in range(len(binedg)-1):
		msk=where(logical_and(lums>=binedg[i],lums<binedg[i+1]))
		wbin=weight[msk[0]]
		wsum += [ma.sum(wbin)]
		#print shape(msk[0]), ma.sum(wbin)

	count=wsum
	
	#a = -0.4364 # From my HIPASS fit.
# 	if survey=='2df':
# 		a=-0.375 # Karachentsev value is -0.375
# 	if survey=='g':
# 		a=-0.60	# Toribio et al 2011
		
	#a = ohi_vars.alpha_ml
		
	#(phi, alpha, Lstar) = (1.61e-2, -1.21, 9.54993e9)
	#(phi, alpha, Lstar) = (ohi_vars.phi_rhoL, ohi_vars.alpha_rhoL, ohi_vars.Lstar_rhoL)
	
	if typ=='avMX':
		pow1 = a
		pow2 = 0.0
	if typ=='avMavX':
		pow1 = a+1
		pow2 = 1.0
		
	sgpbot=sum(count*binc10**pow1)
	Ns=sum(count * binc10**pow2) #=sum(count)

# 	if typ=='avMX':
# 		sgpbot=sum(count*binc10**a)
# 		Ns=len(lums) #=sum(count)
# 
# 	if typ=='avMavX':	
# 		sgpbot=sum(count*binc10**(a+1))
# 		Ns=sum(count * binc10)
	
	# According to MZ:
# 	sgpbot=sum(count*binc10**(a+1))
# 	Ns=sum(count*binc10)
	
	def ldens(L):
		out=L*(phi/Lstar)*(L/Lstar)**alpha*e**(-L/Lstar)
		return out
	
	rho=integrate.quad(ldens,10**5,10**12)[0]
	print 'rho from Ff code =','%.3e' %rho
	#rho_gamma = phi*Lstar*()
	norb=1.8*10**8
	#rho=norb
	
		
	def int2(L):
		out=L*(phi/Lstar)*(L/Lstar)**alpha*e**(-L/Lstar)*L**a
		return out
	 
	top=integrate.quad(int2,10**5,10**12)[0] 
	#top=16760 if using HIPASS fit. =58851 if using kara
	
	F = (top/rho)*(Ns/sgpbot)
	
	print typ, 'factor ML = ', '%.3f' %F

	return F
	
	
########################################
# Stellar mass version
def Ffactor_sm(sm, a, phi1, phi2, alpha1, alpha2, xstar, typ):

	# Plot the histogram to get the count in each bin
	figure()
	out=hist(sm,bins=200)
	
	logsm=log10(sm)
	#out=hist(logsm,bins=200)
	
	count=out[0]
	binedg=out[1] # The lower edge of the bin 
	binw=binedg[1]-binedg[0]
	binc=binedg+binw/2
	binc=binc[0:-1] # The bin centres
	binc10 = binc
	#binc10=10**binc
	
	if typ=='avMX':
		pow1 = a
		pow2 = 0.0
	if typ=='avMavX':
		pow1 = a+1
		pow2 = 1.0
		
	sgpbot=sum(count*binc10**pow1)
	Ns=sum(count * binc10**pow2) #=sum(count)

# 	if typ=='avMX':
# 		sgpbot=sum(count*binc10**a)
# 		Ns=len(sm) #=sum(count)
# 
# 	if typ=='avMavX':	
# 		sgpbot=sum(count*binc10**(a+1))
# 		Ns=sum(count * binc10)

	schecSM = lambda logx: log(10) * xstar * ( phi1 * (10**logx/xstar)**(alpha1+2) + phi2 * (10**logx/xstar)**(alpha2+2)) * e**(-1*10**logx/xstar)
	XaschecSM = lambda logx: (10**logx)**a * log(10) * xstar * ( phi1 * (10**logx/xstar)**(alpha1+2) + phi2 * (10**logx/xstar)**(alpha2+2)) * e**(-1*10**logx/xstar)
	
	rho=integrate.quad(schecSM,6,12)[0]
	Xarho = integrate.quad(XaschecSM,6,12)[0]
	
	F = (Xarho/rho)*(Ns/sgpbot)
	
	print typ, 'factor SM = ', '%.3f' %F

	return F
	
########################################
# Stellar mass version
def Ffactor_sm_weight(weight, sm, a, phi1, phi2, alpha1, alpha2, xstar, typ):

	# Plot the histogram to get the count in each bin
	figure()
	out=hist(sm,bins=200)
	
	logsm=log10(sm)
	#out=hist(logsm,bins=200)
	
	count=out[0]
	binedg=out[1] # The lower edge of the bin 
	binw=binedg[1]-binedg[0]
	binc=binedg+binw/2
	binc=binc[0:-1] # The bin centres
	binc10 = binc
	#binc10=10**binc
	
	wsum=[]
	for i in range(len(binedg)-1):
		msk=where(logical_and(sm>=binedg[i],sm<binedg[i+1]))
		wbin=weight[msk[0]]
		wsum += [ma.sum(wbin)]

	count=wsum
	
	if typ=='avMX':
		pow1 = a
		pow2 = 0.0
	if typ=='avMavX':
		pow1 = a+1
		pow2 = 1.0
		
	sgpbot=sum(count*binc10**pow1)
	Ns=sum(count * binc10**pow2) #=sum(count)

# 	if typ=='avMX':
# 		sgpbot=sum(count*binc10**a)
# 		Ns=len(sm) #=sum(count)
# 
# 	if typ=='avMavX':	
# 		sgpbot=sum(count*binc10**(a+1))
# 		Ns=sum(count * binc10)

	schecSM = lambda logx: log(10) * xstar * ( phi1 * (10**logx/xstar)**(alpha1+2) + phi2 * (10**logx/xstar)**(alpha2+2)) * e**(-1*10**logx/xstar)
	XaschecSM = lambda logx: (10**logx)**a * log(10) * xstar * ( phi1 * (10**logx/xstar)**(alpha1+2) + phi2 * (10**logx/xstar)**(alpha2+2)) * e**(-1*10**logx/xstar)
	
	rho=integrate.quad(schecSM,6,12)[0]
	Xarho = integrate.quad(XaschecSM,6,12)[0]
	
	F = (Xarho/rho)*(Ns/sgpbot)
	
	print typ, 'factor SM = ', '%.3f' %F

	return F	
	
########################################
# H-factor (for <M>/<L> - MZ's version)
# def Hfactor(lums, a, phi, alpha, Lstar):
# 
# 	# Plot the histogram to get the count in each bin
# 	figure()
# 	out=hist(lums,bins=200)
# 	
# 	loglums=log10(lums)
# 	#out=hist(loglums,bins=200)
# 	
# 	count=out[0]
# 	binedg=out[1] # The lower edge of the bin 
# 	binw=binedg[1]-binedg[0]
# 	binc=binedg+binw/2
# 	binc=binc[0:-1] # The bin centres
# 	#binc10=10**binc
# 	binc10=binc
# 	
# 	#a=-0.4364 #This is from my HIPASS fit. 
# 	#a=-0.375 #Karachentsev value is -0.375
# 	#(phi, alpha, Lstar) = (1.61e-2, -1.21, 9.54993e9)
# 
# 	#a = ohi_vars.alpha_ml
# 
# 	#(phi, alpha, Lstar) = (ohi_vars.phi_rhoL, ohi_vars.alpha_rhoL, ohi_vars.Lstar_rhoL)
# 	
# 	#sgpbot=sum(count*binc10**a)
# 	sgpbot=sum(count*binc10**(a+1))
# 	sgptop=sum(count*binc10)
# 	Ns=len(lums) #=sum(count)
# 	
# 	# According to MZ:
# # 	sgpbot=sum(count*binc10**(a+1))
# # 	Ns=sum(count*binc10)
# 	
# 	
# 	def ldens(L):
# 		out=L*(phi/Lstar)*(L/Lstar)**alpha*e**(-L/Lstar)
# 		return out
# 	
# 	rho=integrate.quad(ldens,10**5,10**12)[0]
# 	#norb=1.8*10**8
# 	
# 		
# 	def int_top(L):
# 		out=L*(phi/Lstar)*(L/Lstar)**alpha*e**(-L/Lstar)*L**(a)
# 		return out
# 		
# 	def int_bot(L):
# 		out=L*(phi/Lstar)*(L/Lstar)**alpha*e**(-L/Lstar)
# 		return out
# 	 
# 	top=integrate.quad(int_top,10**5,10**12)[0] 
# 	bot=integrate.quad(int_bot,10**5,10**12)[0] 
# 	#top=16760 if using HIPASS fit. =58851 if using kara
# 	
# 	#F=(top/rho)*(Ns/sgpbot)
# 	H=(top/rho)*(sgptop/sgpbot)
# 	print 'H factor = ', '%.3f' %H
# 	print 'NOTE: Use this for <M>/<L> case only!!!'
# 
# 	return H
	
########################################
# G-factor (for <M>/<L> - my version)
# def Ffactor(lums):
# 
# 	# Plot the histogram to get the count in each bin
# 	figure()
# 	out=hist(lums,bins=200)
# 	
# 	loglums=log10(lums)
# 	#out=hist(loglums,bins=200)
# 	
# 	count=out[0]
# 	binedg=out[1] # The lower edge of the bin 
# 	binw=binedg[1]-binedg[0]
# 	binc=binedg+binw/2
# 	binc=binc[0:-1] # The bin centres
# 	#binc10=10**binc
# 	binc10=binc
# 	
# 	#a=-0.4364 #This is from my HIPASS fit. 
# 	a=-0.375 #Karachentsev value is -0.375
# 	(phi, alpha, Lstar) = (1.61e-2, -1.21, 9.54993e9)
# 	
# 	#sgpbot=sum(count*binc10**a)
# 	sgpbot=sum(count*binc10**(a+1))
# 	sgptop=sum(count*binc10)
# 	Ns=len(lums) #=sum(count)
# 	
# 	# According to MZ:
# # 	sgpbot=sum(count*binc10**(a+1))
# # 	Ns=sum(count*binc10)
# 	
# 	def ldens(L):
# 		out=L*(phi/Lstar)*(L/Lstar)**alpha*e**(-L/Lstar)
# 		return out
# 	
# 	rho=integrate.quad(ldens,10**5,10**12)[0]
# 	norb=1.8*10**8
# 	#rho=norb
# 	
# 		
# 	def int_top(L):
# 		out=L*(phi/Lstar)*(L/Lstar)**alpha*e**(-L/Lstar)*L**(a+1)
# 		return out
# 		
# 	def int_bot(L):
# 		out=L*(phi/Lstar)*(L/Lstar)**alpha*e**(-L/Lstar)*L
# 		return out
# 	 
# 	top=integrate.quad(int_top,10**5,10**12)[0] 
# 	bot=integrate.quad(int_bot,10**5,10**12)[0] 
# 	#top=16760 if using HIPASS fit. =58851 if using kara
# 	
# 	#F=(top/rho)*(Ns/sgpbot)
# 	F=(top/bot)*(sgptop/sgpbot)
# 	print 'G factor = ', '%.3f' %F
# 
# 	return F

# Run the function
#Ffactor(lums)
