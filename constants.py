#!/usr/bin/env python

from numpy import *

#*****************
# Define constants
#*****************
h=1.0
H0=100.*h #km/s/Mpc
G=4.3*10**-9
c=299792.458 #km/s
tH=9.78*10**9/h #yr
DH=(9.26*10.**25.)/h #In metres
Omega_M=0.3
Omega_lambda=0.7
Omega_k=1.-Omega_M-Omega_lambda

freq0=1420.405751786 # Rest freq of HI in MHz
MpctoMet=3.0857*10**22