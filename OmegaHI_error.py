#!/usr/bin/env python

from constants import *

def OHI_error(MHI,MHI_err,rhoL,rhoL_err,Lav,Lav_err=0,OHI):

	rhoHI_err=( (rhoL*MHI_err/Lav)**2+(MHI*rhoL_err/Lav)**2+(-1*MHI*rhoL*Lav_err/Lav**2)**2 )**0.5
	
	OHI_err=Omega_HI=(8*pi*G*rhoHI_err)/(3*H0**2)
	
	return OHI_err, rhoHI_err