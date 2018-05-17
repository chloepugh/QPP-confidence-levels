#!/usr/bin/env python

"""
Script to help calculate 95% and 99% confidence levels for data with a power law 
shape in the power spectrum (e.g. data with red noise).

Reads in a file named 'pow_fit_err.dat' containing a single column of uncertainties
on the power model for each frequency.

Outputs 2 files, one for 95% level, one for 99% level. Note that these levels are
normalised, and defined in linear space, so need to do pow_fit+ALOG10(0.5*cl95*MEAN(10.^(pow-pow_fit))) and 
pow_fit+ALOG10(0.5*cl99*MEAN(10.^(pow-pow_fit))) in IDL, where pow is an array containing the periodogram powers, and pow_fit is the fitted power law model.

Follows method described by Vaughan, S. 2005, A&A, 431, 391.

Created 23-Nov-2016 by Chloe Pugh
"""

import numpy as np
import math
from scipy.integrate import quad
from scipy.optimize import brentq

#Read in file containing uncertainties on fitted log powers
pow_fit_err = np.loadtxt('pow_fit_err.dat', dtype='f8')
pow_fit_err *= math.log(10.) #Convert to S_j, defined in Vaughan 2005

#Pre-define a few variables to speed up calculations
npow = len(pow_fit_err)
sq2pi = math.sqrt(2.*math.pi)
prob95 = 0.05/float(npow) #For the 95% confidence level
prob99 = 0.01/float(npow) #For the 99% confidence level

#The probability is equal to the integral of this function
def mypdf(w,x,sj):
	logw = np.log(w, dtype='f8')/sj
	return (np.exp((-0.5*logw*logw)-(0.5*x*w))/(sq2pi*sj*w))

#Integral of probability density (subtracting prob so the returned quantity
#should be equal to zero, then can use brentq to solve):
def myprob(x, sj, prob, w1, w2):
	result = quad(mypdf, w1, w2, args=(x, sj))
	#Using logs here in case numbers get too extreme
	return np.log(result[0])-np.log(prob)

#Solve integral for each model power uncertainty
def calc_cl(sj, prob):
	#If the uncertainty is very small then restrict the limits of integration
	w1,w2 = (0., np.inf) if sj > 0.01 else (1.-6.*sj, 1.+6.*sj)
	#Use brentq to find x
	return brentq(myprob, 0., 500., args=(sj,prob,w1,w2))

cl95 = np.array([calc_cl(sj, prob95) for sj in pow_fit_err])
cl99 = np.array([calc_cl(sj, prob99) for sj in pow_fit_err])

#Save confidence levels to files:
np.savetxt('cl95.dat', cl95, fmt="%f0")
np.savetxt('cl99.dat', cl99, fmt="%f0")

print('Python: Done.')
