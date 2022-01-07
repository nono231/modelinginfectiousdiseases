#!/usr/bin/env python

####################################################################
###    This is the PYTHON version of program 7.1 from page 241 of  #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the SIR epidemic in a metapopulation with "animal-like"  #
### movements of infected or susceptible individuals across the    #
### network.                                                       #
###																   #
### The default model employs global coupling between ALL          #
### subpopulations                                                 #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

n=5;
beta=1.0*np.ones(n);
gamma=0.1*np.ones(n);
nu=0.0001*np.ones(n);
mu=0.0001*np.ones(n);
X0=0.1*np.ones(n);
Y0=0.0*np.ones(n); Y0[0]=0.0001;
m=0.001*np.ones((n,n)); m=m-np.diag(np.diag(m));
ND=MaxTime=2910.0;
TS=1.0

INPUT=np.hstack((X0,Y0))

def diff_eqs(INP,t):  
	'''The main set of equations'''
	Y=np.zeros((2*n))
	V = INP   
	for i in range(n):
		Y[i] = nu[i] - beta[i]*V[i]*V[n+i] - mu[i]*V[i]; 
		Y[n+i] = beta[i]*V[i]*V[n+i] - mu[i]*V[n+i] - gamma[i]*V[n+i]
		for j in range(n):
			Y[i]+=m[i][j]*V[j] - m[j][i]*V[i];
			Y[n+i]+=m[i][j]*V[n+j] - m[j][i]*V[n+i];
	return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

print RES

#Ploting
pl.subplot(211)
for i in range(n):
	pl.plot(t_range/365.0, RES[:,i], color=(0.0,0.3,i/5))
pl.xlabel('Time (Years)')
pl.ylabel('Susceptibles')
pl.subplot(212)
for i in range(n):
	pl.plot(t_range/365.0, RES[:,i+n], color=(0.7,0.0,i/5))
pl.ylabel('Infected')
pl.xlabel('Time (Years)')

pl.show()