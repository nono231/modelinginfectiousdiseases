#!/usr/bin/env python

####################################################################
###    This is the PYTHON version of program 5.4 from page 186 of  #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is a model for Rabit Hemorrhagic Disease, in which both     #
### transmission rate and birth rates can be seasonally forced.    #
### Note that we are using numbers of individuals in this formulation.#
### Bifurcation plots are not possible with this code.             #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

beta0=0.936;
beta1=0.1;
gamma=0.025;
alpha0=0.02;
alpha1=0.1;
mu=0.01;
m=0.475;
K=10000;

X0=0.5;
Y0=0.01;
N0=0.6;
Years = 60
ND=MaxTime=Years*365;
TS=1.0;

INPUT=np.array((X0, Y0, N0))

def diff_eqs(INP,t):  
	'''The main set of equations'''
	Y=np.zeros((3))
	V = INP   
	t=np.mod(t,365)
	alpha = alpha0*(1+alpha1*np.sin(2*np.pi*t/365))
	beta = beta0*(1+beta1*np.sin(2*np.pi*t/365))
	Y[0] = alpha * V[2]- beta*V[0]*V[1] - (mu+V[2]/K) * V[0] # dX/dt
	Y[1] = beta*V[0]*V[1] - (mu + m + gamma + V[2]/K) * V[1] # dY/dt
	Y[2] = (alpha - mu - V[2]/K) * V[2] - m * V[1]
	return Y   # For odeint


t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

print RES

t=(np.arange(ND)/365.)
#Ploting
pl.subplot(311)
pl.plot(t,RES[1:,0], 'g', label='S')
pl.xlabel ('Time (years)')
pl.ylabel ('Susceptibles')
pl.subplot(312)
pl.plot(t,RES[1:,1], 'r', label='I')
pl.xlabel ('Time (years)')
pl.ylabel ('Infectious')
pl.subplot(313)
pl.plot(t,1-(RES[1:,0]+RES[1:,1]), 'k', label='R')
pl.xlabel ('Time (years)')
pl.ylabel ('Recovereds')

pl.show()
