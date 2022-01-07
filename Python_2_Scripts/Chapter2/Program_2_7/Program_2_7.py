#!/usr/bin/env python

####################################################################
###    This is the PYTHON version of program 2.7 from page 44 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the SICR which includes a carrier class.		           #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

beta=0.2
epsilon=0.1
gamma=0.01
Gamma=0.001
mu=1/(50*365.0)
q=0.4
S0=0.1
I0=1e-4
C0=1e-3
ND=60*365
TS=1.0
INPUT = (S0, I0, C0)

def diff_eqs(INP,t):  
	'''The main set of equations'''
	Y=np.zeros((3))
	V = INP    
	Y[0] = mu - beta * V[0] * (V[1] + epsilon * V[2]) - mu * V[0]
	Y[1] = beta * V[0] * (V[1] + epsilon * V[2]) - gamma * V[1] - mu * V[1]
	Y[2] = q * gamma * V[1] - Gamma * V[2] - mu * V[2]
	return Y   # For odeint



t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

Rec=1. - (RES[:,0]+RES[:,1]+RES[:,2])
print RES

#Ploting
pl.subplot(311)
pl.plot(RES[:,0], '-g', label='Susceptibles')
pl.title('Program_2_7.py')
pl.xlabel('Time')
pl.ylabel('Susceptibles')
pl.subplot(312)
pl.plot(RES[:,1], '-r', label='Infectious')
pl.xlabel('Time')
pl.ylabel('Infected')
pl.subplot(313)
pl.plot(RES[:,1], '-m', label='Carriers')
pl.xlabel('Time')
pl.ylabel('Carriers')
pl.show()
