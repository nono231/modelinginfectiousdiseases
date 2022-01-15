#!/usr/bin/env python3

####################################################################
###    This is the PYTHON version of program 2.3 from page 35 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the SIR model with a probability of mortality, and	   #
### unequal births and deaths. This code assumes Density-          #
### Dependent Transmission        								   #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

rho=0.5
nu=mu=1/(70*365.0)
beta=520/365.0
gamma=1/7.0
TS=1.0
ND=1e5
N0=1
X0=0.2
Y0=1e-4
Z0=N0-X0-Y0
INPUT = (X0, Y0, Z0)

def diff_eqs(INP,t):
	'''The main set of equations'''
	Y=np.zeros((3))
	V = INP
	Y[0] = mu - beta * V[0] * V[1] - mu * V[0]
	Y[1] = beta * V[0] * V[1] - (gamma + mu) * V[1]/(1-rho)
	Y[2] = gamma * V[1] - mu * V[2]
	return Y   # For odeint



t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

print (RES)

#Ploting
pl.figure(figsize=(16, 8))
pl.subplot(311)
pl.plot(RES[:,0], '-g', label='Susceptibles')
pl.title('Program_2_3.py')
#pl.xlabel('Time')
pl.ylabel('Susceptibles')
pl.subplot(312)
pl.plot(RES[:,1], '-r', label='Infectious')
#pl.xlabel('Time')
pl.ylabel('Infectious')
pl.subplot(313)
pl.plot(RES[:,2], '-k', label='Recovereds')
pl.plot(sum((RES[:,0],RES[:,1],RES[:,2])), '--k', label='Total Population')
pl.xlabel('Time')
pl.legend(loc=0)
pl.ylabel('Recovereds\nTotal Population')
pl.savefig("Program_2_3.png")
pl.show()
