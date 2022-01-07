#!/usr/bin/env python

####################################################################
###    This is the PYTHON version of program 7.8 from page 285 of  #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the pairwise approximation to the SIS model on a random  #
### network of N individuals, each with n contacts.                #
###																   #
### Initially it is assumed that there is no correlation between   #
### infected and susceptible individuals in the network.           #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

n=4;
tau=0.1;
gamma=0.05;
Y0=1.0;
N0=10000;
ND=MaxTime=100.0;
TS=1.0
X0=N0-Y0
XY0=n*Y0*X0/N0; 

INPUT=np.hstack((X0,XY0))

def diff_eqs(INP,t):  
	'''The main set of equations'''
	Y=np.zeros((2))
	V = INP   
	Y[0] = gamma*(N0-V[0]) - tau*V[1];
	Y[1] = tau*(n-1)*(n*V[0]-V[1])*V[1]/(n*V[0]) + gamma*(n*N0-n*V[0]-V[1]) - tau*V[1]\
	- tau*(n-1)*V[1]*V[1]/(n*V[0]) - gamma*V[1];
	return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)
X = RES[:,0]
XY= RES[:,1]
Y=N0-X


#Ploting
pl.subplot(311)
pl.plot(Y, '-g')
pl.ylabel('Infectious')
pl.subplot(312)
pl.plot(XY, '-r')
pl.ylabel('[XY] pairs')
pl.subplot(313)
pl.plot(XY/((n*X)*(Y/N0)), '-r')
pl.ylabel('relative XY correlation')
pl.xlabel('Time (Years)')
print XY/(n*X*Y/N0)

pl.show()
