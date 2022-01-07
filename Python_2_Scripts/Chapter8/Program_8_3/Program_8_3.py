#!/usr/bin/env python

####################################################################
###    This is the PYTHON version of program 8.3 from page 302 of  #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the SIR epidemic (with equal births and deaths) and      #
### pulsed vaccination. Vaccination starts at time tV, after which #
### a proportion p of all susceptible individuals are vaccinated   #
### every T days            									   #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

beta=520/365.0;
gamma=1/7.0;
mu=1/(70*365.0);
S0=0.1;
I0=1e-4;
p=0.1;
T=2*365;
tV=30*365;
ND=MaxTime=100*365;
TS=1.0
R0=1-S0-I0

INPUT = np.hstack((S0,I0,R0))

def diff_eqs(INP,t):  
	'''The main set of equations'''
	Y=np.zeros((3))
	V = INP   
	Y[0]= mu - beta*V[0]*V[1] - mu*V[0]
	Y[1]= beta*V[0]*V[1] - gamma*V[1] - mu*V[1]
	Y[2]= gamma*V[1] - mu*V[2]
	return Y   # For odeint

t_start = 0.0; t_end = tV; t_inc = TS
t_range1 = np.arange(t_start, t_end+t_inc, t_inc)
t_start = tV; t_end = ND+TS; t_inc = TS
t_range2 = np.arange(tV, t_end, t_inc)
TT = np.hstack((t_range1, t_range2))
RES1 = spi.odeint(diff_eqs,INPUT,t_range1)
i=0
INPUT=RES1[-1]
RES2=np.zeros((3))
while t_range2[i]<ND:
	INPUT[2]=INPUT[2]+INPUT[0]*p;
	INPUT[0]=INPUT[0]*(1-p);
	t_range3 = np.arange(t_range2[i], t_range2[i+T], t_inc)
	tc2 = spi.odeint(diff_eqs,INPUT,t_range3)
	INPUT=tc2[-1]
	RES2= np.vstack((RES2, tc2))
	i+=T
print len (TT)
RES2=RES2[1:,]
S = np.hstack((RES1[:,0],RES2[:,0]))
I = np.hstack((RES1[:,1],RES2[:,1]))
R = np.hstack((RES1[:,2],RES2[:,2]))
TT = np.arange(len(S))
print len(S)

pl.subplot(311)
pl.plot(TT/365.0, S, '-g')
ll=pl.ylim()
tVV=np.repeat([tV/365.],len(ll))
pl.plot(tVV, ll, '--k')
pl.ylabel('Susceptibles')
pl.subplot(312)
pl.plot(TT/365.0, I, '-r')
ll=pl.ylim()
tVV=np.repeat([tV/365.],len(ll))
pl.plot(tVV, ll, '--k')
pl.ylabel('Infected')
pl.subplot(313)
pl.plot(TT/365.0, R, '-k')
ll=pl.ylim()
tVV=np.repeat([tV/365.],len(ll))
pl.plot(tVV, ll, '--k')
pl.ylabel('Recovered')
pl.xlabel('Time (Years)')

pl.show()
