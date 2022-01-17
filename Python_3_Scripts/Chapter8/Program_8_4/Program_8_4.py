#!/usr/bin/env python3

####################################################################
###    This is the PYTHON version of program 8.4 from page 306 of  #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the SIR model with two different risk-groups, each group #
### has an associated birth and vaccination rate.				   #
###																   #
### Note, gamma, mu, p, S and I are all vectors. beta is a matrix. #
### tV is the time at which vaccination starts.    				   #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@gmail.com	  #
###################################

# Environment Preparation
import scipy.integrate as spi
import numpy as np
import pylab as pl

# Parameters
beta=np.array([[1., 0.01],[0.01 ,0.1]]);
gamma=np.array([0.1, 0.1]);
mu=np.array([0.2, 0.8])*5e-5;
p0=np.array([0.4, 0.1]);
tV=50*365;
S0=np.array([0.1, 0.7]);
I0=np.array([1e-5, 1e-5]);
ND=MaxTime=100*365;
TS=1.0

INPUT = np.hstack((S0,I0))

# Model Definition
def diff_eqs(INP,t):
	'''The main set of equations'''
	Y=np.zeros((4))
	V = INP
	MU=sum(mu)
	for i in range(2):
		Y[i]= mu[i]*(1-p[i]) - (beta[i,0]*V[2]+beta[i,1]*V[3])*V[i] - MU*V[i]
		Y[i+2]= (beta[i,0]*V[2]+beta[i,1]*V[3])*V[i] - gamma[i]*V[i+2] - MU*V[i+2]
	return Y   # For odeint

# Model Run
t_start = 0.0; t_end = tV; t_inc = TS
t_range1 = np.arange(t_start, t_end, t_inc)
t_start = tV; t_end = ND; t_inc = TS
t_range2 = np.arange(tV, t_end, t_inc)
T = np.hstack((t_range1, t_range2))
p=np.array([0,0])
RES1 = spi.odeint(diff_eqs,INPUT,t_range1)
p=p0
RES2 = spi.odeint(diff_eqs,RES1[-1],t_range2)

# Print Results
print (RES2)

S1 = np.hstack((RES1[:,0],RES2[:,0]))
S2 = np.hstack((RES1[:,1],RES2[:,1]))
I1 = np.hstack((RES1[:,2],RES2[:,2]))
I2 = np.hstack((RES1[:,3],RES2[:,3]))

TT = np.arange(len(S1))

# Ploting
pl.figure(figsize=(16, 8))
pl.subplot(211)
pl.plot(TT/365.0, S1, '-g')
pl.plot(TT/365.0, S2, '--g')
pl.legend(('High Risk','Low Risk'), loc=0)
ll=pl.ylim()
tVV=np.repeat([tV/365.],len(ll))
pl.title("Program_8_4.py")
pl.plot(tVV, ll, '--k')
pl.ylabel('Susceptibles')
pl.subplot(212)
pl.plot(TT/365.0, I1, '-r')
pl.plot(TT/365.0, I2, '--r')
ll=pl.ylim()
tVV=np.repeat([tV/365.],len(ll))
pl.plot(tVV, ll, '--k')
pl.ylabel('Infected')
pl.legend(('High Risk','Low Risk'), loc=0)
pl.xlabel('Time (Years)')
pl.savefig("Program_8_4.png")
pl.show()
