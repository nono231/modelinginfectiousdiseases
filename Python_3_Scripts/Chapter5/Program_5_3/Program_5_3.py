#!/usr/bin/env python3

####################################################################
###    This is the PYTHON version of program 5.3 from page 184 of  #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the simple SIR epidemic with sinusoidal forcing of the   #
### birth rate.				         							   #
### Note: setting beta1 too high can cause numerical difficulties. #
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
beta=17/13.
gamma=1/13.0
alpha0=1/(50*365.0)
alpha1=([0.25])
S0=1/17.
I0=1e-4
ND=MaxTime=60*365
TS=1.0


### This code can also be used to generate bifurcation diagrams, by setting
### beta1 equal to a vector of seasonality rates. The bifurcation diagram is
### constructed using extrapolated initial conditions. Try:
# (beta,gamma,alpha0, alpha1,S0,I0,ND)=(17/13.,1/13., 1./(50*365), np.arange(0.00,1.0,0.01),1/17., 1e-4, 20*365)

INPUT=np.array((S0,I0, 1-S0-I0))

# Model Definition
def diff_eqs(INP,t):
	'''The main set of equations'''
	Y=np.zeros((3))
	V = INP
	t=np.mod(t,365.)
	alpha=alpha0*(1+alpha1*np.sin(2*np.pi*t/365))
	mu=alpha0
	Y[0] = alpha - beta*V[0]*V[1] - mu*V[0]
	Y[1] = beta*V[0]*V[1] - mu*V[1] - gamma*V[1]
	Y[2] = gamma * V[1] - mu * V[2]
	return Y   # For odeint

# Model Run
if len(alpha1)==1:
	alpha1=alpha1[0]
	t_start = 0.0; t_end = ND; t_inc = TS
	t_range = np.arange(t_start, t_end+t_inc, t_inc)
	RES = spi.odeint(diff_eqs,INPUT,t_range)

	print (RES)

	t=(np.arange(ND)/365.)
	#Ploting
	pl.figure(figsize=(16, 8))
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
	pl.savefig('Program_5_3.png')

else:
	if ND < 3650:
		ND = 3650
	alpha2=alpha1
	Bifur_I=np.zeros((len(alpha2),10))
	for i in range(len(alpha2)):
		alpha1 = alpha2[i]

		t_start = 0.0; t_end = ND; t_inc = TS
		t_range = np.arange(t_start, t_end+t_inc, t_inc)
		RES = spi.odeint(diff_eqs,INPUT,t_range)
		INPUT=RES[-1]

		for j in range(10):
			Bifur_I[i,j]=RES[np.arange(ND)[((ND-j*365)-1)],1]

	### Plotting
	pl.figure(figsize=(16, 8))
	pl.plot (alpha2, np.log10(Bifur_I), '.k')
	### if TeX commands do not work comment comment the next line
	pl.xlabel (r'Seasonality, $\alpha_1$')
	pl.ylabel (r'Level of Infection $(log_{10})$')
	pl.savefig('Program_5_3_bifurcation.png')
	### if TeX commands do not work comment uncomment the next line
#	pl.xlabel ('Seasonality, beta1')
#	pl.ylabel ('Level of Infection (log_10)')
pl.show()
