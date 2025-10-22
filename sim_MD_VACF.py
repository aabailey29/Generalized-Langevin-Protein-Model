# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 23:21:49 2025

@author: alana
"""

import numpy as np 
from numpy import *
from scipy.optimize import fsolve
from scipy import integrate
from scipy.stats import moment
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time as tm

#Number of particles
p = 1000

#Number of trials
R = 1000

#Masses
m = np.full(p,1)
m[0] = 10

#Spring Constant
k = 1

#Rest length
L = 1

#Time
steps = 10000         #Number of time steps
time = 1000         #Total time
h = time/steps     #Step size

#kBT calc
kBT = np.zeros(R)

#VACF Calcs
VACF = np.zeros(steps+1, dtype=float)
#store_runs = np.zeros([R,steps])

#Observables
obs = np.zeros(steps)
vt100 = np.zeros(R)



for jj in range (0,R): #loop over number of runs
 
    vtv0 = np.zeros(steps+1)
    
    # Initial Position
    x = np.arange(0,p,dtype=float) + np.random.normal(loc = 0.0, scale = 1.0, size = p)
    
    obs[0] = x[0]
    
    #Initial Velocity
    v = np.random.normal(loc = 0.0, scale = 1.0, size = p)
    vnot = v[0]
    vtv0[0] = vnot*vnot
    
    
    #Initial Force
    F = np.zeros(p)
    F_new = np.zeros(p)
    
    F[0] = k*(x[1]-x[0]-L)
    for n in range (1,p-1):
        F[n] = k*(x[n+1]-x[n]-L) - k*(x[n]-x[n-1]-L)        
    F[p-1] = -k*(x[p-1]-x[p-2]-L) 
    
    
    #Initial Energy
    #KE = 0.5*np.sum(m * np.square(v))
    #print("kinetic energy:", KE)
    
    # VE = 0
    # for n in range (1, p-1):
    #     VE += 0.5 * k * np.square(x[n] - x[n-1] - L)
    
    #MD Solver
    for j in range (0,steps):

        x += h*v + 0.5*np.square(h)*F/m

        
        obs[j] = x[0]
        
        #Updated force  
        F_new[0] = k*(x[1]-x[0]-L)
        for n in range(1,p-1):
            F_new[n] = k*(x[n+1]-x[n]-L) - k*(x[n]-x[n-1]-L)
        F_new[p-1] = -k*(x[p-1]-x[p-2]-L)   

       
        #Velocity
        v += 0.5*h*(F + F_new)/m
        
        #Store new force
        F = F_new
        
        #Kinetic Energy
        #KE = 0.5*np.sum(m * np.square(v))
    
        #Potential Energy
        #VE = 0
        #for n in range (1, p-1):
        #    VE += 0.5*k*np.square(x[n] - x[n-1] - L)
        
        #VACF
        vtv0[j+1] = v[0] * vnot
        #print(vtv0)
        
   
    #Updated Observables
    obs[steps-1] = x[0]
    vt100[jj] = v[0]
    VACF += vtv0
    kBT[jj] = (1/3)*np.sum(m*v*v)/p

#np.mean(kBT)
#plt.hist(vt100, bins=20)
#plt.show()
#ave = np.mean(store_runs, axis = 1)/v0v0
VACF0 = VACF[0]
VACF = VACF/VACF0   
