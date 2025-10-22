# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 00:06:24 2025

@author: alana
"""
import numpy as np 
from numpy import *
from scipy.optimize import fsolve
from scipy import integrate
from scipy.stats import moment
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#Laplace Transform
tvec = np.linspace(0, time, steps+1)
svec = np.linspace(0, steps*np.pi/time, steps+1)

N = len(svec)
Fvec = np.zeros(N)
for ii in range (0, N):
    s = svec[ii]
    f = VACF*np.exp(-s*tvec)
    Fvec[ii] = np.trapz(f, tvec)

Kvec = (1/Fvec) - svec
ds = (steps*np.pi/time)/steps    

#Fitting to large s

#Final Conditions
M = 3000
K_0 = Kvec[M]
K_1 = (Kvec[M] - Kvec[M-1])/ds
K_2 = (Kvec[M] - 2*Kvec[M-1] + Kvec[M-2])/(ds*ds)
K_3 = (Kvec[M] - 3*Kvec[M-1] + 3*Kvec[M-2] - Kvec[M-3])/(ds*ds*ds)
K_4 = (Kvec[M] - 4*Kvec[M-1] + 6*Kvec[M-2] - 4*Kvec[M-3] + Kvec[M-4])/(ds*ds*ds*ds)


Z_0 = Fvec[M]
Z_1 = (Fvec[M] - Fvec[M-1])/ds
Z_2 = (Fvec[M] - 2*Fvec[M-1] + Fvec[M-2])/(ds*ds)
Z_3 = (Fvec[M] - 3*Fvec[M-1] + 3*Fvec[M-2] - Fvec[M-3])/(ds*ds*ds)
Z_4 = (Fvec[M] - 4*Fvec[M-1] + 6*Fvec[M-2] - 4*Fvec[M-3] + Fvec[M-4])/(ds*ds*ds*ds)


#finalConditions = [K_0, K_1, K_2, K_3, K_4]
finalConditions = [Z_0, Z_1, Z_2, Z_3, Z_4]

#K_0 Kernel Hierarchy (one parameter)
#fitting to K
Kfit0 = np.full(N, K_0)

#fitting to Z
#Kfit0 = np.full(N, (-Z_2/2))

#K_1 Kernel Hierarchy (two parameters)
#fitting to K
a = K_1
gamma = -K_2/(2*K_1)

#fitting to Z
#a = -Z_3/6
#gamma = -Z_4/(4*Z_3)
# Kfit1 = a/(svec + gamma) 

#K_2 Kernel Hierarchy (three parameters)
# Anot = K_1
# A1 = -K_3/(6*K_1)
# gamma2 = -K_4/(4*K_3)
# Kfit2 = Anot/(svec + (A1/(svec + gamma2)))