import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

timev = np.zeros(501)
data = np.zeros(501)
f = open(r"C:\Users\alana\OneDrive\Documents\Thesis\VACF_10ps.txt", "r")
line = f.readline()
timev[0] = float(line.split()[0])
data[0] = float(line.split()[1])
# print(line)
i = 1
while line and i < 501:
    line = f.readline()
    str = line
    #print(str)
    spl = str.split()
    #print(spl)
    timev[i] = float(spl[0])
    data[i] = float(spl[1])
    i += 1
          
f.close()

# plt.plot(timev, data)

#Laplace Transform
svector = np.linspace(0, 500*np.pi/10000, 501)

N = len(svector)
Fdata = np.zeros(N)
for ii in range (0, N):
    s = svector[ii]
    func = data*np.exp(-s*timev)
    Fdata[ii] = np.trapz(func, timev)

Kdata = (1/Fdata) - svector
ds = (500*np.pi/10000)/500


#Fixing corrupt data

Kdata[0] = 0.0186132
Kdata[1] = 0.0164846


# Kdata_0 = Kdata[0]
# Kdata_1 = (Kdata[1] - Kdata[0])/ds
# #K_1 = (-3*Kvec[0] + 4*Kvec[1] - Kvec[2])/(2*ds)
# Kdata_2 = (Kdata[2] - 2*Kdata[1] + Kdata[0])/(ds*ds)

# #K_0 Kernel Hierarchy
# Kdatafit0 = np.full(N, Kdata_0)


# #K_1 Kernel Hierarchy
# gamma = -Kdata_0/Kdata_1
# a = -Kdata_0*Kdata_0/Kdata_1
# print("a=", a, "gamma=", gamma)
# Kdatafit1 = a/(svector + gamma)

# #K_2 Kernel Hierarchy
# Anot = ((-4*Kdata_0*Kdata_1) - np.sqrt(16*Kdata_0*Kdata_0*Kdata_1*Kdata_1 - 8*Kdata_2*Kdata_0*Kdata_0*Kdata_0))/(2*Kdata_2)
# gamma2 = (Kdata_0*Anot)/(Anot*Kdata_1 + Kdata_0*Kdata_0)
# A1 = (Anot*Anot)/(Anot*Kdata_1 + Kdata_0*Kdata_0)
# print("Anot = ", Anot, "gamma = ", gamma2, "A1 = ", A1)
# Kdatafit2 = Anot/(svector + (A1/(svector + gamma2)))


#Fit to large s
#Final Conditions
M = 140
K_0 = Kdata[M]
K_1 = (Kdata[M] - Kdata[M-1])/ds
K_2 = (Kdata[M] - 2*data[M-1] + Kdata[M-2])/(ds*ds)
K_3 = (Kdata[M] - 3*Kdata[M-1] + 3*Kdata[M-2] - Kdata[M-3])/(ds*ds*ds)
K_4 = (Kdata[M] - 4*Kdata[M-1] + 6*Kdata[M-2] - 4*Kdata[M-3] + Kdata[M-4])/(ds*ds*ds*ds)


# Z_0 = Fvec[M]
# Z_1 = (Fvec[M] - Fvec[M-1])/ds
# Z_2 = (Fvec[M] - 2*Fvec[M-1] + Fvec[M-2])/(ds*ds)
# Z_3 = (Fvec[M] - 3*Fvec[M-1] + 3*Fvec[M-2] - Fvec[M-3])/(ds*ds*ds)
# Z_4 = (Fvec[M] - 4*Fvec[M-1] + 6*Fvec[M-2] - 4*Fvec[M-3] + Fvec[M-4])/(ds*ds*ds*ds)


finalConditions = [K_0, K_1, K_2, K_3, K_4]
#finalConditions = [Z_0, Z_1, Z_2, Z_3, Z_4]
print(finalConditions)

#K_0 Kernel Hierarchy
#fitting to K
Kdatafit0 = np.full(N, K_0)


#K_1 Kernel Hierarchy
#fitting to K
a = K_1
gamma = -K_2/(2*K_1)
Kdatafit1 = a/(svector + gamma)

#K_2 Kernel Hierarchy
Anot = K_1
A1 = -K_3/(6*K_1)
gamma2 = -K_4/(4*K_3)
Kdatafit2 = Anot/(svector + (A1/(svector + gamma2)))


#Plots
# plt.loglog(svector, Kdata, label = "K(s)")
# plt.loglog(svector, Kdatafit0, label = "Kfit0")
# plt.loglog(svector, Kdatafit1, label = "Kfit1")
# plt.loglog(svector, Kdatafit2, label = "Kfit2")

plt.plot(svector, Kdata, label = "K(s)")
plt.plot(svector, Kdatafit0, label = "Kfit0")
plt.plot(svector, Kdatafit1, label = "Kfit1")
plt.plot(svector, Kdatafit2, label = "Kfit2")
plt.xlabel("s - frequency")
plt.ylabel("K(s)")
plt.legend()
plt.title("LLNL Protein Data")
plt.grid()
plt.xlim([0, 0.5])
plt.ylim([-0.005, 0.025])
plt.show()


#Check Initial Conditions of Fit Curves
#Actual K(s)
Kdata_0 = Kdata[0]
Kdata_1 = (Kdata[1] - Kdata[0])/ds
Kdata_2 = (2*Kdata[0] - 5*Kdata[1] + 4*Kdata[2] - Kdata[3])/(ds*ds)
print("Initial Conditions for K(s):", Kdata_0, Kdata_1, Kdata_2)

#Kfit0
Kdatafit0_0 = Kdatafit0[0]
Kdatafit0_1 = (Kdatafit0[1] - Kdatafit0[0])/ds
Kdatafit0_2 = (2*Kdatafit0[0] - 5*Kdatafit0[1] + 4*Kdatafit0[2] - Kdatafit0[3])/(ds*ds)
print("Initial Conditions for Kfit0:", Kdatafit0_0, Kdatafit0_1, Kdatafit0_2)

#Kfit1
Kdatafit1_0 = Kdatafit1[0]
Kdatafit1_1 = (Kdatafit1[1] - Kdatafit1[0])/ds
Kdatafit1_2 = (2*Kdatafit1[0] - 5*Kdatafit1[1] + 4*Kdatafit1[2] - Kdatafit1[3])/(ds*ds)
print("Initial Conditions for Kfit1:", Kdatafit1_0, Kdatafit1_1, Kdatafit1_2)

#Kfit2
Kdatafit2_0 = Kdatafit2[0]
Kdatafit2_1 = (Kdatafit2[1] - Kdatafit2[0])/ds
Kdatafit2_2 = (2*Kdatafit2[0] - 5*Kdatafit2[1] + 4*Kdatafit2[2] - Kdatafit2[3])/(ds*ds)
print("Initial Conditions for Kfit2:", Kdatafit2_0, Kdatafit2_1, Kdatafit2_2)


#Try Curve Fit
popt, pcov = curve_fit(Kdatafit2, svector, Kdata)
