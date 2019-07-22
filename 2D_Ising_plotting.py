#Author:  Kellie McGuire      kellie@kelliejensen.com
#Plotting and curve fitting for 2D_Ising_model.cpp

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit

#Import data
data128 = np.genfromtxt('test_L128.txt')
data64 = np.genfromtxt('test_L64.txt')
data32 = np.genfromtxt('test_L32.txt')
data16 = np.genfromtxt('test_L16.txt')
data8 = np.genfromtxt('test_L8.txt')
data4 = np.genfromtxt('test_L4.txt')
data3 = np.genfromtxt('test_L3.txt')
data2 = np.genfromtxt('test_L2.txt')

#Initialize data arrays
temps,mags3,mags128,mags64,mags16,mags8,mags4,mags32,mags2,\
E2,E128,E64,E16,E8,E4,E32,C2,C4,C8,C16,C32,C64,C128,Chi8,Chi16,\
Chi2, Chi8, Chi4, Chi32, Chi64, Chi128= [np.zeros(len(data2)) for _ in range(31)]



#Reshape data
for i in range(len(data2)):
    temps[i]=data2[i][0]

    mags16[i]=data16[i][1]
    mags8[i]=data8[i][1]
    mags4[i]=data4[i][1]
    mags32[i]=data32[i][1]
    mags128[i]=data128[i][1]
    mags64[i]=data64[i][1]
    mags2[i]=data2[i][1]
#    mags3[i]=data3[i][1]

    E16[i]=data16[i][2]
    E8[i]=data8[i][2]
    E4[i]=data4[i][2]
    E32[i]=data32[i][2]
    E128[i]=data128[i][2]
    E64[i]=data64[i][2]
    E2[i]=data2[i][2]

    C16[i]=data16[i][3]
    C8[i]=data8[i][3]
    C4[i]=data4[i][3]
    C32[i]=data32[i][3]
    C128[i]=data128[i][3]
    C64[i]=data64[i][3]
    C2[i]=data2[i][3]

    Chi2[i]=data2[i][4]
    Chi16[i]=data16[i][4]
    Chi8[i]=data8[i][4]
    Chi4[i]=data4[i][4]
    Chi64[i]=data64[i][4]
    Chi128[i]=data128[i][4]





#Define function for fitting magnetization curves
def mag_fit(T, A, B, C, a, b, c, d):
#   return (c + (np.exp(-T+t)-np.exp(T+t))/(np.exp(-T+t)+np.exp(T+t)))
   return C+A*np.tanh(a*T+b)+B*np.tanh(c*T+d)
# def mag_fit(T, A, B, C, D):
#     return (C + D*np.tanh(-T*A + B ))

mag_params, mag_param_covars = curve_fit(mag_fit, temps, mags128, p0 = [0.5 ,2.5, 1.25, 1., 1., 1., 1.])
#mag_params, mag_param_covars = curve_fit(mag_fit, temps, mags128, p0 = [0.5 ,2.5, 1.25, 0])

T = np.arange(0.5,5,0.01)





#Plot magnetization
plt.figure(figsize=(10,7))
plt.rcParams.update({'font.size': 13})
plt.scatter(temps, mags2, label='L=2', marker='x', color='r')
#plt.scatter(temps, mags3, label='L=3', marker='4', color='m')
plt.scatter(temps, mags4, label='L=4', marker='+', color='b')
plt.scatter(temps, mags8, label='L=8', marker='|', color='k')
plt.scatter(temps, mags16, label='L=16', marker='1', color ='c')
plt.scatter(temps, mags32, label='L=32', marker='2', color='g')
plt.scatter(temps, mags64, label='L=64', marker='3', color='m')
plt.scatter(temps, mags128, label='L=128', marker='4', color='k')
#plt.plot(T, mag_fit(T, *mag_params))
plt.title('Average Absolute Magnetization per spin ($\\langle$|M|$\\rangle$/N) vs. Temperature (T)')
plt.xlabel("Temperature (T)")
plt.ylabel("$\\langle$|M|$\\rangle$/N")
plt.yticks(np.arange(0, 1.1, 0.1))
plt.xticks(np.arange(0, 5, 0.5))
plt.legend(loc=1, prop={'size': 13})
plt.grid(True)
plt.savefig("Magnetization.png", bbox_inches="tight", pad_inches=0)
plt.show()




#Plot energy
plt.figure(figsize=(10,7))
plt.rcParams.update({'font.size': 13})
plt.scatter(temps, E2, label='L=2', marker='x', color='r')
plt.scatter(temps, E4, label='L=4', marker='+', color='b')
plt.scatter(temps, E8, label='L=8', marker='|', color='k')
plt.scatter(temps, E16, label='L=16', marker='1', color ='c')
plt.scatter(temps, E32, label='L=32', marker='2', color='g')
plt.scatter(temps, E64, label='L=64', marker='3', color='m')
plt.scatter(temps, E128, label='L=128', marker='4', color='k')
plt.title('Energy per spin vs Temp (T)')
plt.xlabel("Temperature (T)")
plt.ylabel("Energy/N")
plt.xticks(np.arange(0, 5, 0.5))
#plt.yticks(np.arange())
plt.legend(loc=4, prop={'size': 13})
plt.grid(True)
plt.savefig("Energy.png", bbox_inches="tight", pad_inches=0)
plt.show()


#Plot magnetic susceptibility
plt.figure(figsize=(10,7))
plt.rcParams.update({'font.size': 13})
# plt.scatter(temps, Chi2, label='L=2', marker='x', color='r')
# plt.scatter(temps, Chi4, label='L=4', marker='+', color='b')
# plt.scatter(temps, Chi8, label='L=8', marker='|', color='k')
# plt.scatter(temps, Chi16, label='L=16', marker='1', color ='c')
plt.scatter(temps, Chi32, label='L=32', marker='2', color='g')
plt.scatter(temps, Chi64, label='L=64', marker='3', color='m')
plt.scatter(temps, Chi128, label='L=128', marker='4', color='k')
plt.title('Magnetic Susceptibility per spin vs Temp (T)')
plt.xticks(np.arange(0, 5, 0.5))
plt.xlabel("Temperature (T)")
plt.ylabel("Susceptibility")
#plt.yticks(np.arange())
plt.legend(loc=4, prop={'size': 13})
plt.grid(True)
plt.savefig("susceptibility.png", bbox_inches="tight", pad_inches=0)
plt.show()

#Plot heat capacity
plt.figure(figsize=(10,7))
plt.rcParams.update({'font.size': 13})
plt.scatter(temps, C2, label='L=2', marker='x', color='r')
plt.scatter(temps, C4, label='L=4', marker='+', color='b')
plt.scatter(temps, C8, label='L=8', marker='|', color='k')
plt.scatter(temps, C16, label='L=16', marker='1', color ='c')
plt.scatter(temps, C32, label='L=32', marker='2', color='g')
plt.scatter(temps, C64, label='L=64', marker='3', color='m')
#plt.scatter(temps, C128, label='L=128', marker='4', color='k')
plt.title('Heat Capacity per spin (C/N) vs Temp (T)')
#plt.ylim(0,4)
plt.xlabel("Temperature (T)")
plt.ylabel("Heat Capacity per Spin (C/N)")
plt.xticks(np.arange(0, 5, 0.5))
#plt.yticks(np.arange())
#plt.xlim(2,3)
plt.legend(loc=4, prop={'size': 13})
plt.grid(True)
plt.savefig("Heat_capacity.png", bbox_inches="tight", pad_inches=0)
plt.show()
