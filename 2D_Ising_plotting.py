import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import os
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
import numpy as np


data32 = np.genfromtxt('test_L32.txt')
data16 = np.genfromtxt('test_L16.txt')
data8 = np.genfromtxt('test_L8.txt')
data4 = np.genfromtxt('test_L4.txt')
data = np.genfromtxt('test_L2.txt')

temps16 = np.zeros(len(data))
mags16 = np.zeros(len(data))

temps8 = np.zeros(len(data))
mags8 = np.zeros(len(data))

temps4 = np.zeros(len(data))
mags4 = np.zeros(len(data))

temps32 = np.zeros(len(data))
mags32 = np.zeros(len(data))

temps = np.zeros(len(data))
mags = np.zeros(len(data))


for i in range(len(data)):
    temps16[i]=data16[i][0]
    mags16[i]=data16[i][1]
    temps8[i]=data8[i][0]
    mags8[i]=data8[i][1]
    temps4[i]=data4[i][0]
    mags4[i]=data4[i][1]
    temps32[i]=data32[i][0]
    mags32[i]=data32[i][1]

    temps[i]=data[i][0]
    mags[i]=data[i][1]



plt.figure(figsize=(10,7))
plt.scatter(temps4, mags4, label='N=4', marker='*')
plt.scatter(temps, mags, label='N=2', marker='.')
plt.scatter(temps8, mags8, label='N=8', marker='.')
plt.scatter(temps16, mags16, label='N=16', marker='+')
plt.scatter(temps32, mags32, label='N=32', marker='1')
#plt.xlim([-20,2000])
#plt.ylim([0,1500])
plt.title('Absolute Magnetization/N')
plt.xlabel("Temperature")
plt.ylabel("$\\langle |M| \\rangle$/N")
plt.legend(loc=1, prop={'size': 20})
plt.savefig("Magnetization.png", bbox_inches="tight", pad_inches=0)
plt.show()
