#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 14:46:08 2020

@author: cloud
"""

import netCDF4 as nc
#import math
import numpy as np
#import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm

file   = "/data/cloud/ucantseeme/hw08/CRH_timeseries.nc"
root   = nc.Dataset(file)
CRH_PDF = root.variables['crh_pdf'][:]
CRH_CDF = np.zeros((444,100))    #from 1979 to 2015, monthly data
ENSO = np.loadtxt("/data/cloud/PDK/ENSO/index.txt")
ENSO_reshape = np.zeros((1,504))

for i in range(42):
    print(i)
    ENSO_reshape[0,12*i:12*(i+1)] = ENSO[i,1:13]

for i in range(444):
    for j in range(100):
        CRH_CDF[i,j] = np.sum(CRH_PDF[i,0:j])


# In[]
t = np.linspace(1,120,120)
ONIp  = np.zeros((120))+1
ONIn = np.zeros((120))-1
fig, ax = plt.subplots(nrows=1, ncols=6, sharey='row', figsize=(12,7))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
for i in range(3):
    ax[2*i].plot(ENSO_reshape[0,120*i:120*(i+1)], t, 'k', linewidth=1.2)
    ax[2*i].plot(ONIp, t, 'r--', linewidth=0.9)
    ax[2*i].plot(ONIn,t, 'b--', linewidth=0.9)
    ax[2*i].grid()
    ax[2*i].set_ylim(1,120)
    ax[2*i].set_xlim(-2.5,2.5)
    ax[2*i].set_xlabel('ONI index')
    ax[2*i].set_yticks(np.linspace(1,121,11))
    ax[2*i].set_yticklabels(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
    ax[2*i].set_title(str(1979+10*i)+' ~ '+str(1979+10*i+9))
    w = ax[2*i+1].pcolormesh(CRH_CDF[120*i:120*(i+1),:], cmap='jet')
    ax[2*i+1].grid()
    ax[2*i+1].set_ylim(1,120)
    ax[2*i+1].set_xlim(0,100)
    ax[2*i+1].set_xlabel('CDF of CRH')
    ax[2*i+1].set_xticks(np.linspace(0,100,3))
    ax[0].set_ylabel('Time (year)')
fig.colorbar(w, ax=ax)
plt.savefig('ONIandCDF_1.png')    


# In[]
p = 0
Elnino_CDF = np.zeros((100))
for i in range(444):
    if ENSO_reshape[0,i]>1:
        p = p + 1
        Elnino_CDF[:] = Elnino_CDF[:] + CRH_CDF[i,:]
Elnino_CDF[:] = Elnino_CDF[:]/p

n = 0
Lanina_CDF = np.zeros((100))
for i in range(444):
    if ENSO_reshape[0,i]<-1:
        n = n + 1
        Lanina_CDF[:] = Lanina_CDF[:] + CRH_CDF[i,:]
Lanina_CDF[:] = Lanina_CDF[:]/n

all_CDF = np.mean(CRH_CDF, axis=0)

# In[]
x = np.linspace(1,100,100)
plt.plot(x,all_CDF, 'k')
plt.plot(x, Elnino_CDF, 'r')
plt.plot(x, Lanina_CDF, 'b')
plt.xlim(0,100)
plt.ylim(0,1)
plt.xlabel('CRH',fontsize=16)
plt.ylabel('CDF',fontsize=16)
plt.legend(['Mean','ONI > 1.0','ONI < -1.0'],fontsize=16)
#plt.show()
plt.savefig('CDFcomparison_1.png')


# In[]
p = 0
Elnino_PDF = np.zeros((100))
for i in range(444):
    if ENSO_reshape[0,i]>1:
        p = p + 1
        Elnino_PDF[:] = Elnino_PDF[:] + CRH_PDF[i,:]
Elnino_PDF[:] = Elnino_PDF[:]/p

n = 0
Lanina_PDF = np.zeros((100))
for i in range(444):
    if ENSO_reshape[0,i]<-1:
        n = n + 1
        Lanina_PDF[:] = Lanina_PDF[:] + CRH_PDF[i,:]
Lanina_PDF[:] = Lanina_PDF[:]/n

all_PDF = np.mean(CRH_PDF, axis=0)

# In[]
x = np.linspace(1,100,100)
#plt.plot(x, all_PDF, 'k', linewidth=0.85)
plt.plot(x, Elnino_PDF, 'r', linewidth=1.1)
plt.plot(x, Lanina_PDF, 'b', linewidth=1.1)
plt.xlim(0,100)
plt.ylim(0,0.02)
plt.xlabel('CRH',fontsize=16)
plt.ylabel('PDF',fontsize=16)
plt.legend(['ONI > 1.0','ONI < -1.0'],fontsize=16)
#plt.show()
plt.savefig('PDFcomparison_nomean.png')

