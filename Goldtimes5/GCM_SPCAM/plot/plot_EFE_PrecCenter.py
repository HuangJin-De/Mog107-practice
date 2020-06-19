import csv
import numpy as np
import numpy.matlib
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt
from scipy.stats import t as t_dist
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy
import sys
sys.path.append('/data/cloud/Goldtimes5/code/')
from PYfunctions.statistics.smooth import smooth121 as sm
from PYfunctions.statistics.smooth import running_mean as rm

data_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/'

data_file= data_path+'EFE_cor.pk'
with open(data_file, 'rb') as f:
    EFE= pickle.load(f)

data_file= data_path+'Prec_Center.pk'
with open(data_file, 'rb') as f:
    PC= pickle.load(f)

EFEall= EFE.reshape((EFE.size))
PCall= PC.reshape((PC.size))

EFEs= rm(EFEall,15)
EFEs1= EFEs.reshape((5,365))[1,:]

EFE1= EFE[1,:]
PC1= PC[1,:]

#  5 yr
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,6,1/365) ,np.zeros(EFEall.size),'k',linewidth=0.5)
plt.plot(np.arange(1,6,1/365) ,EFEall,'k',linewidth=.8, label= 'Energy flux equator')
plt.plot(np.arange(1,6,1/365) ,EFEs  ,'r',linewidth=1, label= 'Energy flux equator (15day)')
plt.plot(np.arange(1,6,1/365) ,PCall ,'b',linewidth=1, label= 'Precipitation maxima')
ax=plt.gca()
ax.set_xlim(1, 6)
ax.set_ylim(-20, 30)
plt.legend(loc='upper right', fontsize=8)
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/EFE_PC/All_5yr.png'
plt.savefig(picsave)
plt.show()

# 1 yr
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,5*365+1,1) ,np.zeros(EFEall.size),'k',linewidth=0.5)
plt.plot(np.arange(1,366) ,EFE1  ,'k',linewidth=.8, label= 'Energy flux equator')
plt.plot(np.arange(1,366) ,EFEs1 ,'r',linewidth=1, label= 'Energy flux equator (15day)')
plt.plot(np.arange(1,366) ,PC1   ,'b',linewidth=1, label= 'Precipitation maxima')
ax=plt.gca()
ax.set_xlim(1, 366)
ax.set_ylim(-20, 30)
plt.legend(loc='upper right', fontsize=8)
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/EFE_PC/yr2.png'
plt.savefig(picsave)
plt.show()

