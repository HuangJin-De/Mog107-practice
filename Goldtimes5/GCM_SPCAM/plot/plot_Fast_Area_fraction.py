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
from PYfunctions.statistics.smooth import running_mean
from PYfunctions.statistics.smooth import panta_mean as panta

#=== path
data_path_map= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/'
filesave  = '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/Area_fraction/'

#=== parameter
Lv= 2.477*10**6 # J/kg
rho= 1000 # kg/m3

#=== map
GCM_temp= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(GCM_temp,'r') as D:
  lat=D.variables['lat'][:] # level in hPa
  lon=D.variables['lon'][:] # level in hPa

lat_mask = (lat> -30) & (lat< 30)
lat      = lat[lat_mask]
A        = np.cos (lat * np.pi /180)
A        = np.tile(A, [1, len(lon)] )
A_r   = A.reshape(A.size) # Area of each grid (cosine latitude)

#=== load data
data_file= data_path_map+'PRECT.pk'
with open(data_file, 'rb') as f:
    tmp= pickle.load(f)

tmp       = tmp.reshape(tmp.shape[0]*tmp.shape[1], tmp.shape[2], tmp.shape[3])
tmp       = panta(tmp)
#tmp_smooth= running_mean(tmp,6)
data      = tmp #- tmp_smooth+ np.nanmean(tmp_smooth, axis=0) # 5 - 30 days
data      = data [:, lat_mask, :] # tropics
data      = data [3:-3, :, :] # remove nans due to smoothing

Nd, Y, X = data.shape
print (data.shape)
A_frac   = np.zeros([3,Nd]) # 3 crit ratio: 70, 80 90 %
P_tot    = np.zeros([3,Nd])
P_check  = np.zeros([3,Nd])
for dd in range(Nd):
  tmp   = data[dd,:,:]
  tmp_r = tmp.reshape((tmp.size))
  #== sort by prec (large -> small)
  idx   = np.argsort(tmp_r)[::-1]
  tmp_s = tmp_r [idx]
  A_s   = A_r   [idx]
  #=== accumulated prec
  P_acum= 0 # count
  for i in range(Y*X):  # 70%
    P_acum = P_acum + tmp_s[i]
    if P_acum > 0.7 * np.sum(tmp_s):
      crit_i = i+1
      P_check [0, dd] = np.sum(tmp_s [0:i+1]) / np.sum(tmp_s[:])
      A_frac  [0, dd] = np.sum(A_s   [0:i+1]) / np.sum(A_s[:])
      break
  for i in range(crit_i, Y*X):
    P_acum = P_acum + tmp_s[i]
    if P_acum > 0.8 * np.sum(tmp_s):
      crit_i = i+1
      P_check [1, dd] = np.sum(tmp_s [0:i+1]) / np.sum(tmp_s)
      A_frac  [1, dd] = np.sum(A_s   [0:i+1]) / np.sum(A_s)
      break
  for i in range(crit_i, Y*X):
    P_acum = P_acum + tmp_s[i]
    if P_acum > 0.9 * np.sum(tmp_s):
      crit_i = i+1
      P_check [2, dd] = np.sum(tmp_s [0:i+1]) / np.sum(tmp_s)
      A_frac  [2, dd] = np.sum(A_s   [0:i+1]) / np.sum(A_s)
      break
  
# Rad decompose
plt.figure(figsize=(8,3))
plt.plot(np.arange(4,363)*5 ,np.zeros(Nd),'k',linewidth=0.5)
h1= plt.plot(np.arange(4,363)*5 , A_frac[2,:],'r',linewidth=2, label= '90%')
h1= plt.plot(np.arange(4,363)*5 , A_frac[1,:],'k',linewidth=2, label= '80%')
h1= plt.plot(np.arange(4,363)*5 , A_frac[0,:],'b',linewidth=2, label= '70%')
ax=plt.gca()
ax.set_xlim(20, 363*5)
ax.set_ylim(0, 0.6)
plt.legend(loc='upper right', fontsize=8)
picsave= filesave+ 'A_frac.png'
plt.savefig(picsave)

# Rad decompose
plt.figure(figsize=(8,3))
plt.plot(np.arange(4,363)*5 ,np.zeros(Nd),'k',linewidth=0.5)
h1= plt.plot(np.arange(4,363)*5 , P_check[2,:],'r',linewidth=2, label= '90%')
h1= plt.plot(np.arange(4,363)*5 , P_check[1,:],'k',linewidth=2, label= '80%')
h1= plt.plot(np.arange(4,363)*5 , P_check[0,:],'b',linewidth=2, label= '70%')
ax=plt.gca()
ax.set_xlim(20, 363*5)
ax.set_ylim(0, 1)
plt.legend(loc='upper right', fontsize=8)
picsave= filesave+ 'P_check.png'
plt.savefig(picsave)
#plt.show()
