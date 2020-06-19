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

data_path_map= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/'
filesave  = '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/Fast/LargeVar_region/'

#=== parameters
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
A        = A.reshape((A.size, 1))
A        = np.tile(A, [1, len(lon)] )
A_r   = A.reshape(A.size) # Area of each grid (cosine latitude)

vars= ['PRECT']#, 'CAPE','MSE_sur']
#=== load data
data= {}
for vv in vars:
  data_file= data_path_map+vv+'.pk'
  with open(data_file, 'rb') as f:
    tmp= pickle.load(f)
    tmp= tmp.reshape(tmp.shape[0]*tmp.shape[1], tmp.shape[2], tmp.shape[3])
    tmp= panta(tmp)
    tmp_smooth= running_mean(tmp,6)
    tmp_f = tmp - tmp_smooth;  
    tmp_f = tmp_f[:, lat_mask, :] # tropics
    tmp_f = tmp_f[3:-3,:,:] # remove nans due to running mean
  data[vv]= tmp_f

#=== convert unit
P= data['PRECT']
P= P * Lv * rho # convert to W/m2

#=== prepare matrix
Nd, Y, X= P.shape
Pstd= np.std(P, axis=0, ddof=1)
Pstd_mask= (Pstd>120) # large variance regions
Pstdr= Pstd_mask.reshape(Pstd_mask.size)
A_frac=  np.sum(A_r [Pstdr]) / np.sum(A_r)
print(A_frac)

#=== time series for prec
P_trop  = np.zeros(Nd) # tropical mean
P_lv    = np.zeros(Nd) # large-variance region
for dd in range(Nd):
  tmp = P[dd, :, :]
  nan_mask = ~np.isnan(tmp)
  P_trop [dd]=  np.sum(tmp[nan_mask] * A[nan_mask] ) / np.sum(A[nan_mask])
  nan_Pstd = nan_mask & Pstd_mask
  P_lv [dd]=  np.sum(tmp[nan_Pstd] * A[nan_Pstd] ) / np.sum(A[nan_mask])

#=== plot P std
clevel= np.arange(20,200,20) 
Pstd[(Pstd < 120)] = np.nan #clevel[0]
Pstd[(Pstd > clevel[-1])] = clevel[-1]
plt.figure(figsize=(8,6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
cs1= plt.contourf(lon,lat, Pstd, levels=clevel ,cmap='RdYlBu_r')
#plt.clim(,cmax)
plt.colorbar(shrink=0.6, orientation='horizontal').set_ticks(clevel)
ax.set_xlim(-180, 180)
ax.set_ylim(-30, 30)
picsave= filesave+ 'P_std.png'
plt.savefig(picsave)
plt.show()

cor = np.corrcoef(P_trop, P_lv)[0,1]
print(cor)
plt.figure(figsize=(8,3))
plt.plot(np.arange(4, 4+Nd)*5 ,np.zeros(Nd),'k',linewidth=0.5)
plt.plot(np.arange(4, 4+Nd)*5 ,P_trop,'k',linewidth=1, label= 'Precipitation')
plt.plot(np.arange(4, 4+Nd)*5 ,P_lv  ,'r',linewidth=1, label= 'Precipitation (Large variance region)')
ax=plt.gca()
ax.set_xlim(1+365, 366+365)
#ax.set_ylim(50, 150)
plt.legend(loc='upper right', fontsize=8)
picsave= filesave+ 'LVregion_trop_mean_PRECT.png'
plt.savefig(picsave)
plt.show()




#=== 
#for 



