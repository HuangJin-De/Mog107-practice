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
filesave  = '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/Fast/P_sort/'

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

vars= ['PRECT', 'CAPE', 'MSE', 'OMEGA500']
#=== load data
data= {}
for vv in vars:
  data_file= data_path_map+vv+'.pk'
  with open(data_file, 'rb') as f:
    tmp= pickle.load(f)
    if vv == 'MSE':
      tmp = tmp [:, :, -1, :, :]
    tmp= tmp.reshape(tmp.shape[0]*tmp.shape[1], tmp.shape[2], tmp.shape[3])
    tmp= panta(tmp)
    tmp_smooth= running_mean(tmp,6)
    tmp_f = tmp - tmp_smooth;  
    tmp_f = tmp_f[:, lat_mask, :] # tropics
    tmp_f = tmp_f[3:-3,:,:] # remove nans due to running mean
    tmp_f=  tmp_f.reshape(tmp_f.shape[0], tmp_f.shape[1]*tmp_f.shape[2]) # data [time, space]
  data[vv]= tmp_f

#=== Prepare precipitation matrix
data['PRECT'] = data['PRECT'] * Lv * rho # convert to W/m2
P= data['PRECT']
Nd, Ns= P.shape
P_trop  = np.zeros(Nd) # tropical mean
for dd in range(Nd):
  tmp = P[dd, :]
  nan_mask = ~np.isnan(tmp)
  P_trop [dd]=  np.sum(tmp[nan_mask] * A_r[nan_mask] ) /np.sum(A_r[nan_mask]) #tropical mean

P_trop = P_trop - np.mean(P_trop)
P_pos_mask  = (P_trop >  5)
P_neg_mask  = (P_trop < -5)
P_pos = P [P_pos_mask, :]
P_neg = P [P_neg_mask, :]

#=== P axis
N= 10 # number of intervals
inters   = np.zeros(N+1)
inters[0] = 0
inters[N] = Ns+1
for i in range(1,N):
  inters [i] = round(Ns *i / N)

inters= np.array(inters, dtype='int')

#=== time series for prec
data_pos={}
data_neg={}
for vv in vars:
  print(vv)
  tmp   = data [vv]
  tmp_pos = tmp[P_pos_mask, :]
  tmp_neg = tmp[P_neg_mask, :]
  # sorted array
  tmp_pos_s = np.zeros([tmp_pos.shape[0], N]) 
  tmp_neg_s = np.zeros([tmp_neg.shape[0], N])
  #== For P > 5
  for dd in range(tmp_pos.shape[0]):
    tmp2   = tmp_pos [dd, :]
    P_tmp2 =   P_pos [dd, :]
    idx    = np.argsort(P_tmp2)
    tmp2_s = tmp2[idx]
    for nn in range(N):
      tmp_pos_s [dd, nn] = np.nanmean( tmp2_s [ inters[nn]:inters[nn+1] ] )
  #== For P < -5
  for dd in range(tmp_neg.shape[0]):
    tmp2   = tmp_neg [dd, :]
    P_tmp2 =   P_neg [dd, :]
    idx    = np.argsort(P_tmp2)
    tmp2_s = tmp2[idx]
    for nn in range(N):
      tmp_neg_s [dd, nn] = np.nanmean( tmp2_s [ inters[nn]:inters[nn+1] ] )
#  data_pos[vv] = tmp_pos_s
#  data_neg[vv] = tmp_pos_s
  pos_mean = np.mean(tmp_pos_s, axis=0)
  neg_mean = np.mean(tmp_neg_s, axis=0)
  pos_std  = np.std (tmp_pos_s, axis=0, ddof=1)
  neg_std  = np.std (tmp_pos_s, axis=0, ddof=1)
  #=== plot
  fig= plt.figure(figsize=(8,3))
  plt.plot([0,100] , [0,0],'k',linewidth=0.5)
  plt.plot(np.arange(5,105,10), pos_mean+pos_std,'r--', linewidth=0.5)
  plt.plot(np.arange(5,105,10), pos_mean-pos_std,'r--', linewidth=0.5)
  plt.plot(np.arange(5,105,10), pos_mean     ,'ro', linewidth=1.5, label= 'Positive Prec.')
  plt.plot(np.arange(5,105,10), neg_mean+neg_std,'b--', linewidth=0.5)
  plt.plot(np.arange(5,105,10), neg_mean-neg_std,'b--', linewidth=0.5)
  plt.plot(np.arange(5,105,10), neg_mean     ,'bo', linewidth=1.5, label= 'Negative Prec.')
  ax=plt.gca()
  ax.set_xlim(0, 100)
  ax.set_xticks(range(0,110,10))
#  ax.set_ylim(0, 0.6)
  plt.legend(loc='upper left', fontsize=10)
  picsave= filesave+ vv + '.png'
  plt.savefig(picsave)
  plt.show()
