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

data_path= '/data/cloud/Goldtimes5/data/SPCAM/CPL64/'
filesave_path= ''

#=== parameters
Lv= 2.477*10**6 # J/kg
rho= 1000 # kg/m3

#=== load data
GCM_temp= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(GCM_temp,'r') as D:
    lat=D.variables['lat'][:]
    lon=D.variables['lon'][:]

lat= lat[40:64]
lon= lon[24:73]
maxGY= lat.size
maxGX= lon.size

#=== load data
data_file= data_path+'Cloud_Type.pk'
with open(data_file, 'rb') as f:
    Cloud_data= pickle.load(f)

data_file= data_path+'PRECT.pk'
with open(data_file, 'rb') as f:
    data= pickle.load(f)
    data= data * Lv * rho # m/s to W/m2

data_file= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/CRM_region_mean/PRECT.pk'
with open(data_file, 'rb') as f:
    prect_mean= pickle.load(f)

data_sm= prect_mean.reshape(prect_mean.size)
data_tm= data.reshape((data.shape[0]*data.shape[1], data.shape[2], data.shape[3]))
data_tm= np.mean(data_tm, axis=0)

#levels= np.zeros((6,10))
#levels[0,:]= np.linspace(0,8  *.9,10);
#levels[1,:]= np.linspace(0,3  *.9,10);
#levels[2,:]= np.linspace(0,3  *.9,10);
#levels[3,:]= np.linspace(0,1  *.9,10);
#levels[4,:]= np.linspace(0,1  *.9,10);
#levels[5,:]= np.linspace(0,1  *.9,10);
level      = np.arange  (1e7, 1.1e8, 1e7);
level_p    = np.arange  (-225, 325, 50);
data_tm[ data_tm< level_p [0] ] = level_p [0]
data_tm[ data_tm> level_p [-1]]= level_p [-1]

cloud_type=['Low','Mid','High','SC','DC','MH']
for ct in np.arange(0,6):
  ctype= cloud_type[ct]
  cdata= Cloud_data['rain'][:,:,:,:,ct]
#  level= levels[ct,:]
  cdata_sm= cdata.reshape((cdata.shape[0], cdata.shape[1], cdata.shape[2]* cdata.shape[3]))
  cdata_sm= np.sum(cdata_sm, axis=2)
  cdata_sm= cdata_sm.reshape(cdata_sm.size)
  cdata_tm= cdata.reshape((cdata.shape[0]* cdata.shape[1], cdata.shape[2], cdata.shape[3]))
  cdata_tm= np.mean(cdata_tm, axis=0)
  cdata_tm[(cdata_tm>level[-1])] = level[-1]
  cdata_tm[(cdata_tm<level[0])] = level[0]

  fig= plt.figure(figsize=(8,3))
  ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
  ax.coastlines()
  cs1= ax.contourf(lon,lat,cdata_tm,levels=level, cmap='RdYlBu_r')
  cs2= ax.contour(lon,lat,data_tm -200 ,levels=level_p, colors= 'k')
  cs1.set_clim(level[0], level[-1])
  plt.colorbar(cs1, orientation= 'horizontal', shrink=0.6).set_ticks(level)
  ax.set_xlim(60, 180)
  ax.set_ylim(-15, 30)
  picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Cloud_type_rain/Spatial_'+ctype+'.png'
  plt.savefig(picsave)
  plt.show()

#  plt.figure(figsize=(8,3))
#  ax = plt.axes(projection=ccrs.PlateCarree())
#  ax.coastlines()
#  cs1= plt.contourf(lon,lat,cdata_tm,levels=level, cmap='RdYlBu_r')
#  plt.clim(level[0], level[-1])
#  plt.colorbar(cs1, orientation= 'horizontal', shrink=0.6).set_ticks(level)
#  ax.set_xlim(60, 180)
#  ax.set_ylim(-15, 30)
#  picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Cloud_type/Spatial_'+ctype+'.png'
#  plt.savefig(picsave)
#  plt.show()
#  cor= np.corrcoef(data_sm, cdata_sm)[0,1]
#  print(ctype,', ',str(cor))
#  plt.figure(figsize=(8,3))
#  plt.plot(data_sm, cdata_sm, 'k.', linestyle='none')
#  picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Cloud_type/PRECT_'+ctype+'.png'
#  plt.savefig(picsave)
#  plt.show()
#  plt.close()

#plt.figure(figsize=(8,3))
#time= np.arange(1,6,1/365)
#plt.plot(time, np.zeros(time.size), 'k', linewidth=0.5)
#colors= np.array([[27,158,119],[102,166,30],[230,171,2],[217,95,2],[231,41,138],[117,112,179]])/256
#lgnd=['Low','Mid','High','SC','DC','Something']
#for ct in np.arange(0,6):
#  ctype= cloud_type[ct]
#  color= colors[ct,:]
#  cdata= Cloud_data['rain'][:,:,:,:,ct]
#  cdata_sm= cdata.reshape((cdata.shape[0], cdata.shape[1], cdata.shape[2]* cdata.shape[3]))
#  cdata_sm= np.sum(cdata_sm, axis=2)
#  cdata_sm= cdata_sm.reshape(cdata_sm.size)
#  plt.plot(time, cdata_sm, color=color, linewidth=1, label= lgnd[ct])
#cdata= Cloud_data['rain']
#cdata= np.sum(cdata,4)
#cdata_sm= cdata.reshape((cdata.shape[0], cdata.shape[1], cdata.shape[2]* cdata.shape[3]))
#cdata_sm= np.sum(cdata_sm, axis=2)
#cdata_sm= cdata_sm.reshape(cdata_sm.size)
#plt.plot(time, cdata_sm, 'k', linewidth=1, label='Total')
#plt.legend(fontsize=12)
#picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Cloud_type/Rain_cloud_type.png'
#plt.savefig(picsave)
#plt.show()


