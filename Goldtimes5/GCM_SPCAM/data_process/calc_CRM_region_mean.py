import csv
import numpy as np
import numpy.matlib
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy
import sys
sys.path.append('/data/cloud/Goldtimes5/code/')
from PYfunctions.statistics.calc_areamean import areamean


data_file= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/RCE.pk'
datasave_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/CRM_region_mean/'
var_names= np.array(['PRECT'])

# map
GCM_temp= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(GCM_temp,'r') as D:
  lat=D.variables['lat'][:] # level in hPa
  lon=D.variables['lon'][:] # level in hPa
maxY= lat.size
maxX= lon.size

Lv= 2.266*10**6 # J/kg
rho= 1000 # kg/m3

#=== load data and plot
data_m={}
with open(data_file, 'rb') as f:
  data= pickle.load(f)

Nyr= data['FLNT'].shape[0]
Nd=  data['FLNT'].shape[1]

for vv in np.arange(0,var_names.shape[0]):
  var_name= var_names[vv]
  print(var_name)
  temp_all= np.zeros((Nyr, Nd))
  for yr in np.arange(0,Nyr):
    for day in np.arange(0,Nd):
      temp=  data[var_name][yr,day,:,:]
      tempm= areamean(temp, lat, lon, [59.5, 180.5, -15, 30])
      temp_all[yr, day]= tempm

  filesave= datasave_path+ var_name +'.pk'
  with open(filesave, 'wb') as fs:
    pickle.dump(temp_all, fs, protocol = 4)
