import csv
import numpy as np
import numpy.matlib
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
#import cartopy.crs as ccrs
#import cartopy
#import sys
#sys.path.append('/data/cloud/Goldtimes5/code/')
#from PYfunctions.statistics.calc_areamean import areamean

data_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/'
var_names= np.array(['LHFLX','FLNT','FLNS','FSNTOA','FSNS','SHFLX'])

mapfile= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(mapfile,'r') as D:
    lat=D.variables['lat'][:]
    lon=D.variables['lon'][:]
    maxGY= lat.size
    maxGX= lon.size

weight= np.cos(lat*np.pi/180)
A= 510100000 * 10**6
A_w= A * weight / np.sum(weight)


#=== load data and plot
data_all={}
for vv in np.arange(0,var_names.shape[0]):
  var_name= var_names[vv]
  print(var_name)
  data_file= data_path+var_name+'.pk'
  with open(data_file, 'rb') as f:
    data= pickle.load(f)
    data_all[var_name]= data

#=== Net energy fluxes
TOA_en = data_all['FSNTOA'] - data_all['FLNT'] 
SUR_en = -data_all['FSNS'] + data_all['FLNS'] + data_all['LHFLX'] + data_all['SHFLX']
ATM_en = TOA_en + SUR_en
data_all['TOA_en']= TOA_en
data_all['SUR_en']= SUR_en
data_all['ATM_en']= ATM_en 

#=== Energy transport
ATMz= np.




filesave= data_path+ 'RCE.pk'
with open(filesave, 'wb') as fs:
  pickle.dump(data_all, fs, protocol = 4)
