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
var_names= np.array(['MSE','V'])

mapfile= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(mapfile,'r') as D:
    lat=D.variables['lat'][:]
    lon=D.variables['lon'][:]
    ilev= D.variables['ilev'][:]
    lev= D.variables['lev'][:]
    maxY= lat.size
    maxX= lon.size
    maxZ = lev.size

PS_file= data_path+'PS.pk'
with open(PS_file, 'rb') as f:
  PS= pickle.load(f)

maxYr, maxD= PS.shape[0:2]

#=== dx
#weight= np.cos(lat*np.pi/180)
R= 6371*1000
#A= 510100000 * 10**6
#A_w= A * weight / np.sum(weight)
dx= 2*np.pi*R* np.cos(lat*np.pi/180) / lon.size
dx3= dx.reshape((1,1,1, dx.size, 1))
dx3= np.tile(dx3, [maxYr, maxD, maxZ, 1, maxX])

#=== dp
dp= np.diff(ilev)
dp3= dp.reshape((1, 1, maxZ, 1, 1))
dp3= np.tile(dp3, [maxYr, maxD, 1, maxY, maxX])

PS= PS.reshape((maxYr, maxD, 1, maxY, maxX))
PS= np.tile(PS, [1, 1, maxZ, 1, 1])
dp3= dp3/1000 * PS # Pa
g= 9.8 # m/s2
print('Pressure field done !!! ')

#=== load data and plot
data_all={}
for vv in np.arange(0,var_names.shape[0]):
  var_name= var_names[vv]
  print(var_name)
  data_file= data_path+var_name+'.pk'
  with open(data_file, 'rb') as f:
    data= pickle.load(f)
    data_all[var_name]= data

print('Data loaded')

#=== Net energy fluxes
V= data_all['V']
MF= np.sum(np.sum(V*dx3*dp3,axis=4),axis=2) / np.sum(np.sum(dx3*dp3,axis=4),axis=2)
MF= MF.reshape(( maxYr, maxD, 1, maxY, 1))
MF= np.tile(MF, [1, 1, maxZ, 1, maxX])
V_cor= V - MF

VT= V * data_all['MSE']
VT= VT * dx3 * dp3 / g # W

filesave= data_path+ 'VMSE_cor.pk'
with open(filesave, 'wb') as fs:
  pickle.dump(VT, fs, protocol = 4)

VTz= np.nansum(VT, axis=4)
VTz= np.nansum(VTz, axis=2)

filesave= data_path+ 'VMSE_cor_inte.pk'
with open(filesave, 'wb') as fs:
  pickle.dump(VTz, fs, protocol = 4)
