# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import time

path = '/data/dadm1/reanalysis/ECMWF/ITM/daily/'

# In[] read the longitude and latitude
filename = 'T/daily_interim_T_1979.nc'

file = nc.MFDataset(path+filename)
lon = file.variables['longitude'][:]
lat = file.variables['latitude'][:]
level = file.variables['level'][:]
sth = -30; nth = 30; # the boundary of tropics (30S~30N)
sth_y = np.where(lat==sth); nth_y = np.where(lat==nth); 
sth_y = int(sth_y[0]); nth_y = int(nth_y[0]);


# In[]

years = [str(i) for i in range(1979,2018)]

filename = 'T/daily_interim_T_*.nc'
files = nc.MFDataset(path+filename)
T = files.variables['t']

filename = 'Z/daily_interim_Z_*.nc'
files = nc.MFDataset(path+filename)
Z = files.variables['z']

filename = 'Q/daily_interim_Q_*.nc'
files = nc.MFDataset(path+filename)
Qv = files.variables['q']

# In[] MSE
Cp = 1004; Lv = 2.26E6; g = 9.8
#T_mean = np.mean(T[:,:,nth_y:sth_y,:], axis=(2,3))
MSE = Cp*np.mean(T[:,:,nth_y:sth_y,:], axis=(2,3)) \
    + np.mean(Z[:,:,nth_y:sth_y,:], axis=(2,3)) \
    + Lv*np.mean(Qv[:,:,nth_y:sth_y,:], axis=(2,3))
#end = time.time()
#print(end - start)
#for i in range(0, len(years)): 
#    print(years[i])
#    filename = 'T/daily_interim_T_' + years[i] + '.nc'
#    file = Dataset(path+filename)
#    start = time.time()
#    T = file.variables['t'][:]
##    if i == 0:
##        T_dd = T[:,:,nth_y:sth_y,:]
##    else:
##        np.concatenate((T_dd,T[:,:,nth_y:sth_y,:]))
#    end = time.time()
#    print(end - start)

# In[]
filename = 'CWV_daily.nc'
f=nc.Dataset(filename)
cwv=f.variables['cwv']
cwv2=cwv[,,:,:]