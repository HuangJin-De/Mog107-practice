# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import netCDF4 as nc
#import math
import numpy as np
#import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#import datetime
#import inspect

"""
print('1979')
w_global = np.zeros((37,241,480))
w_in     = np.zeros((37,241,480)) 
file_w   = "/data/dadm1/reanalysis/ECMWF/ITM/daily/W/daily_interim_W_1979.nc"
w_root   = nc.Dataset(file_w)
w_global = np.mean(w_root.variables['w'][:], axis=0)
lon      = w_root.variables['longitude'][:]
lat      = w_root.variables['latitude'][:]

for i in range(1979:2018):
    print(i)
    # open data
    file_w = "/data/dadm1/reanalysis/ECMWF/ITM/daily/W/daily_interim_W_" + str(i) + ".nc"
    w_root = nc.Dataset(file_w)
    w_in = np.mean(w_root.variables['w'][:], axis=0)
    w_global = (w_in + w_global)/2
""" 


###############################################################################
# calculate areal mean MSE (30S ~ 30N) 
"""
T       = np.zeros((39,37,80,480))
q       = np.zeros((39,37,80,480))
z       = np.zeros((39,37,80,480))
T_tro   = np.zeros((37))
q_tro   = np.zeros((37))
z_tro   = np.zeros((37))
MSE_tro = np.zeros((37))

y = 0
for i in range(1979,2018):
    print(i)
    # open data
    file_T = "/data/dadm1/reanalysis/ECMWF/ITM/daily/T/daily_interim_T_" + str(i) + ".nc"
    T_root = nc.Dataset(file_T) 
    file_q = "/data/dadm1/reanalysis/ECMWF/ITM/daily/Q/daily_interim_Q_" + str(i) + ".nc"
    q_root = nc.Dataset(file_q)
    file_z = "/data/dadm1/reanalysis/ECMWF/ITM/daily/Z/daily_interim_Z_" + str(i) + ".nc"
    z_root = nc.Dataset(file_z)
    
    # read variables
    T[y,:,:,:]   = np.mean(T_root.variables['t'][:,:,80:160,:], axis=0)
    q[y,:,:,:]   = np.mean(q_root.variables['q'][:,:,80:160,:], axis=0)
    z[y,:,:,:]   = np.mean(z_root.variables['z'][:,:,80:160,:], axis=0)
    y = y + 1

# deal with variables in region 30S ~ 30N
T_tro[:] = np.mean(np.mean(np.mean(T[:,:,:,:], axis=3), axis=2), axis=0)
q_tro[:] = np.mean(np.mean(np.mean(q[:,:,:,:], axis=3), axis=2), axis=0)
z_tro[:] = np.mean(np.mean(np.mean(z[:,:,:,:], axis=3), axis=2), axis=0)

# define variables
Cp = 1004       # unit = J/kg*k
g  = 9.8        # unit = m/s^2
Lv = 2.5*10**6  # unit = J/kg

# start to calculate MSE
MSE_tro[:] = (Cp*T_tro + z_tro + Lv*q_tro)/Cp
H = z_tro/g/1000
lev = T_root.variables['level'][:]

# save variables
np.savetxt('MSE_climatology.txt',MSE_tro)
np.savetxt('H_climatology.txt',H)
np.savetxt('lev.txt',lev)

# plot mean MSE in the area of 30S to 30N
plt.plot(MSE_tro[0:29],lev[0:29],'k-')
plt.xlabel('MSE (K)',fontsize=12)
plt.ylabel('Pressure (hPa)',fontsize=12)
plt.yscale('log')
plt.gca().invert_yaxis()
plt.title('Mean MSE in Tropical (1979 ~ 2017)',fontsize=12)
#plt.ylim([0,20])
plt.savefig('meanMSE_P.png')
"""

###############################################################################
# calculate zonal vertical profile of MSE
"""
T_ver = np.zeros((160,37))
q_ver = np.zeros((160,37))
z_ver = np.zeros((160,37))

T_ver = np.mean(np.mean(T[:], axis=3), axis=0)
q_ver = np.mean(np.mean(q[:], axis=3), axis=0)
z_ver = np.mean(np.mean(z[:], axis=3), axis=0)
MSE_ver = (Cp*T_ver + z_ver +Lv*q_ver)/Cp
lon_tro = np.arange(-30,30,0.75)
LT,LV = np.meshgrid(lon_tro,H)

plt.contourf(lon_tro,H,MSE_ver,levels= np.linspace(300,420,13),cmap='jet')
plt.ylim([0,20])
plt.colorbar(orientation='vertical')
plt.xlabel('Latitude',fontsize=12)
plt.ylabel('Height (Km)',fontsize=12)
plt.title('Time averaged zonal mean MSE (K)',fontsize=12)
plt.savefig('verticalMSE.png')
"""

###############################################################################
# calculate time series of annual average
"""
T_ann = np.zeros((37,39))
q_ann = np.zeros((37,39))
z_ann = np.zeros((37,39))

T_ann = np.mean(np.mean(T[:], axis=3) ,axis=2).transpose()
q_ann = np.mean(np.mean(q[:], axis=3) ,axis=2).transpose()
z_ann = np.mean(np.mean(z[:], axis=3) ,axis=2).transpose()
MSE_ann = (Cp*T_ann + z_ann + Lv*q_ann)/Cp

year = np.arange(1979,2018)
plt.contourf(year,H,MSE_ann,levels= np.linspace(300,420,13),cmap='jet')
plt.ylim([0,20])
plt.colorbar(orientation='vertical')
plt.xlabel('Time (year)',fontsize=12)
plt.ylabel('Height (Km)',fontsize=12)
plt.title('Annualy averaged vertical profile of MSE (K)',fontsize=12)
plt.savefig('annualyMSE.png')
"""

###############################################################################
# calculate annualy anomaly MSE
"""
MSE_ano = np.zeros((37,39))
for L in range(0,37):
    MSE_ano[L,:] = MSE_ann[L,:] - MSE_tro[L]

plt.contourf(year,H,MSE_ano,levels= np.linspace(-3,3,7),cmap='coolwarm')
plt.ylim([0,20])
plt.colorbar(orientation='vertical')
plt.xlabel('Time (year)',fontsize=12)
plt.ylabel('Height (Km)',fontsize=12)
plt.title('Annualy Anomaly of MSE (K)',fontsize=12)
plt.savefig('anomalyMSE.png')
"""

###############################################################################
"""
for i in range(1979,2018):
    print(i)
    # open data
    file_T = "/data/dadm1/reanalysis/ECMWF/ITM/daily/T/daily_interim_T_" + str(i) + ".nc"
    T_root = nc.Dataset(file_T) 
    file_q = "/data/dadm1/reanalysis/ECMWF/ITM/daily/Q/daily_interim_Q_" + str(i) + ".nc"
    q_root = nc.Dataset(file_q)
    file_z = "/data/dadm1/reanalysis/ECMWF/ITM/daily/Z/daily_interim_Z_" + str(i) + ".nc"
    z_root = nc.Dataset(file_z)
    
    # read variables
    
    T[y,:,:,:]   = np.mean(T_root.variables['t'][:,:,80:160,:], axis=0)
    q[y,:,:,:]   = np.mean(q_root.variables['q'][:,:,80:160,:], axis=0)
    z[y,:,:,:]   = np.mean(z_root.variables['z'][:,:,80:160,:], axis=0)
"""




###############################################################################
# BELOW IS FOR PRACTICING
"""
# open data
T_root = nc.Dataset("/data/dadm1/reanalysis/ECMWF/ITM/daily/T/daily_interim_T_2017.nc") 
q_root = nc.Dataset("/data/dadm1/reanalysis/ECMWF/ITM/daily/Q/daily_interim_Q_2017.nc") 
z_root = nc.Dataset("/data/dadm1/reanalysis/ECMWF/ITM/daily/Z/daily_interim_Z_2017.nc") 


# read variables
lon = T_root.variables['longitude'][:]
lat = T_root.variables['latitude'][:]
lev = T_root.variables['level'][:]
T   = T_root.variables['t'][:]
q   = q_root.variables['q'][:]
z   = z_root.variables['z'][:]


# define variables
Cp = 1004       # unit = J/kg*k
g  = 9.8        # unit = m/s^2
Lv = 2.5*10**6  # unit = J/kg


# define the region from 30S to 30N
T_tro = np.mean(np.mean(np.mean(T[:,:,80:160,:], axis=0), axis=2), axis=1)
q_tro = np.mean(np.mean(np.mean(q[:,:,80:160,:], axis=0), axis=2), axis=1)
z_tro = np.mean(np.mean(np.mean(z[:,:,80:160,:], axis=0), axis=2), axis=1)


# start to calculate MSE
MSE = (Cp*T_tro + z_tro + Lv*q_tro)/Cp
H = np.mean(np.mean(np.mean(z[:], axis=0), axis=2), axis=1)/g/1000


# plot mean MSE in the area of 30S to 30N
plt.plot(MSE[1:28],H[1:28],'k-')
plt.xlabel('MSE (K)',fontsize=12)
plt.ylabel('Height (Km)',fontsize=12)
plt.title('Areal Mean MSE in Tropical (2017)',fontsize=12)
plt.ylim([0,20])
plt.savefig('2017MSE.png')
"""

