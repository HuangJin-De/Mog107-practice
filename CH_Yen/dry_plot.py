#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 19:13:11 2020

@author: cloud
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import time
import calendar

path = '/data/cloud/ucantseeme/hw01/'
filename = 'dry_daily.nc'
file = nc.Dataset(path+filename)
DRY = file.variables['dry'][:]*100

filename = 'sub_frac_daily.nc'
file = nc.Dataset(path+filename)
SUB = file.variables['sub_frac'][:]*100

DRY = np.ma.masked_where(DRY == 0, DRY)
SUB = np.ma.masked_where(SUB == 0, SUB)

init_yr = 1979; finl_yr = 2017
yy = finl_yr-init_yr+1
dd = 366
cc = DRY.shape[2] # number of criteria
kk = 37

# In[] plot time series
# daily
Time_dd = np.linspace(init_yr, finl_yr, yy*366)
DRY_dd = np.reshape(DRY, (yy*dd, cc))
SUB_dd = np.reshape(SUB, (yy*dd, kk))

# pentad
DRY_pp = np.zeros((yy,int(365/5), cc))
for p in range(int(365/5)):
	DRY_pp[:,p,:] = np.mean(DRY[:,p*5:(p+1)*5,:], axis=1)
SUB_pp = np.zeros((yy,int(365/5), kk))
for p in range(int(365/5)):
	SUB_pp[:,p,:] = np.mean(SUB[:,p*5:(p+1)*5,:], axis=1)
DRY_pp = np.reshape(DRY_pp, (yy*(p+1), cc))
SUB_pp = np.reshape(SUB_pp, (yy*(p+1), kk))
Time_pp = np.linspace(1979,2017,yy*(p+1))

# monthly
Time_mm = np.linspace(1979,2017,yy*12)
DRY_mm = np.zeros((yy, 12, cc))
for y in range(yy):
    dd = 0
    for m in range(12):
        days = calendar.monthrange(init_yr+y,m+1)[1]
        DRY_mm[y,m,:] = np.mean(DRY[y,dd:dd+days,:], axis=0)
        dd += days
DRY_mm = np.reshape(DRY_mm, (yy*12, cc))

SUB_mm = np.zeros((yy, 12, kk))
for y in range(yy):
    dd = 0
    for m in range(12):
        days = calendar.monthrange(init_yr+y,m+1)[1]
        SUB_mm[y,m,:] = np.mean(SUB[y,dd:dd+days,:], axis=0)
        dd += days
SUB_mm = np.reshape(SUB_mm, (yy*12, kk))

# annually
Time_yy = np.linspace(1979,2017,yy)
DRY_yy = np.zeros((yy,cc))
for y in range(yy):
	DRY_yy[y,:] = np.mean(DRY_mm[y*12:(y+1)*12, :], axis=0)
SUB_yy = np.zeros((yy,kk))
for y in range(yy):
	SUB_yy[y,:] = np.mean(SUB_mm[y*12:(y+1)*12, :], axis=0)

# In[] plot
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(15,5), constrained_layout=True)
ax1.plot(Time_dd, DRY_dd[:,3], Time_pp, DRY_pp[:,3], Time_mm, DRY_mm[:,3], Time_yy, DRY_yy[:,3])
ax1.set_xlim(init_yr, finl_yr)
ax1.set_ylabel('%')
ax1.set_title('Areal fraction of dry region (CWV < 40 mm/day)')

ax2.plot(Time_dd, SUB_dd[:,15], Time_pp, SUB_pp[:,15], Time_mm, SUB_mm[:,15], Time_yy, SUB_yy[:,15])
ax2.legend(['daily','pentad','monthly','annually'], loc='upper right')
ax2.set_xlim(init_yr, finl_yr)
ax2.set_xlabel('year')
ax2.set_ylabel('%')
ax2.set_title('Areal fraction of 500 hPa subsidence')

fig.savefig('dry_sub_series.png', dpi=300)
# In[] plot the standardized time series
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
DRY_dd_std = scaler.fit_transform(DRY_dd)
DRY_dd_std = np.ma.masked_where(DRY_dd_std < -5, DRY_dd_std)
DRY_pp_std = scaler.fit_transform(DRY_pp)
DRY_mm_std = scaler.fit_transform(DRY_mm)
DRY_yy_std = scaler.fit_transform(DRY_yy)
# In[]
SUB_dd_std = scaler.fit_transform(SUB_dd)
SUB_dd_std = np.ma.masked_where(SUB_dd_std < -5, SUB_dd_std)
SUB_pp_std = scaler.fit_transform(SUB_pp)
SUB_mm_std = scaler.fit_transform(SUB_mm)
SUB_yy_std = scaler.fit_transform(SUB_yy)


# In[] plot time series (standardized)
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(15,5), constrained_layout=True)
ax1.plot(Time_dd, DRY_dd_std[:,3], Time_pp, DRY_pp_std[:,3], Time_mm, DRY_mm_std[:,3], Time_yy, DRY_yy_std[:,3])
ax1.set_xlim(init_yr, finl_yr)
ax1.set_ylim([-3,3])
ax1.set_ylabel('%')
ax1.set_title('Areal fraction of dry region (CWV < 40 mm/day)')

ax2.plot(Time_dd, SUB_dd_std[:,15], Time_pp, SUB_pp_std[:,15], Time_mm, SUB_mm_std[:,15], Time_yy, SUB_yy_std[:,15])
ax2.legend(['daily','pentad','monthly','annually'], loc='lower right')
ax2.set_xlim(init_yr, finl_yr)
ax2.set_ylim([-3,3])
ax2.set_ylabel('%')
ax2.set_title('Areal fraction of 500 hPa subsidence')
fig.savefig('dry_sub_series_std.png', dpi=300)

# In[]
DRY_dd = np.zeros((366, cc))
DRY_dd = np.mean(DRY, axis=0)
DRY_pp = np.zeros((int(365/5), cc))
for p in range(int(365/5)):
	DRY_pp[p,:] = np.mean(DRY_dd[p*5:(p+1)*5,:], axis=0)
DRY_mm = np.zeros((12, cc))
dd = 0
for m in range(12):
    days = calendar.monthrange(2017,m+1)[1] # because 2017 has 365 days
    DRY_mm[m,:] = np.mean(DRY_dd[dd:dd+days,:], axis=0)
    dd += days

SUB_dd = np.zeros((366, kk))
SUB_dd = np.mean(SUB, axis=0)
SUB_pp = np.zeros((int(365/5), kk))
for p in range(int(365/5)):
	SUB_pp[p,:] = np.mean(SUB_dd[p*5:(p+1)*5,:], axis=0)
SUB_mm = np.zeros((12, kk))
dd = 0
for m in range(12):
    days = calendar.monthrange(2017,m+1)[1] # because 2017 has 365 days
    SUB_mm[m,:] = np.mean(SUB_dd[dd:dd+days,:], axis=0)
    dd += days

Time_mm = np.linspace(1,12,12)   
Time_pp = np.linspace(1,12,p+1)
Time_dd = np.linspace(1,12,366)

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(15,5), constrained_layout=True)
ax1.plot(Time_dd, DRY_dd[:,3], Time_pp, DRY_pp[:,3], Time_mm, DRY_mm[:,3])
ax1.grid()
ax1.set_xlim([1,12])
ax1.set_xticks(ticks=Time_mm)
ax1.set_ylabel('%')
ax1.set_title('Areal fraction of dry region (CWV < 40 mm/day) (climatology)')

ax2.plot(Time_dd, SUB_dd[:,15], Time_pp, SUB_pp[:,15], Time_mm, SUB_mm[:,15])
ax2.grid()
ax2.legend(['daily','pentad','monthly','annually'], loc='lower right')
ax2.set_xlim([1,12])
ax2.set_xticks(ticks=Time_mm)
ax2.set_ylabel('%')
ax2.set_title('Areal fraction of 500 hPa subsidence (climatology)')
fig.savefig('dry_sub_clmt.png', dpi=300)

# In[] coherence
DRY_dd = np.reshape(DRY, (yy*dd, cc))
SUB_dd = np.reshape(SUB, (yy*dd, kk))

fig, axs = plt.subplots(2,1, figsize=(15,5), constrained_layout=True)
Time = np.linspace(1979,2018,yy*dd)
axs[0].plot(Time, DRY_dd[:,3], Time, SUB_dd[:,15])
axs[0].grid()

axs[1].cohere(DRY_dd[:,3], SUB_dd[:,15],256,1/86400)
# In[] correlation between different layers of subsidence and dry region 
Cor_dd = np.zeros(kk); Cor_pp = np.zeros(kk); Cor_mm = np.zeros(kk); Cor_yy = np.zeros(kk)
for k in range(kk):
    COV = np.corrcoef(DRY_dd[:,1], SUB_dd[:,k], rowvar=False)
    Cor_dd[k] = COV[1,0]
    COV = np.corrcoef(DRY_pp[:,1], SUB_pp[:,k], rowvar=False)
    Cor_pp[k] = COV[1,0]
    COV = np.corrcoef(DRY_mm[:,1], SUB_mm[:,k], rowvar=False)
    Cor_mm[k] = COV[1,0]
    COV = np.corrcoef(DRY_yy[:,1], SUB_yy[:,k], rowvar=False)
    Cor_yy[k] = COV[1,0]

# In[] plot the profiles of correlation coefficients
plt.plot(Cor_dd, Lev, Cor_pp, Lev, Cor_mm, Lev, Cor_yy, Lev)
plt.axvline(x=0., color='k')
plt.legend(['daily', 'pentad', 'monthly', 'anually'])
plt.gca().invert_yaxis()
plt.grid()
plt.title('correlation coefficient between \n areal fractions of subsidence and dry region')
plt.savefig('correcof.png', dpi=300)
# In[] to-do:  varainces of different time scales, correlation with different layers of subsidence, to find which layer has the closest relation between CWV and subsidence