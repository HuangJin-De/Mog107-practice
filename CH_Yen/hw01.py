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
import calendar

path = '/data/cloud/ucantseeme/hw01/'

# In[] read MSE
filename = 'mse_daily.nc'
Cp = 1004
file = nc.Dataset(path+filename)
Lev = file.variables['level'][:]
kk  = len(Lev)
init_yr = 1979; finl_yr = 2017
yy = finl_yr-init_yr+1

MSE = file.variables['mse'][:]/Cp
HGT = file.variables['height'][:]
MSE = np.ma.masked_where(MSE < 10, MSE)
HGT = np.ma.masked_where(HGT < 10, HGT)
Hgt = np.mean(np.reshape(HGT,(yy*366, kk)), axis=0)
dd  = 366

# In[]
#trop = int(np.where(Lev==50)[0])

vmin = 320; vmax = 350;
levels = np.linspace(vmin, vmax, (vmax-vmin)/5+1)
cmap = 'coolwarm'
# In[] plot time series
# daily
Time = np.linspace(init_yr, finl_yr, yy*366)
MSE_dd = np.reshape(MSE, (yy*dd, kk))
plt.figure(figsize=(15,5))
plt.contourf(Time, Hgt/1E3, MSE_dd.T, levels=levels, vmin=vmin, vmax=vmax, extend='max', cmap=cmap)
plt.colorbar()
plt.ylim(0, 20)
plt.grid(color='w')
plt.xlabel('year')
plt.ylabel('height (km)')
plt.title('Daily MSE (K)')
plt.savefig('mse_daily.png', dpi=300)
plt.clf()
del MSE_dd

# In[] plot pentad series
MSE_pp = np.zeros((yy,int(365/5), kk))
for p in range(int(365/5)):
	MSE_pp[:,p,:] = np.mean(MSE[:,p*5:(p+1)*5,:], axis=1)
MSE_pp = np.reshape(MSE_pp, (yy*(p+1), kk))
Time = np.linspace(1979,2017,yy*(p+1))
plt.figure(figsize=(15,5))
plt.contourf(Time, Hgt/1E3, MSE_pp.T, levels=levels, vmin=vmin, vmax=vmax, extend='max', cmap=cmap)
plt.colorbar()
plt.ylim(0, 20)
plt.grid(color='w')
plt.xlabel('year')
plt.ylabel('height (km)')
plt.title('Pentad MSE (K)')
plt.savefig('mse_pentad.png', dpi=300)
plt.clf()
del MSE_pp

# In[] plot monthly aeries
MSE_mm = np.zeros((yy,12,kk))
for y in range(yy):
    dd = 0
    for m in range(12):
        days = calendar.monthrange(init_yr+y,m+1)[1]
        MSE_mm[y,m,:] = np.mean(MSE[y,dd:dd+days,:], axis=0)
        dd += days
MSE_mm = np.reshape(MSE_mm, (yy*12, kk))
Time = np.linspace(1979,2017,yy*12)
plt.figure(figsize=(15,5))
plt.contourf(Time, Hgt/1E3, MSE_mm.T, levels=levels, vmin=vmin, vmax=vmax, extend='max', cmap=cmap)
plt.colorbar()
plt.ylim(0, 20)
plt.grid(color='w')
plt.xlabel('year')
plt.ylabel('height (km)')
plt.title('Monthly MSE (K)')
plt.savefig('mse_month.png', dpi=300)
plt.clf()

# In[] plot annual series
MSE_yy = np.zeros((yy,kk))
for y in range(yy):
	MSE_yy[y,:] = np.mean(MSE_mm[y*12:(y+1)*12, :], axis=0)

Time = np.linspace(1979,2017,yy)
plt.figure(figsize=(15,5))
plt.contourf(Time, Hgt/1E3, MSE_yy.T, levels=levels, vmin=vmin, vmax=vmax, extend='max', cmap=cmap)
plt.colorbar()
plt.ylim(0, 20)
plt.grid(color='w')
plt.xlabel('year')
plt.ylabel('height (km)')
plt.title('Annually MSE (K)')
plt.savefig('mse_annual.png', dpi=300)
plt.clf()

# In[] modify the levels of colorbar
levels = np.linspace(vmin, vmax, (vmax-vmin)/2+1)
# In[] plot the climatology of daily MSE
MSE_dd = np.zeros((366, kk))
MSE_dd = np.mean(MSE, axis=0)
Time = np.linspace(1,12,366)
plt.contourf(Time, Hgt/1E3, MSE_dd.T, levels=levels, vmin=vmin, vmax=vmax, extend='max', cmap=cmap)
plt.colorbar()
plt.ylim(0, 20)
plt.grid(color='w')
plt.xlabel('month')
plt.ylabel('height (km)')
plt.title('Daily MSE Climatology (K)')
plt.savefig('mse_daily_clmt.png', dpi=300)
plt.clf()

# In[] plot the climatology of pentad MSE
MSE_pp = np.zeros((int(365/5), kk))
for p in range(int(365/5)):
	MSE_pp[p,:] = np.mean(MSE_dd[p*5:(p+1)*5,:], axis=0)
    
Time = np.linspace(1,12,p+1)
plt.contourf(Time, Hgt/1E3, MSE_pp.T, levels=levels, vmin=vmin, vmax=vmax, extend='max', cmap=cmap)
plt.colorbar()
plt.ylim(0, 20)
plt.grid(color='w')
plt.xlabel('month')
plt.ylabel('height (km)')
plt.title('Pentad MSE Climatology (K)')
plt.savefig('mse_pentad_clmt.png', dpi=300)
plt.clf()


# In[] plot the climatology of monthly MSE
MSE_mm = np.zeros((12, kk))
dd = 0
for m in range(12):
    days = calendar.monthrange(2017,m+1)[1] # because 2017 has 365 days
    MSE_mm[m,:] = np.mean(MSE_dd[dd:dd+days,:], axis=0)
    dd += days
    
Time = np.linspace(1,12,12)
plt.contourf(Time, Hgt/1E3, MSE_mm.T, levels=levels, vmin=vmin, vmax=vmax, extend='max', cmap=cmap)
plt.colorbar()
plt.ylim(0, 20)
plt.grid(color='w')
plt.xlabel('month')
plt.ylabel('height (km)')
plt.title('Monthly MSE Climatology (K)')
plt.savefig('mse_month_clmt.png', dpi=300)
plt.clf()
