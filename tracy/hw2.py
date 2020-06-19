#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 16:46:41 2020

@author: cloud
"""

"""
GOAL
1. average cloud size
2. average rain rate
3. cloud frequency
4. distribution of certain criteria of cloud size
"""

import netCDF4 as nc
#import math
import numpy as np
#import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm


###############################################################################
"""
file   = "/data/dadm1/obs/TRMM/TRMM3B42/3B42.2015.3hr.nc"
root   = nc.Dataset(file)
lon    = root.variables['longitude'][:]
lat    = root.variables['latitude'][:]


cloud_tro = np.zeros((18,161,241))
L = 0
for i in range(1998,2016):
    print(i)
    file  = "/data/dadm1/obs/TRMM/TRMM3B42size/TRMMsize_3hrs_" + str(i) + ".nc"
    root  = nc.Dataset(file)
    cloudsize = root.variables['objsize'][:]
    cloudtemp = cloudsize[:,160:321,1120:1361]
    cloudtemp[(cloudtemp>0)]  = 1
    cloudtemp[(cloudtemp<=0)] = 0
    cloud_tro[L,:,:] = np.sum(cloudtemp, axis=0)
    L = L + 1


t = 52592
cloud_fre = 100 * np.sum(cloud_tro, axis=0)/t
lon_tro = lon[1120:1361]
lat_tro = lat[160:321]
m = Basemap(projection='cyl',llcrnrlat=-10, llcrnrlon=100,urcrnrlat=30, urcrnrlon=160)
m.drawcoastlines()
parallels = np.arange(-10,31,10)
m.drawparallels(parallels) # draw parallels 畫緯度線
meridians = np.arange(100,161,20)
m.drawmeridians(meridians) # draw meridians 畫經度線
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=11)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=11)
LG,LT = np.meshgrid(lon_tro,lat_tro)
cx,cy =m(LG,LT)
CS = m.pcolormesh(cx,cy,cloud_fre,cmap = ('RdYlBu_r'),shading='gouraud')
plt.colorbar(CS,orientation='vertical',pad=0.08)
plt.title('Cloud Frequency (%)\n(1998 ~ 2015)',fontsize=14)
plt.savefig('CldFre.png')
"""

###############################################################################
# averaged rain rate (mm/hr)
"""
file   = "/data/dadm1/obs/TRMM/TRMM3B42/3B42.2015.3hr.nc"
root   = nc.Dataset(file)
lon    = root.variables['longitude'][:]
lat    = root.variables['latitude'][:]

L = 0
pcp_tro = np.zeros((18,161,241))
for i in range(1998,2016):
    print(i)
    file   = "/data/dadm1/obs/TRMM/TRMM3B42/3B42." + str(i) + ".3hr.nc"
    root   = nc.Dataset(file)
    pcp    = root.variables['pcp'][:,160:321,1120:1361]
    tt     = np.size(pcp,0)
        
    pcp[(pcp==0)] = np.nan
    pcp_tro[L,:,:] = np.nanmean(pcp[:], axis = 0)
    L = L + 1
pcp_timemean = np.nanmean(pcp_tro[:], axis = 0)


lon_tro = lon[1120:1361]
lat_tro = lat[160:321]
m = Basemap(projection='cyl',llcrnrlat=-10, llcrnrlon=100,urcrnrlat=30, urcrnrlon=160)
m.drawcoastlines()
parallels = np.arange(-10,31,10)
m.drawparallels(parallels) # draw parallels 畫緯度線
meridians = np.arange(100,161,20)
m.drawmeridians(meridians) # draw meridians 畫經度線
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=11)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=11)
LG,LT = np.meshgrid(lon_tro,lat_tro)
cx,cy =m(LG,LT)
CS = m.pcolormesh(cx,cy,pcp_timemean,cmap=cm.jet,shading='gouraud')
plt.colorbar(CS,orientation='vertical',pad=0.08)
plt.title('Averaged Rain Rate (mm/hr)\n(1998 ~ 2015)',fontsize=14)
plt.savefig('AveRainRate.png')
"""


###############################################################################
# practice

file   = "/data/dadm1/obs/TRMM/TRMM3B42/3B42.2015.3hr.nc"
root   = nc.Dataset(file)
lon    = root.variables['longitude'][:]
lat    = root.variables['latitude'][:]
pcp    = root.variables['pcp'][:]
# print(root.variables['pcp'])

"""
pcp_tro = pcp[1,160:320,1120:1360]
lon_tro = lon[1120:1360]
lat_tro = lat[160:320]


m = Basemap(projection='cyl',llcrnrlat=-10, llcrnrlon=100,urcrnrlat=30, urcrnrlon=160)
m.drawcoastlines()
LG,LT = np.meshgrid(lon_tro,lat_tro)
cx,cy =m(LG,LT)
CS = m.pcolormesh(cx,cy,pcp_tro,cmap=cm.jet,shading='gouraud')
plt.colorbar(CS,orientation='horizontal')
"""

