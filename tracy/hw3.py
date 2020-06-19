#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 00:47:31 2020

@author: cloud
"""

import netCDF4 as nc
#import math
import numpy as np
#import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm

"""
file   = "/data/dadm1/reanalysis/ECMWF/ITM/daily/Q/daily_interim_Q_2015.nc"
root   = nc.Dataset(file)
lon    = root.variables['longitude'][:]
lat    = root.variables['latitude'][:]
lev    = root.variables['level'][:]

g = 9.8

CWV_time = np.zeros((18,366))
s = np.zeros((18))
y = 0


#for i in range(1998,2000):
for i in range(1998,2016):
    file   = "/data/dadm1/reanalysis/ECMWF/ITM/daily/Q/daily_interim_Q_" + str(i) + ".nc"
    root   = nc.Dataset(file)
    q      = root.variables['q'][:,:,80:160,:]
    s[y] = np.size(q, axis=0)
    n    = np.size(q, axis=0)
    print(i)

    CWV = np.zeros((n,37,80,480))
    for l in range(36):
        CWV[:,l,:,:] = 0.5*(q[:,l,:,:]+q[:,l+1,:,:])*100*(lev[l]-lev[l+1])/g
        #print(l)
    CWV[:,36,:,:] = q[:,36,:,:]*lev[36]
    CIWV = np.sum(CWV[:], axis=1)
    CIWV_m = CIWV
    CIWV_m[(CIWV_m<40)]  = 0
    CIWV_m[(CIWV_m>=40)]  = 1
    print(l)
    
    for t in range (n):
        CWV_time[y,t] = 100*np.sum(np.sum(CIWV[t,:,:], axis=1), axis=0)/(80*480)
    y = y + 1


p = np.sum(s, axis=0)
CWV_t = np.zeros((int(p)))
CWV_t[0:365] = CWV_time[0,0:365]
for i in range(1,18):
    s1 = int(s[i-1])
    s2 = int(np.sum(s[0:i+1], axis=0))
    s3 = int(s[i])
    CWV_t[s2-s3:s2] = CWV_time[i,0:s3]

x = np.linspace(1,6574,6574)
plt.plot(x,CWV_t)
plt.xlim([0,6574])
plt.xlabel('Time (day)')
plt.ylabel('Fraction (%)')
plt.title('Areal Fraction of Dry Region\n(CWV<40mm, 1998~2015)')
plt.savefig('DryRegion.png')
"""

# practice

file   = "/data/dadm1/reanalysis/ECMWF/ITM/daily/Q/daily_interim_Q_2015.nc"
root   = nc.Dataset(file)
lon    = root.variables['longitude'][:]
lat    = root.variables['latitude'][:]
lev    = root.variables['level'][:]
q      = root.variables['q'][:]

lon_tro = lon[:]
lat_tro = lat[80:160]
q       = q[:,:,80:160,:]

g = 9.8

CWV = np.zeros((365,37,80,480))
for l in range(36):
    CWV[:,l,:,:] = 0.5*(q[:,l,:,:]+q[:,l+1,:,:])*100*(lev[l]-lev[l+1])/g
    #print(l)
CWV[:,36,:,:] = q[:,36,:,:]*lev[36]
CIWV = np.sum(CWV[:], axis=1)

CIWV_m = CIWV
CIWV_m[(CIWV_m<40)]  = 0
CIWV_m[(CIWV_m>=40)]  = 1

CWV_time = np.zeros((365))
for t in range (365):
    CWV_time[t] = 100*np.sum(np.sum(CIWV[t,:,:], axis=1), axis=0)/(80*480)




# plot CIWV for testing

m = Basemap(projection='cyl',llcrnrlat=-30, llcrnrlon=0,urcrnrlat=30, urcrnrlon=360)
m.drawcoastlines()
parallels = np.arange(-30,31,15)
m.drawparallels(parallels) # draw parallels 畫緯度線
meridians = np.arange(0,360,60)
m.drawmeridians(meridians) # draw meridians 畫經度線
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=11)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=11)
LG,LT = np.meshgrid(lon_tro,lat_tro)
cx,cy =m(LG,LT)
CS = m.pcolormesh(cx,cy,CIWV[0,:,:])
plt.colorbar(CS,orientation='vertical')

plt.show()

