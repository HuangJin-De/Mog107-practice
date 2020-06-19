#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 20:29:49 2020

@author: cloud
"""

import netCDF4 as nc
#import math
import numpy as np
#import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm

file   = "/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc"
root   = nc.Dataset(file)
#print(root.variables['FLNT'])
lon =  root.variables['lon'][:]
lat =  root.variables['lat'][:]
"""
OLR = np.zeros((3650,96,144))

earstr  = ['01','02','03','04','05','06','07','08','09','10']
monstr   = ['01','02','03','04','05','06','07','08','09','10','11','12']
day      = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
daystr   = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'\
            ,'17','18','19','20','21','22','23','24','25','26','27','28','29','30','31']

cal = 0
yy = 0
for yn in range (10):
    yy = yy + 1
    print(yy)
    mm = 0
    for mn in range (12):
        mm = mm + 1
        print(mm)
        for dn in range (1,day[mn]+1):
            file   = "/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.00" + '{:0>2}'.format(yy) + "-" + '{:0>2}'.format(mm) + "-" + '{:0>2}'.format(dn) + "-00000.nc"
            root   = nc.Dataset(file)
            OLR[cal,:,:] = np.mean(root.variables['FLNT'][:], axis=0)
            cal = cal + 1
"""

"""
CRH = np.zeros((3650,96,144))

earstr  = ['01','02','03','04','05','06','07','08','09','10']
monstr   = ['01','02','03','04','05','06','07','08','09','10','11','12']
day      = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
daystr   = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'\
            ,'17','18','19','20','21','22','23','24','25','26','27','28','29','30','31']

cal = 0
yy = 0
for yn in range (10):
    yy = yy + 1
    print(yy)
    mm = 0
    for mn in range (12):
        mm = mm + 1
        print(mm)
        for dn in range (1,day[mn]+1):
            file   = "/data/cloud/ucantseeme/hw06/CRH/CRH_" + '{:0>2}'.format(yy) + "-" + '{:0>2}'.format(mm) + "-" + '{:0>2}'.format(dn) + ".nc"
            root   = nc.Dataset(file)
            CRH[cal,:,:] = np.mean(root.variables['cwv'][:], axis=0)
            cal = cal + 1

lon = root.variables['lon'][:]
lat = root.variables['lat'][:]
"""


CorrCO = np.zeros((96,144))
for i in range(96):
    for j in range(144):
        CorrmatrixCRHOLR = np.corrcoef(CRH[:,i,j], OLR[:,i,j])
        CorrCO[i,j] = CorrmatrixCRHOLR[0,1]
     

m = Basemap(projection='cyl',llcrnrlat=-90, llcrnrlon=0,urcrnrlat=90, urcrnrlon=360.1)
m.drawcoastlines()
parallels = np.arange(-90,91,30)
m.drawparallels(parallels) # draw parallels 畫緯度線
meridians = np.arange(0,361,60)
m.drawmeridians(meridians) # draw meridians 畫經度線
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=11)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=11)
LG,LT = np.meshgrid(lon,lat)
cx,cy =m(LG,LT)
CS = m.pcolormesh(cx,cy,CorrCO,cmap = ('RdYlBu_r'),shading='flat')
plt.colorbar(CS,orientation='vertical')
CS.set_clim(-0.8,0.8)
plt.title('Correlation between CRH and OLR')
plt.savefig('corrCRHOLR.png')



