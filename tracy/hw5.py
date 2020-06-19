#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 21:09:38 2020

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
file   = "/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc"
root   = nc.Dataset(file)
lon =  root.variables['lon'][:]
lat =  root.variables['lat'][:]
Cp = 1004
Lv = 2.5*10**6

gZ  = np.zeros((10,12,96,144))
CpT = np.zeros((10,12,96,144))
Lvq = np.zeros((10,12,96,144))

earstr  = ['01','02','03','04','05','06','07','08','09','10']
monstr   = ['01','02','03','04','05','06','07','08','09','10','11','12']
day      = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
daystr   = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'\
            ,'17','18','19','20','21','22','23','24','25','26','27','28','29','30','31']
yy = 0
for yn in range (10):
    yy = yy + 1
    print(yy)
    mm = 0
    for mn in range (12):
        mm = mm + 1
        print(mm)
        geop_temp = np.zeros((day[mn],96,144))
        q_temp    = np.zeros((day[mn],96,144))
        T_temp    = np.zeros((day[mn],96,144))
        for dn in range (day[mn]):
            geop_temp[dn,:,:] = np.sum(np.mean(root.variables['Z3'][:], axis=0),axis=0)
            q_temp[dn,:,:]    = np.sum(np.mean(root.variables['Q'][:], axis=0),axis=0)
            T_temp[dn,:,:]    = np.sum(np.mean(root.variables['T'][:], axis=0),axis=0)
        gZ[yn,mn,:,:]  = np.mean(geop_temp, axis=0)
        CpT[yn,mn,:,:] = Cp * np.mean(T_temp, axis=0)
        Lvq[yn,mn,:,:] = Lv * np.mean(q_temp, axis=0)
"""

MSE_g = np.mean(np.mean((gZ + CpT + Lvq)/Cp, axis=0), axis=0)
for j in range(96):
    MSE_global[j,:] = MSE_g[95-j,:]
m = Basemap(projection='cyl',llcrnrlat=-90, llcrnrlon=0,urcrnrlat=90, urcrnrlon=360.5)
m.drawcoastlines()
parallels = np.arange(-90,91,30)
m.drawparallels(parallels) # draw parallels 畫緯度線
meridians = np.arange(0,361,60)
m.drawmeridians(meridians) # draw meridians 畫經度線
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=11)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=11)
LG,LT = np.meshgrid(lon,lat)
cx,cy =m(LG,LT)
CS = m.pcolormesh(cx,cy,MSE_global[:],cmap = ('RdYlBu_r'),shading='flat')
CS.set_clim(6000,7000)
plt.colorbar(CS,orientation='vertical')
plt.title('CPL64 MSE Global Distribution (K)')
plt.savefig('MSE.png')



   
        


# MSE zonal mean
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
lon =  root.variables['lon'][:]
lat =  root.variables['lat'][:]
lev =  root.variables['lev'][:]
Cp = 1004
Lv = 2.5*10**6

gZ  = np.zeros((10,12,26,96))
CpT = np.zeros((10,12,26,96))
Lvq = np.zeros((10,12,26,96))

yearstr  = ['01','02','03','04','05','06','07','08','09','10']
monstr   = ['01','02','03','04','05','06','07','08','09','10','11','12']
day      = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
daystr   = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'\
            ,'17','18','19','20','21','22','23','24','25','26','27','28','29','30','31']
yy = 0
for yn in range (10):
    yy = yy + 1
    print(yy)
    mm = 0
    for mn in range (12):
        mm = mm + 1
        print(mm)
        geop_temp = np.zeros((day[mn],26,96))
        q_temp    = np.zeros((day[mn],26,96))
        T_temp    = np.zeros((day[mn],26,96))
        for dn in range (day[mn]):
            file = nc.Dataset('/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.00'+yearstr[yn]+'-'+monstr[mn]+'-'+daystr[dn]+'-00000.nc')
            geop_temp[dn,:,:] = np.mean(np.mean(root.variables['Z3'][:], axis=0), axis=2)
            q_temp[dn,:,:]    = np.mean(np.mean(root.variables['Q'][:], axis=0), axis=2)
            T_temp[dn,:,:]    = np.mean(np.mean(root.variables['T'][:], axis=0), axis=2)
        gZ[yn,mn,:,:]  = np.mean(geop_temp, axis=0)
        CpT[yn,mn,:,:] = Cp * np.mean(T_temp, axis=0)
        Lvq[yn,mn,:,:] = Lv * np.mean(q_temp, axis=0)
"""
"""
MSE_global_zonal = np.mean(np.mean((CpT + gZ +Lvq)/Cp, axis=0), axis=0)
LL,HH = np.meshgrid(lat,lev)
CS = plt.pcolormesh(LL,HH,MSE_global_zonal,cmap = ('RdYlBu_r'),shading='flat')
#plt.gca().invert_yaxis()
CS.set_clim(200,320)
plt.colorbar(CS,orientation='vertical')
plt.xlabel('Latitude',fontsize=13.5)
plt.ylabel('Sigma Coordinate',fontsize=13.5)
plt.title('CPL64 MSE Zonal Mean',fontsize=13.5)
plt.savefig('MSEzonal.png')
"""
"""
MSE_global_JJA = np.mean(np.mean((CpT[:,5:9,:,:] + gZ[:,5:9,:,:] +Lvq[:,5:9,:,:])/Cp, axis=0), axis=0)
LL,HH = np.meshgrid(lat,lev)
CS = plt.pcolormesh(LL,HH,MSE_global_JJA,cmap = ('RdYlBu_r'),shading='flat')
#plt.gca().invert_yaxis()
CS.set_clim(200,320)
plt.colorbar(CS,orientation='vertical')
plt.xlabel('Latitude',fontsize=13.5)
plt.ylabel('Sigma Coordinate',fontsize=13.5)
plt.title('CPL64 MSE Zonal Mean (JJA)',fontsize=13.5)
plt.savefig('MSEzonal_JJA.png')
"""
"""
MSE_global_JF = np.mean(np.mean((CpT[:,0:2,:,:] + gZ[:,0:2,:,:] +Lvq[:,0:2,:,:])/Cp, axis=0), axis=0)
MSE_global_D = np.mean((CpT[:,11,:,:] + gZ[:,11,:,:] +Lvq[:,11,:,:])/Cp, axis=0)
MSE_global_DJF = (MSE_global_JF + MSE_global_D)/2
LL,HH = np.meshgrid(lat,lev)
CS = plt.pcolormesh(LL,HH,MSE_global_DJF,cmap = ('RdYlBu_r'),shading='flat')
#plt.gca().invert_yaxis()
CS.set_clim(200,320)
plt.colorbar(CS,orientation='vertical')
plt.xlabel('Latitude',fontsize=13.5)
plt.ylabel('Sigma Coordinate',fontsize=13.5)
plt.title('CPL64 MSE Zonal Mean (DJF)',fontsize=13.5)
plt.savefig('MSEzonal_DJF.png')
"""


# below is for practice
"""
file   = "/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc"
root   = nc.Dataset(file)

lon =  root.variables['lon'][:]
lat =  root.variables['lat'][:]
lev =  root.variables['lev'][:]
geop = root.variables['Z3'][:]
q   =  root.variables['Q'][:]
T   =  root.variables['T'][:]
Cp = 1004
Lv = 2.5*10**6

# MSE = (CpT + gZ + Lvq)/Cp
CpT = Cp * np.mean(np.mean(T[:,:,:,:], axis = 0), axis = 2)
gZ = np.mean(np.mean(geop[:,:,:,:], axis = 0), axis = 2)
Lvq = Lv * np.mean(np.mean(q[:,:,:,:], axis = 0), axis = 2)
MSE_global_zonal = (CpT + gZ + Lvq)/Cp
LL,HH = np.meshgrid(lat,lev)
CS = plt.pcolormesh(LL,HH,MSE_global_zonal,cmap = ('RdYlBu_r'),shading='flat')
#plt.gca().invert_yaxis()
CS.set_clim(200,340)
plt.colorbar(CS,orientation='vertical')
"""
