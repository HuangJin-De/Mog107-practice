#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:49:15 2020

@author: cloud
"""
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import time
from io import StringIO


from read_file import read_file

m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360)
parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)



import glob

path2 = '/data/cloud/ucantseeme/hw06/CRH'


flist2 = glob.glob(path2+'/CRH_01*.nc')
CRH_year = np.zeros((24*len(flist2),96,144))*np.nan
f_idx = 0
for filename in flist2:
    nc2 = Dataset(filename)
    CRH = nc2.variables['cwv'][:]
    CRH_year[24*f_idx:24*(f_idx+1),:,:] = CRH
    f_idx = f_idx + 1

############### CRH plot
lat1d = nc2.variables['lat'][:]
lon1d = nc2.variables['lon'][:]
CRH_std = np.std(CRH_year,axis=0)


LG,LT = np.meshgrid(lon1d,lat1d)
cx,cy =m(LG,LT)

plt.figure(1)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt1 = m.pcolormesh(cx,cy,CRH_std,cmap=plt.cm.jet)
plt.colorbar()

plt.figure(2)

m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt1 = m.pcolormesh(cx,cy,np.mean(CRH_year,axis=0),cmap=plt.cm.jet)
plt.colorbar()



######################### OLR plot

path = '/data/dadm1/model_output/SPCAM/CTRL64'
flist = glob.glob(path+'/CTRL64.cam.h0.0001*-00000.nc')
FLNT_year = np.zeros((24*len(flist),96,144))*np.nan

f_idx =0

for filename in flist:
    nc = Dataset(filename)
    FLNT = nc.variables['FLNT'][:]
    FLNT_year[24*f_idx:24*(f_idx+1),:,:] = FLNT
    f_idx = f_idx + 1

lat1d = nc.variables['lat'][:]
lon1d = nc.variables['lon'][:]

FLNT_std = np.std(FLNT_year,axis=0)


LG,LT = np.meshgrid(lon1d,lat1d)
cx,cy =m(LG,LT)


plt.figure(3)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt1 = m.pcolormesh(cx,cy,FLNT_std,cmap=plt.cm.jet)
plt.colorbar()


plt.figure(4)

m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt1 = m.pcolormesh(cx,cy,np.mean(FLNT_year,axis=0),cmap=plt.cm.jet)
plt.colorbar()



### conbine CRH OLR plot
#plt.figure(5)
#plt.plot(CRH_std[48-16:48+16,:],FLNT_std[48-16:48+16,:],'b.')
## CRH frequency plot
FLNT_disc = np.arange(0,400,10)
CRH_disc = np.arange(0,100,10)
X,Y = np.meshgrid(FLNT_disc,CRH_disc,indexing='ij')

disc_array = np.zeros((FLNT_disc.size,CRH_disc.size))

FLNT_f = np.mean(FLNT_year[:,48-16:48+16,:],axis=0).flatten()
CRH_f = np.mean(CRH_year[:,48-16:48+16,:],axis=0).flatten()
freq_idx = 0
for STD in CRH_f:
    disc_array[int(FLNT_f[freq_idx]/10)-1,int(CRH_f[freq_idx]/10)-1] = disc_array[int(FLNT_f[freq_idx]/10)-1,int(CRH_f[freq_idx]/10)-1] + 1
    freq_idx = freq_idx + 1
plt.figure(5)
plt.pcolor(X,Y,disc_array,cmap='jet')
plt.colorbar()

## CRH frequency plot
CRH_disc = np.arange(0,100,10)
CRH_std_disc = np.arange(0,30,2)
X,Y = np.meshgrid(CRH_disc,CRH_std_disc,indexing='ij')

disc_array = np.zeros((CRH_disc.size,CRH_std_disc.size))

CRH_std_f = CRH_std[48-16:48+16,:].flatten()
CRH_f = np.mean(CRH_year[:,48-16:48+16,:],axis=0).flatten()
freq_idx = 0
for STD in CRH_std_f:
    disc_array[int(CRH_f[freq_idx]/10)-1,int(CRH_std_f[freq_idx]/2)-1] = disc_array[int(CRH_f[freq_idx]/10)-1,int(CRH_std_f[freq_idx]/2)-1] + 1
    freq_idx = freq_idx + 1
plt.figure(10)
plt.pcolor(X,Y,disc_array,cmap='jet')
plt.colorbar()
plt.show()
'''
for filename in flist:
    with open(os.path.join(os.cwd(), filename), 'r') as f: # open in readonly mode

for date in range(1,29):
    if date < 10:
        day = '0' + str(date)
    else:
        day = str(date)
    path = '/data/dadm1/model_output/SPCAM/CTRL64'
    nc = Dataset(path+'/CTRL64.cam.h0.0004-05-' + day + '-00000.nc')
    FLNT = nc.variables['FLNT'][:]
''' 
    




AA
##########################
# reset ciwv to start at -180
'''
ciwv_array_r = np.zeros((365,6,36))
ciwv_array_r[:,:,0:18] = ciwv_array[:,:,18:]
ciwv_array_r[:,:,18:] = ciwv_array[:,:,0:18]
ciwv_array = ciwv_array_r
'''
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

lon1d = np.arange(-180,181,20)
lat1d = np.arange(-30,31,20)
print(lon1d)
LG,LT = np.meshgrid(lon1d,lat1d)
cx,cy =m(LG,LT)
#CS = m.pcolormesh(cx,cy,np.nanmean(cld_array,axis=0),cmap = ('RdYlBu_r'),shading='gouraud')

print(cx[1,:],cy[:,1])
print(np.nanmean(cld_array,axis=0)[1,:])
m.pcolormesh(cx,cy,np.nanmean(cld_array,axis=0),cmap=plt.cm.jet)
plt.colorbar()
plt.show()
AAAA
plt.figure(5)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)



m.pcolormesh(cx,cy,np.nanmean(ciwv_array,axis=0),cmap=plt.cm.jet)
plt.colorbar()

plt.figure(2)
#plt.plot(np.nanmean(np.nanmean(cld_array,axis=1),axis=1))

'''
C = cld_array[:,0,18]
V = ciwv_array[:,0,18]

CD = C[~np.isnan(C)]
CV = V[~np.isnan(C)]

cor = np.corrcoef(CD,CV)
'''
cor = np.zeros((len(years),6,36))
for year in years:
    for lat_88 in range(0,6):
        for lon_88 in range(0,36):
            C = cld_array_all[year-years[0],:,lat_88,lon_88]
            V = ciwv_array_all[year-years[0],:,lat_88,lon_88]
            
            CD = C[~np.isnan(C)]
            CV = V[~np.isnan(C)]
            
            cor[year-years[0],lat_88,lon_88] = np.corrcoef(CD,CV)[0,1]

cor_mean = np.nanmean(cor,axis=0)
#/np.sqrt(np.sum(CD**2))/np.sqrt(np.sum(CV**2))
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
m.pcolormesh(cx,cy,cor_mean,cmap=plt.cm.jet,vmin=-0.5, vmax=0.5)
plt.colorbar()

R1 = np.random.rand(1,200)
R2 = np.random.rand(1,200)

cor2 = np.corrcoef(R1[0,:],R2[0,:])# / np.sqrt(np.sum(R1**2))/np.sqrt(np.sum(R2**2))
#W_mean = W_mean_sum /(years.size+1)

#lat2d,lon2d = np.meshgrid(lat,lon,indexing='ij')
#level3d = np.ones((level.size,lat.size,lon.size))*level.reshape(level.size,1,1)
#cwv = -np.trapz(2260000*q_mean+1005*T_mean+Z_mean,x=level3d*100,axis=0)/9.8/1005

#print(lat2d.size,W_mean.size)


#plt.contourf(lon2d.T,lat2d.T,W_mean.T)
#plt.colorbar(orientation='horizontal')
plt.show()