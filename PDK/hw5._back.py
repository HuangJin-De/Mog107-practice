#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 18:11:07 2020

@author: cloud
"""


import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import time


from io import StringIO

m = Basemap(projection='cyl',llcrnrlat=-30,urcrnrlat=30,llcrnrlon=-180,urcrnrlon=180)
parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)

#############################################################################################
# generate column water vapor data
############################################################################################## 
years = np.arange(1999,2015)
path = "/data/dadm1/reanalysis/ECMWF/ITM/daily"

nc  = Dataset(path+'/Q/daily_interim_Q_2001.nc')
q = np.roll(nc.variables['q'][:,0:27,80:160,:],240,axis=3)

tt = np.arange(0,20)
zz = np.arange(0,27)
#yy = np.arange(80,161)

lon_nc=nc.variables['longitude'][:]
lat_nc=nc.variables['latitude'][80:160]
level=nc.variables['level'][zz]

level3d = np.ones((level.size,80,480))*level.reshape(level.size,1,1)
level4d = np.zeros((365,level.size,80,480))
for jjj in range(0,365):
	level3d = np.ones((level.size,80,480))*level.reshape(level.size,1,1)
	level4d[jjj,:,:,:] = level3d 

ciwv = -np.trapz(q,x=level4d*100,axis=1)/9.8
ciwv_lon,ciwv_lat = np.meshgrid(lon_nc,lat_nc,indexing='ij')
AAA
##################################################################################
plt.figure(6)
plt.pcolor(ciwv[0,:,:])
#plt.show()

path2 = '/data/dadm1/obs/TRMM/TRMM3B42'
nc2 = Dataset(path2+'/3B42.1998.3hr.nc')


lat1d = nc2.variables['latitude'][:]
lon1d = nc2.variables['longitude'][:]
print(lon1d)





f = open("/data/cloud/ucantseeme/hw03/cloud_info_daily/cloud_2001.txt","r", encoding='utf-8')
 
raw = f.read()


data =  StringIO(raw)
#aa = np.genfromtxt(data)
data_array = np.loadtxt(data)

#######################################
# calculating CWV

############
# take out specific lat,lon cloud

print(data_array.shape)

cld_array = np.zeros((365,6,36))
ciwv_array = np.zeros((365,6,36))
#cld_array = np.zeros((75,10,36))
#cld_array = np.zeros((10,10,36))

#plt.pcolor(lon2d.T,lat2d.T,.T)
for times in range(1,2921,8): #2921
	for lat in range(30,90,10):
		for lon in range(10,370,10):
			cld_lat_lon = data_array[ ( data_array[:,4] >= (lat-10)*4) & (data_array[:,4] <= lat*4)\
			& (data_array[:,3] >= (lon-10)*4) & (data_array[:,3] <= (lon)*4) & (data_array[:,0] == times) ,:]
            			
			#print(cld_lat_lon.shape)
			#print(cld_lat_lon)
			###########
			# work out SCAI
			# in one region
			num = int(cld_lat_lon.size/5)
			d_geo_mean = 0
			if (num != 0) & (num  != 1):
				#print(num)
				for num_i in range(0,num):
					#print(d_geo_mean)
					for num_j in range(num_i+1,num):

						#d_geo_mean = d_geo_mean + np.log(np.sqrt(np.square( \
						#cld_lat_lon[num_i,4]-cld_lat_lon[num_j,4])\
						#+np.square(cld_lat_lon[num_i,3]-cld_lat_lon[num_j,3])))
						d_geo_mean = d_geo_mean + np.sqrt(np.square( \
                                                cld_lat_lon[num_i,4]-cld_lat_lon[num_j,4])\
                                                +np.square(cld_lat_lon[num_i,3]-cld_lat_lon[num_j,3]))\
						-np.sqrt(cld_lat_lon[num_i,2]) - np.sqrt(cld_lat_lon[num_j,2])

						#print(d_geo_mean)
						#print('aa',np.sqrt(np.square(cld_lat_lon[num_i,4]-cld_lat_lon[num_j,4])\
				                #+np.square(cld_lat_lon[num_i,3]-cld_lat_lon[num_j,3])))
						#print('SSSS',cld_lat_lon[num_i],cld_lat_lon[num_j])
				#d_geo_mean = 25*nuncompatiblem*np.exp((d_geo_mean/num/(num-1)*2))/360/1000*1000
				d_geo_mean = 25*num*d_geo_mean/num/(num-1)*2/360/1000*1000/8
			else:
				d_geo_mean = np.nan
			#print(d_geo_mean)
			#print(int(times-1),int((lat-30)/10),int((lon-10)/10),d_geo_mean)
			if (d_geo_mean<0):
				d_geo_mean = 0

			cld_array[int((times-1)/8),int((lat-30)/10),int((lon-10)/10)] = d_geo_mean
			lat_index = np.where((lat_nc >= lat-60)&(lat_nc <= lat-50))[0]
			lon_index = np.where((lon_nc >= lon-10)&(lon_nc <= lon))[0]
			ciwv_array[int((times-1)/8),int((lat-30)/10),int((lon-10)/10)] = \
			np.mean((ciwv[int((times-1)/8),lat_index])[:,lon_index] < 40)				
			
		#print(int((times-1)/8),int((lat+40)/10),int((lon-10)/10))
print(d_geo_mean)


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

lon1d = np.arange(-180,181,10)
lat1d = np.arange(-30,31,10)
print(lon1d)
LG,LT = np.meshgrid(lon1d,lat1d)
cx,cy =m(LG,LT)
#CS = m.pcolormesh(cx,cy,np.nanmean(cld_array,axis=0),cmap = ('RdYlBu_r'),shading='gouraud')

print(cx[1,:],cy[:,1])
print(np.nanmean(cld_array,axis=0)[1,:])
m.pcolormesh(cx,cy,np.nanmean(cld_array,axis=0),cmap=plt.cm.jet)
plt.colorbar()

plt.figure(5)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)



m.pcolormesh(cx,cy,np.nanmean(ciwv_array,axis=0),cmap=plt.cm.jet)
plt.colorbar()

plt.figure(2)
#plt.plot(np.nanmean(np.nanmean(cld_array,axis=1),axis=1))


C = cld_array[:,0,18]
V = ciwv_array[:,0,18]

CD = C[~np.isnan(C)]
CV = V[~np.isnan(C)]

cor = np.corrcoef(CD,CV)

cor = np.zeros((6,36))
for lat_88 in range(0,6):
    for lon_88 in range(0,36):
        C = cld_array[:,lat_88,lon_88]
        V = ciwv_array[:,lat_88,lon_88]
        
        CD = C[~np.isnan(C)]
        CV = V[~np.isnan(C)]
        
        cor[lat_88,lon_88] = np.corrcoef(CD,CV)[0,1]
#/np.sqrt(np.sum(CD**2))/np.sqrt(np.sum(CV**2))
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
m.pcolormesh(cx,cy,cor,cmap=plt.cm.jet,vmin=-1, vmax=1)
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

