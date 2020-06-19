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

#m = Basemap(projection='cyl',lon_0=180,lat_0=0,llcrnrlat=-30,urcrnrlat=-30)
m = Basemap(projection='cyl',llcrnrlat=-30,urcrnrlat=30,llcrnrlon=-180,urcrnrlon=180)


parallels = np.arange(-90,90,30)
meridians = np.arange(-180,180,60)






#cx,cy =m(LG,LT)
#CS = m.pcolormesh(cx,cy,pcp_timemean,cmap=cm.jet,shading='gouraud')

path = "/data/dadm1/reanalysis/ECMWF/ITM/daily"
nc  = Dataset(path+"/Q/daily_interim_Q_2002.nc")
q = nc.variables['q'][:,0:27,80:160,:]

nc = Dataset(path+'/T/daily_interim_T_1979.nc')
T = nc.variables['t'][:,0:27,80:160,:]

nc = Dataset(path+'/Z/daily_interim_Z_1979.nc')
Z = nc.variables['z'][:,0:27,80:160,:]


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

#print(nc.variables['q'])
ciwv = -np.trapz(q,x=level4d*100,axis=1)/9.8
ciwv_lon,ciwv_lat = np.meshgrid(lon_nc,lat_nc,indexing='ij')

#print(ciwv_lat.shape)



plt.figure(6)
plt.pcolor(ciwv[0,:,:])
#plt.show()

path2 = '/data/dadm1/obs/TRMM/TRMM3B42'
nc2 = Dataset(path2+'/3B42.1998.3hr.nc')


lat1d = nc2.variables['latitude'][:]
lon1d = nc2.variables['longitude'][:]
print(lon1d)

years = np.arange(1999,2015)



f = open("/data/cloud/ucantseeme/hw03/cloud_info_daily/cloud_2002.txt","r", encoding='utf-8')
 
#lines = f.readlines() 

raw = f.read()


data =  StringIO(raw)
#aa = np.genfromtxt(data)
data_array = np.loadtxt(data)

#######################################
# calculating CWV





#print(aa.shape,aa[1:10,:])


#for line in f:
#	text.append(line)
#	print(line)
#	c = StringIO(line)
#	aa = np.genfromtxt(c,missing_values=999)
#	#print(np.loadtxt(c))
#	#print(aa[1])
	
############
# take out specific lat,lon cloud

print(data_array.shape)

cld_array = np.zeros((365,3,18))
ciwv_array = np.zeros((365,3,18))
#cld_array = np.zeros((75,10,36))
#cld_array = np.zeros((10,10,36))


#print(data_array[ np.logical_and( data_array[:,2] > 30, data_array[:,2] < 60) ,2])
#print(range(1,10))

#lat = np.arange(-)

#lat2d,lon2d = np.meshgrid(lat,lon,indexing='ij')

#plt.pcolor(lon2d.T,lat2d.T,.T)
for times in range(1,2921,8): #2921
	for lat in range(30,90,20):
		for lon in range(10,370,20):
			cld_lat_lon = data_array[ ( data_array[:,4] >= (lat-20)*4) & (data_array[:,4] <= lat*4)\
			& (data_array[:,3] >= (lon-10)*4) & (data_array[:,3] <= (lon+10)*4) & (data_array[:,0] == times) ,:]
            			
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
						#print(num_i,num_j)
						#print(num_i,num_j,'num')
						#time.sleep(0.003)


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

			cld_array[int((times-1)/8),int((lat-30)/20),int((lon-10)/20)] = d_geo_mean
			lat_index = np.where((lat_nc >= lat-60)&(lat_nc <= lat-40))[0]
			lon_index = np.where((lon_nc >= lon-10)&(lon_nc <= lon+10))[0]
			ciwv_array[int((times-1)/8),int((lat-30)/20),int((lon-10)/20)] = \
			np.mean((ciwv[int((times-1)/8),lat_index])[:,lon_index] < 40)				
			
		#print(int((times-1)/8),int((lat+40)/10),int((lon-10)/10))
print(d_geo_mean)


##########################
# reset ciwv to start at -180

ciwv_array_r = np.zeros((365,3,18))
ciwv_array_r[:,:,0:9] = ciwv_array[:,:,9:]
ciwv_array_r[:,:,9:] = ciwv_array[:,:,0:9]
ciwv_array = ciwv_array_r

#print(prec_per_size)
#plt.figure(0)
#plt.plot(np.arange(0,12)*20000,prec_per_size/prec_per_size[0],'.')

#plt.figure(1)
#plt.plot(400*np.power(2,np.arange(0,12)),prec_per_size_log/prec_per_size_log[0],'.')
#plt.xscale('log')

#plt.figure(2)
#plt.plot(400*np.power(2,np.arange(0,12)),num_per_size/num_per_size[0],'.')
#plt.xscale('log')

#plt.figure(3)

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

plt.figure(5)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)



m.pcolormesh(cx,cy,np.nanmean(ciwv_array,axis=0),cmap=plt.cm.jet)
plt.colorbar()

plt.figure(2)
#plt.plot(np.nanmean(np.nanmean(cld_array,axis=1),axis=1))


C = cld_array[:,0,9]
V = ciwv_array[:,0,9]

CD = C[~np.isnan(C)]
CV = V[~np.isnan(C)]

cor = np.corrcoef(CD,CV)

cor = np.zeros((3,18))
for lat_88 in range(0,3):
    for lon_88 in range(0,18):
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

