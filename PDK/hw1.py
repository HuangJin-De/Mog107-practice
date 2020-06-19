import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

m = Basemap(projection='cyl',lon_0=180,lat_0=0)

parallels = np.arange(-90,90,30)
meridians = np.arange(0,360,60)


m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)


years = np.arange(1980,2018)


path='/data/dadm1/reanalysis/ECMWF/ITM/daily'
nc = Dataset(path+'/Q/daily_interim_Q_1979.nc')

tt = np.arange(0,365)
zz = np.arange(0,27)

lon=nc.variables['longitude'][:]
lat=nc.variables['latitude'][:]
level=nc.variables['level'][zz]

print(tt)
q_mean_sum = np.mean(nc.variables['q'][tt,zz,:,:],axis=0)

print(nc.variables['q'][0:2,0:5,:,:].shape)

nc = Dataset(path+'/T/daily_interim_T_1979.nc')
T_mean_sum = np.mean(nc.variables['t'][tt,zz,:,:],axis=0)

nc = Dataset(path+'/Z/daily_interim_Z_1979.nc')
Z_mean_sum = np.mean(nc.variables['z'][tt,zz,:,:],axis=0)
for year in years:

	
	nc = Dataset(path+'/Q/daily_interim_Q_' + str(year) + '.nc')
	q_mean_sum = q_mean_sum + np.mean(nc.variables['q'][tt,zz,:,:],axis=0)
	nc = Dataset(path+'/T/daily_interim_T_' + str(year) + '.nc')
	T_mean_sum = T_mean_sum + np.mean(nc.variables['t'][tt,zz,:,:],axis=0)	
	nc = Dataset(path+'/Z/daily_interim_Z_' + str(year) + '.nc')
	Z_mean_sum = Z_mean_sum + np.mean(nc.variables['z'][tt,zz,:,:],axis=0)

print(level)

q_mean = q_mean_sum /(years.size+1)
T_mean = T_mean_sum /(years.size+1)
Z_mean = Z_mean_sum /(years.size+1)

lat2d,lon2d = np.meshgrid(lat,lon,indexing='ij')
level3d = np.ones((level.size,lat.size,lon.size))*level.reshape(level.size,1,1)
cwv = -np.trapz(2260000*q_mean+1005*T_mean+Z_mean,x=level3d*100,axis=0)/9.8/1005

print(lat2d.size,cwv.size)


plt.contourf(lon2d.T,lat2d.T,cwv.T)
plt.colorbar(orientation='horizontal')
plt.show()

