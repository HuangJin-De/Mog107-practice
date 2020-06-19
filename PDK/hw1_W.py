import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

m = Basemap(projection='cyl',lon_0=180,lat_0=0,urcrnrlat=60)

parallels = np.arange(-30,30,30)
meridians = np.arange(0,360,60)


m.drawcoastlines()
m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)


years = np.arange(1980,2018)


path='/data/dadm1/reanalysis/ECMWF/ITM/daily'
nc = Dataset(path+'/W/daily_interim_W_1979.nc')

tt = np.arange(0,20)
zz = 15#np.arange(0,27)
yy = np.arange(80,161)

lon=nc.variables['longitude'][:]
lat=nc.variables['latitude'][yy]
level=nc.variables['level'][zz]


print(tt)
print(level)
print(lat)
W_mean_sum = np.mean(nc.variables['w'][tt,zz,yy,:],axis=0)

for year in years:

	
	nc = Dataset(path+'/W/daily_interim_W_' + str(year) + '.nc')
	W_mean_sum = W_mean_sum + np.mean(nc.variables['w'][tt,zz,yy,:],axis=0)



W_mean = W_mean_sum /(years.size+1)

lat2d,lon2d = np.meshgrid(lat,lon,indexing='ij')
level3d = np.ones((level.size,lat.size,lon.size))*level.reshape(level.size,1,1)
#cwv = -np.trapz(2260000*q_mean+1005*T_mean+Z_mean,x=level3d*100,axis=0)/9.8/1005

print(lat2d.size,W_mean.size)


plt.contourf(lon2d.T,lat2d.T,W_mean.T)
plt.colorbar(orientation='horizontal')
plt.show()

