import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

#m = Basemap(projection='cyl',lon_0=180,lat_0=0,urcrnrlat=60)

#parallels = np.arange(-30,30,30)
#meridians = np.arange(0,360,60122
#m.drawcoastlines()
#m.drawparallels(parallels,labels=[1,1,0,0],fontsize=10)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)


years = np.arange(1999,2015)



# work out SCAI
# in one region
N
lat_series
lon_series
for num_i in range(0,num-1):
	for num_j in range(1,num-1):
		d_geo_mean = d_geo_mean * np.sqrt(np.square(lat_series[num_i]-lat_series[num_j])\
		+np.square(lon_series[num_i]-lon_series[num_j]))



	print(prec_per_size)
	plt.figure(0)
	plt.plot(np.arange(0,12)*20000,prec_per_size/prec_per_size[0],'.')
	
	plt.figure(1)
	plt.plot(400*np.power(2,np.arange(0,12)),prec_per_size_log/prec_per_size_log[0],'.')
	plt.xscale('log')
	
	plt.figure(2)
	plt.plot(400*np.power(2,np.arange(0,12)),num_per_size/num_per_size[0],'.')
	plt.xscale('log')
	plt.show()




#W_mean = W_mean_sum /(years.size+1)

#lat2d,lon2d = np.meshgrid(lat,lon,indexing='ij')
#level3d = np.ones((level.size,lat.size,lon.size))*level.reshape(level.size,1,1)
#cwv = -np.trapz(2260000*q_mean+1005*T_mean+Z_mean,x=level3d*100,axis=0)/9.8/1005

#print(lat2d.size,W_mean.size)


#plt.contourf(lon2d.T,lat2d.T,W_mean.T)
#plt.colorbar(orientation='horizontal')
#plt.show()

