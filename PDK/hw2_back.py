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


path='/data/dadm1/obs/TRMM/TRMM3B42size'
nc = Dataset(path+'/TRMMsize_3hrs_1998.nc')

path2 = '/data/dadm1/obs/TRMM/TRMM3B42'
nc2 = Dataset(path2+'/3B42.1998.3hr.nc')

#print (nc2.variables)
cldsize = nc.variables['objsize'][0:2920,:,:]
hqp =  nc2.variables['pcp'][0:2920,:,:]

#MASK = np.ma.array(hqp, mask=hqp < 20)

prec_size_num_array = 0*cldsize[0,:,:]
print(prec_size_num_array.size)
AAA
prec_per_size = np.arange(0,12)
prec_per_size_log = np.arange(0,12)
num_per_size = np.arange(0,12)


#plt.contourf(cldsize)
#plt.colorbar(orientation='horizontal')
#print(np.sum(hqp[hqp < 10))


#tt = np.arange(0,20)
#zz = 15#np.arange(0,27)
#yy = np.arange(80,161)

#lon=nc.variables['longitude'][:]
#lat=nc.variables['latitude'][yy]



#print(tt)
#print(level)
#print(lat)
#W_mean_sum = np.mean(nc.variables['w'][tt,zz,yy,:],axis=0)

for year in years:
	nc = Dataset(path+'/TRMMsize_3hrs_' + str(year) +'.nc')
	nc2 = Dataset(path2+'/3B42.'+str(year)+'.3hr.nc')	
	cldsize = nc.variables['objsize'][0:2920,:,:]
	hqp =  nc2.variables['pcp'][0:2920,:,:]
	for jj in range(0,12):
		prec_per_size[jj] = prec_per_size[jj] + np.sum(hqp[cldsize > 20000*jj])	
		prec_per_size_log[jj] = prec_per_size_log[jj] + np.sum(hqp[cldsize > 400*np.power(2,jj)])	
		print(hqp[cldsize > 400*np.power(2,jj)].size,np.power(2,jj))
		print(hqp[cldsize > 400*np.power(2,jj)])
		num_per_size[jj] = num_per_size[jj] + hqp[cldsize > 400*np.power(2,jj)].size + 0.0
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

