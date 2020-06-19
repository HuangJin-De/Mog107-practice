#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 15:39:03 2020

@author: cloud
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# In[]
path = '/data/cloud/ucantseeme/hw03/'
filename = 'cld_obj_1998.nc' 

#path = '/data/dadm1/obs/TRMM/'
#filename = 'TRMM3B42size/TRMMsize_3hrs_1998.nc' 


file      = nc.Dataset(path+filename)
lon = file.variables['longitude'][:]
lat = file.variables['latitude'][:]
SZE       = file.variables['objsize'][:]

sth = lat[0]; nth = lat[-1]; wst = lon[0]; est = lon[-1];
sth_y = int(np.where(lat==sth)[0]); nth_y = int(np.where(lat==nth)[0]); 
wst_x = int(np.where(lon==wst)[0]); est_x = int(np.where(lon==est)[0]); 

# In[]
#fig = plt.figure(figsize=[9,3])
m = Basemap(projection='cyl', llcrnrlat=sth, urcrnrlat=nth,llcrnrlon=wst, urcrnrlon=est)
m.drawcoastlines()
parallels = np.arange(lat[0],lat[-1],10.)
m.drawparallels(parallels,fontsize=10)
meridians = np.arange(lon[0],lon[-1],30.)
m.drawmeridians(meridians,fontsize=10)
ny = SZE.shape[1]; nx = SZE.shape[2]
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.
cs = m.contourf(x,y,SZE[285,:,:])
cbar = m.colorbar(cs,location='right',pad="5%")
plt.savefig('test.png', dpi=300)