#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 23:56:51 2020

@author: cloud
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# In[]
init_yr = 2001
finl_yr = 2015
yy = finl_yr-init_yr+1

path = '/data/cloud/ucantseeme/hw13/'
filename = 'rce_budget_{:4d}.nc'.format(init_yr)
file = nc.Dataset(path+filename)
net_flux_ave = file.variables['net_flux_ave'][:]

for y in range(1,yy):
    filename = 'rce_budget_{:4d}.nc'.format(init_yr+y)
    file = nc.Dataset(path+filename)
    net_flux_ave = np.concatenate((net_flux_ave, file.variables['net_flux_ave'][:]), axis=1)

# In[]
label = ['Pacific', 'meridional', 'zonal']
color = ['r', 'g', 'b']
fig,ax = plt.subplots()

for i in range(3):
    ax.hist(net_flux_ave[i,:], bins=20, density=True, histtype='step', label=label[i], color=color[i])
    
ax.grid()
ax.legend()
ax.set_xlim(-200,200)
ax.set_xlabel('W/m^2')
ax.set_ylabel('PDF')
ax.set_title('RCE imbalance (daily, 2001~2015)')
plt.savefig('RCEimbalance.png', dpi=300, transparent=False)