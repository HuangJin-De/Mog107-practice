#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 09:22:15 2020

@author: cloud
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

path = '/data/cloud/ucantseeme/hw08/'
filename = 'vertical_mass_flux.nc'

file = nc.Dataset(path+filename)
crh = file.variables['crh_scale'][:]
pre = file.variables['pre_scale'][:]
MF = file.variables['mass_flux'][:]
freq = file.variables['freq'][:]

# In[]
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots()
vmin = -4; vmax = 4
levels = np.arange(vmin, vmax, 1)
scale = 1E5
m = ax.pcolormesh(crh, pre/100, MF*scale, cmap='BrBG', vmin=vmin, vmax=vmax)#, levels=levels, extend='both')
ax.contour(crh, pre/100, MF*scale, levels=levels, colors='k')
ax.invert_yaxis()
ax.grid()
ax.set_xlim(0, 100)
fig.colorbar(m)
ax.set_ylabel('height (hPa)')
ax.set_title('vertical mass flux (*10^(-5))')

ax.xaxis.set_tick_params(labelbottom=False)
divider = make_axes_locatable(ax)
ax2 = divider.append_axes('bottom', 0.5, pad=0.1, sharex=ax)
ax2.set_xticks(np.arange(0,110,20))

ax2.plot(crh, freq/np.sum(freq))
ax2.set_xlabel('CRH (%)')
ax2.set_ylabel('%')
ax2.grid()
ax2.set_xlim(0, 100)
fig.savefig('vertical_mass_flux.png', dpi=300)

# In[]
filename = 'CRH_timeseries.nc'
file = nc.Dataset(path+filename)
crh = file.variables['crh_scale'][:]
crh_pdf = file.variables['crh_pdf'][:]*100

crh_pdf_clim = np.mean(crh_pdf, 0)
init_y = 1979; finl_y = 2015
yy = finl_y-init_y+1
mnth = np.arange(0,12*yy+1,1)
for d in range(3):
    fig, ax = plt.subplots(figsize=(10,10))
    m1 = d*10*12; m2 = m1+10*12
    y1 = init_y+d*10; y2 = y1+10
    mnth = np.linspace(y1, y2, 10*12)
    m = ax.pcolormesh(crh, mnth, crh_pdf[m1:m2,:], cmap='jet', vmin=0, vmax=2.25)
    ax.grid()
    ax.set_yticks(np.arange(y1, y2+1, 1))
    ax.set_xlabel('CRH (%)')
    fig.colorbar(m, shrink=0.5, label='%')
    fig.savefig('crh_'+str(y1)+'_'+str(y2)+'.png', dpi=300)
    
#m1 = (d+1)*10*12; m2 = m1+6*12
#y1 = init_y+d*10; y2 = y1+6
#mnth = np.linspace(y1, y2, 10*12)
#    m = ax.contourf(crh, mnth, crh_pdf[m1:m2,:], cmap='inferno')
#    ax.grid()
#    ax.set_yticks(np.arange(y1, y2+1, 1))
#    ax.set_xlabel('CRH (%)')
#    fig.colorbar(m, shrink=0.5, label='%')
#    fig.savefig('crh_'+str(y1)+'_'+str(y2)+'.png', dpi=300)