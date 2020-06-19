#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 22:30:28 2020

@author: cloud
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

path = '/data/cloud/ucantseeme/hw07/'
filename = 'CRH_OLR_freq.nc'

file = nc.Dataset(path+filename)
Lon = file.variables['lon'][:]
Lat = file.variables['lat'][:]
CRH = file.variables['crh'][:]
OLR = file.variables['olr'][:]
FRQ = file.variables['frequency'][:]
map_a = file.variables['map_a'][:]
map_b = file.variables['map_b'][:]
map_c = file.variables['map_c'][:]
map_d = file.variables['map_d'][:]
pdf_crh = file.variables['pdf_crh'][:]
pdf_olr = file.variables['pdf_olr'][:]
w_2pdf = file.variables['omega_2pdf'][:]
pdf_crh_olr = file.variables['pdf_crh_olr'][:]
pdf_w_olr = file.variables['pdf_w_olr'][:]

w_2pdf[w_2pdf == 0]=np.nan
# In[]
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib import ticker

YlGnBu = cm.get_cmap('YlGnBu', 256)
newcolors = YlGnBu(np.linspace(0, 1, 256))
grey = np.array([128/256, 128/256, 128/256, 1])
newcolors[:25, :] = grey
newcmp = ListedColormap(newcolors)

fig, ax = plt.subplots(figsize=(8,8))
divider = make_axes_locatable(ax)
ax_crh = divider.append_axes('bottom', 1, pad=0.1, sharex=ax)
ax_olr = divider.append_axes('right', 1, pad=0.1, sharey=ax)
ax_cb = divider.append_axes('right', 0.2, pad=0.1)
ax.xaxis.set_tick_params(labelbottom=False)
ax_olr.yaxis.set_tick_params(labelleft=False)

#m = ax.pcolormesh(CRH, OLR, FRQ*100, cmap=newcmp)#'YlGnBu')#, edgecolors='w')
vmin = 50; vmax = 400
levels = np.arange(vmin, vmax, 25)
m2 = ax.contourf(CRH, OLR, w_2pdf, cmap='Spectral_r', levels=levels, vmin=vmin, vmax=vmax, extend='both')
ax.contour(CRH, OLR, w_2pdf, levels=[250], colors='black')
levels = np.zeros(10)
vmin = -7; vmax = -1; n = 10
levels = np.linspace(vmin,vmax,n)
for i in range(n):
    levels[i] = 10**(levels[i])
m=ax.contour(CRH, OLR, FRQ*100, locator=ticker.LogLocator(), cmap='cool', levels=levels)#, vmin=levels[0], vmax=levels[-1])
ax.clabel(m)
ax.grid()
ax.set_ylabel('omega (Pa/s)')
ax.set_title('CRH, midlevel omega joint PDF (%) & OLR distribution')
#fig.colorbar(m, cax=ax_cb)
fig.colorbar(m2, cax=ax_cb, label='W/m^2')
ax_crh.plot(CRH, pdf_crh/np.sum(pdf_crh)*100)
ax_crh.grid()
ax_crh.set_xlim(CRH[0], CRH[99])
ax_crh.set_xlabel('CRH (%)')
ax_olr.plot(pdf_olr/np.sum(pdf_olr)*100, OLR)
ax_olr.grid()
ax_olr.set_ylim(OLR[0], OLR[174])
plt.savefig('CRH_omega_freq_olr.png', dpi=300)

# In[]
region = ['A', 'B', 'C', 'D']
MAP = [map_a, map_b, map_c, map_d]
for i in range(len(region)):
    fig, ax = plt.subplots(figsize=(6, 2))
    m = Basemap(projection='cyl',llcrnrlat=-30, llcrnrlon=0,urcrnrlat=30, urcrnrlon=360.5)
    m.drawcoastlines()
    parallels = np.arange(-30,31,10)
    m.drawparallels(parallels) # draw parallels 畫緯度線
    meridians = np.arange(0,361,60)
    m.drawmeridians(meridians) # draw meridians 畫經度線
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=11)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=11)
    LON,LAT = np.meshgrid(Lon,Lat)
    cx,cy = m(LON,LAT)
    CS = m.pcolormesh(cx,cy,MAP[i]*100,cmap = 'YlGnBu', shading='flat', vmin=0, vmax=0.1)
    m.colorbar(CS)
    plt.title(region[i]+' spatial frequency')
    plt.savefig(region[i]+'.png', dpi=300)
    plt.clf()
