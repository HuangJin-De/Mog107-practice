#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:57:44 2020

@author: cloud
"""

# hw02
# built time: 2020/03/18

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# In[]
path = '/data/dadm1/obs/TRMM/'
filename_sze = 'TRMM3B42size/TRMMsize_3hrs_1998.nc' 
filename_pcp = '/TRMM3B42/3B42.1998.3hr.nc'

file      = nc.Dataset(path+filename_pcp)
lon = file.variables['longitude'][:]
lat = file.variables['latitude'][:]
pcp       = file.variables['pcp']

file      = nc.Dataset(path+filename_sze)
sze       = file.variables['objsize']

# In[]
sth = -10.125; nth = 30.125; wst = 100.125; est = 160.125;
sth_y = int(np.where(lat==sth)[0]); nth_y = int(np.where(lat==nth)[0]); 
wst_x = int(np.where(lon==wst)[0]); est_x = int(np.where(lon==est)[0]); 

# In[]
#fig = plt.figure(figsize=[9,3])
m = Basemap(projection='cyl', llcrnrlat=sth, urcrnrlat=nth,llcrnrlon=wst, urcrnrlon=est)
m.drawcoastlines()
parallels = np.arange(sth,nth,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
meridians = np.arange(wst,est,30.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ny = pcp.shape[1]; nx = pcp.shape[2]
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.
cs = m.contourf(x,y,np.flip(np.squeeze(sze[0,:,:]),axis=0), cmap='Reds')
cbar = m.colorbar(cs,location='right',pad="5%")

# In[] size distribution 
def plot_loghist(x, binsN, scale):
  # https://stackoverflow.com/questions/47850202/plotting-a-histogram-on-a-log-scale-with-matplotlib?rq=1
  lnrbins = np.linspace(10,1500,binsN+1)
  logbins = np.logspace(np.log10(lnrbins[0]),np.log10(lnrbins[-1]),len(lnrbins))
  if scale == 'linear':
      hist, bins, patches = plt.hist(x, bins=lnrbins)
      lower = 0 # don't give a shit for linear cases
      
  elif scale == 'log':
      hist, bins, patches = plt.hist(x, bins=logbins)
      lower = np.where(abs(bins-300) == min(abs(bins-300)))[0][0]
      
  hmax_i = np.where(hist == hist[lower:].max())[0][0]
  patches[hmax_i].set_fc('r')   
  plt.xscale(scale)
  plt.text(bins[hmax_i], hist[hmax_i], '{:.0f}'.format(0.5*(bins[hmax_i]+bins[hmax_i+1]))+' km')
  return hist, bins

sze2 = np.sqrt(np.reshape(sze[:, sth_y:nth_y, wst_x:est_x], -1))
# linear-scale and log scale for binning
scale = ['linear', 'log']
BinsN = [10, 50, 100, 200, 300]
Hist = [[]]*len(BinsN); Bins = [[]]*len(BinsN)
for j in range(len(scale)):
    for i in range(len(BinsN)):
        Hist[i], Bins[i] = plot_loghist(sze2, binsN=BinsN[i], scale=scale[j])
        plt.xlabel('size (km)')
        plt.ylabel('counts')
        plt.title('cloud size distribution ('+str(BinsN[i])+' bins)')
        plt.savefig('size_hist_'+scale[j]+'_b'+str(BinsN[i])+'.png', dpi=300)
        plt.clf()

# In[] size distribution (normalized to number)
def plot_loghist(x, binsN, scale):
  # https://stackoverflow.com/questions/47850202/plotting-a-histogram-on-a-log-scale-with-matplotlib?rq=1
  lnrbins = np.linspace(10,1500,binsN+1)
  logbins = np.logspace(np.log10(lnrbins[0]),np.log10(lnrbins[-1]),len(lnrbins))
  if scale == 'linear':
      hist, bins = np.histogram(x, bins=lnrbins)
      
      lower = 0 # don't give a shit for linear cases
      plt.plot(lnrbins[1:], hist/lnrbins[1:], 'o', markersize=1)
  elif scale == 'log':
      hist, bins = np.histogram(x, bins=logbins)
      
      lower = np.where(abs(bins-300) == min(abs(bins-300)))[0][0]
      plt.plot(logbins[1:], hist/logbins[1:], 'o', markersize=1)
  #hmax_i = np.where(hist == hist[lower:].max())[0][0]
  #patches[hmax_i].set_fc('r')   
  plt.xscale('log')
  plt.yscale('log')
  #plt.text(bins[hmax_i], hist[hmax_i], '{:.0f}'.format(0.5*(bins[hmax_i]+bins[hmax_i+1]))+' km')
  return hist, bins

sze2 = np.sqrt(np.reshape(sze[:, sth_y:nth_y, wst_x:est_x], -1))
# linear-scale and log scale for binning
scale = ['linear', 'log']
BinsN = [10, 50, 100, 200, 300]
Hist = [[]]*len(BinsN); Bins = [[]]*len(BinsN)
for j in range(len(scale)):
    for i in range(len(BinsN)):
        Hist[i], Bins[i] = plot_loghist(sze2, binsN=BinsN[i], scale=scale[j])
        plt.ylim(1E1,1E4)
        plt.xlabel('size (km)')
        plt.ylabel('counts')
        plt.title('cloud size distribution ('+str(BinsN[i])+' bins)')
        plt.savefig('size_hist_normed_'+scale[j]+'_b'+str(BinsN[i])+'.png', dpi=300)
        plt.clf()


# In[] seasonal cycle (2D)
import calendar
import seaborn as sb

binsN =[10, 50, 100, 200, 300]
for i in range(len(binsN)):
    SZE_mm = np.zeros((2,12,binsN[i]-1))
    Bins = [[]]*2
    Bins[0] = np.linspace(10,1500,binsN[i])
    Bins[1] = np.logspace(np.log10(Bins[0][0]),np.log10(Bins[0][-1]),binsN[i])
    for j in range(len(scale)):
        accum = 0
        for m in range(12):
            days = calendar.monthrange(1998,m+1)[1]
            hrs  = 8*days
            sze3 = np.sqrt(np.reshape(sze[accum:accum+hrs, sth_y:nth_y, wst_x:est_x], -1))
            SZE_mm[j,m,:], bins = np.histogram(sze3, bins=Bins[j]) 
            #SZE_mm[j,m,:] /= sum(SZE_mm[j,m,:])
            accum += hrs
        mm = [i+1 for i in range(12)]
        #sb.heatmap(SZE_mm[1,:,:])
        plt.contourf(Bins[j][1:], mm, np.squeeze(SZE_mm[j,:,:]))
        plt.grid()
        plt.xscale(scale[j])
        if i >=2:
            plt.xticks(ticks=Bins[j][::int(binsN[i]/10)], labels=['{:.0f}'.format(i) for i in Bins[j][::int(binsN[i]/10)]])
        plt.yticks(mm, labels=[str(i) for i in mm])
        plt.xlabel('size (km)')
        plt.ylabel('month')
        plt.colorbar()
        plt.title('seasonal cycle of cloud size distribution ('+str(binsN[i])+' bins)')
        plt.savefig('size_seasonal_'+scale[j]+'_b'+str(binsN[i])+'.png', dpi=300)
        plt.clf()
# 3D histogram?
