import csv
import numpy as np
import numpy.matlib
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt
from scipy.stats import t as t_dist
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy
import sys
sys.path.append('/data/cloud/Goldtimes5/code/')
from PYfunctions.statistics.smooth import running_mean
from PYfunctions.statistics.smooth import panta_mean as panta

data_path_map= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/'
data_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/Trop_mean/'
data_path_P= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/EOF_PRECT.pk'
filesaveT  = '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/Fast/Timeseries/' 
filesaveS  = '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/Fast/Corr/'

data_file= data_path+'RCE.pk'
with open(data_file, 'rb') as f:
    data= pickle.load(f)

data_file= data_path_map+'RCE.pk'
with open(data_file, 'rb') as f:
    data_map= pickle.load(f)

GCM_temp= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(GCM_temp,'r') as D:
    lat=D.variables['lat'][:]
    lon=D.variables['lon'][:]
    Y= lat.size
    X= lon.size

vars= ['PRECT']
data_filt={}
data_map_filt={}
self_corr= {}
self_corr_mask= {}
for vv in vars:
  #== trop mean data
  tmp= data[vv]
  tmp= tmp.reshape(tmp.size)
  tmp= panta(tmp)
  tmp_smooth= running_mean(tmp,6)
  tmp_f = tmp - tmp_smooth;  data_filt[vv]= tmp_f
  #== map data
  tmp= data_map[vv]
  tmp= tmp.reshape(tmp.shape[0]*tmp.shape[1], tmp.shape[2], tmp.shape[3])
  tmp= panta(tmp)
  tmp_smooth= running_mean(tmp,6)
  tmp_map_f= tmp - tmp_smooth; data_map_filt[vv]= tmp_map_f
  #== correlation
  Y,X= tmp_map_f.shape[1:3]
  cor= np.zeros((Y,X))
  cor_mask= np.zeros((Y,X))
  for yy in range(Y):
    for xx in range(X):
      tmp_map_f_loc= tmp_map_f[:, yy, xx]
      prec_mask= (tmp_map_f_loc > 10e-1)
      cor      [yy, xx]= np.corrcoef(tmp_map_f_loc[3:-4], tmp_f[3:-4])[0,1]
      cor_mask [yy, xx]= np.corrcoef(tmp_map_f_loc[prec_mask], tmp_f[prec_mask])[0,1]
  self_corr     [vv] = cor
  self_corr_mask[vv] = cor_mask

pr= self_corr['PRECT']
clevel= np.arange(-1,1.2,0.2)
plt.figure(figsize=(8,6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
cs1= plt.contourf(lon,lat, pr, levels=clevel ,cmap='RdYlBu_r')
#plt.clim(,cmax)
plt.colorbar(shrink=0.6, orientation='horizontal').set_ticks(clevel)
ax.set_xlim(-180, 180)
ax.set_ylim(-30, 30)
picsave= filesaveS+ 'PRECT_self.png'
plt.savefig(picsave)
plt.show()

pr= self_corr_mask['PRECT']
clevel= np.arange(-1,1.2,0.2)
plt.figure(figsize=(8,6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
cs1= plt.contourf(lon,lat, pr, levels=clevel ,cmap='RdYlBu_r')
#plt.clim(,cmax)
plt.colorbar(shrink=0.6, orientation='horizontal').set_ticks(clevel)
ax.set_xlim(-180, 180)
ax.set_ylim(-30, 30)
picsave= filesaveS+ 'PRECT_self_no_prec_masked.png'
plt.savefig(picsave)
plt.show()


#for ss in np.arange(0,2):
#    if ss==0:
#       smask= DJF
#       season= 'DJF'
#    elif ss==1:
#       smask= JJA    
#
#    plt.figure(figsize=(7,3))
#    ax = plt.axes(projection=ccrs.PlateCarree())
#    ax.coastlines()
#    cs1= plt.contourf(lon,lat,cCv,levels=np.arange(-1,1.2,0.2),cmap='RdYlBu_r')
#    plt.clim(-1,1)
#    plt.colorbar().set_ticks((-1,-0.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1))
#    ax.set_xlim(60, 180)
#    ax.set_ylim(-15, 30)
#    picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Corr/'+var_name+'_varCCW_'+season+'.png'
#    plt.savefig(picsave)
#    plt.show()

#    plt.figure(figsize=(7,3))
#    ax = plt.axes(projection=ccrs.PlateCarree())
#    ax.coastlines()
#    cs1= plt.contourf(lon,lat,cCa,levels=np.arange(-1,1,0.2),cmap='RdYlBu_r')
#    plt.clim(-1,1)
#    plt.colorbar().set_ticks((-1,-0.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1))
#    ax.set_xlim(60, 180)
#    ax.set_ylim(-15, 30)
#    picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Corr/'+var_name+'_aggCCW2_'+season+'.png'
#    plt.savefig(picsave)
#    plt.show()

#    plt.figure(figsize=(7,3))
#    ax = plt.axes(projection=ccrs.PlateCarree())
#    ax.coastlines()
#    cs1= plt.contourf(lon,lat,cCmCa,levels=np.arange(-1,1,0.2),cmap='RdYlBu_r')
#    plt.clim(-1,1)
##    plt.colorbar().set_ticks((-1,-0.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1))
#    ax.set_xlim(60, 180)
#    ax.set_ylim(-15, 30)
#    picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Corr/CCWmean_aggCCW2_'+season+'.png'
#    plt.savefig(picsave)
#    plt.show()

#    plt.figure(figsize=(7,3))
#    ax = plt.axes(projection=ccrs.PlateCarree())
#    ax.coastlines()
#    cs1= plt.contourf(lon,lat,cCaCv,levels=np.arange(-1,1.2,0.2),cmap='RdYlBu_r')
#    plt.clim(-1,1)
#    plt.colorbar().set_ticks((-1,-0.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1))
#    ax.set_xlim(60, 180)
#    ax.set_ylim(-15, 30)
#    picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Corr/CCWvar_aggCCW2_'+season+'.png'
#    plt.savefig(picsave)
#    plt.show()

#plt.figure(figsize=(7,3))
#cs1= plt.contourf(lon,lat,cCa,levels=np.arange(-1,1.2,0.2),cmap='RdYlBu_r')
#plt.clim(-1,1)
#plt.colorbar().set_ticks((-1,-0.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1))
#picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Corr/Colorbar.png'
#plt.savefig(picsave)
#plt.show()
