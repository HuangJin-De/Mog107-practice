import csv
import numpy as np
import numpy.matlib
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
sys.path.append('/data/cloud/Goldtimes5/code/')
from PYfunctions.statistics.calc_areamean import areamean
from PYfunctions.plot.plot_func import add_colorbar

data_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/RCE.pk'
data_file= data_path
with open(data_file, 'rb') as f:
    data_all= pickle.load(f)

CRM_temp= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(CRM_temp,'r') as D:
    lat=D.variables['lat'][:]
    lon=D.variables['lon'][:]
    maxGY= lat.size
    maxGX= lon.size

# number index for DJF and JJA
day_num= np.zeros(365)
DJFm= day_num.copy()
JJAm= day_num.copy()
DJFm[0:59]=1
DJFm[334:]=1
JJAm[151:243]=1
DJF= (DJFm==1)
JJA= (JJAm==1)

#scmax= 450; scmin=50;
#scinter= (scmax-scmin)/10; scmax= scmax+scinter
slevel= np.arange(60, 550, 60)

#var_names= np.array(['RCE','Rad_cool','PRECT','FLNT','FLNS','FSNTOA','FSNS'])
var_names= np.array(['Rad_cool','Rad_cool_SW','Rad_cool_LW'])
#=== load data and plot
for vv in np.arange(0,var_names.shape[0]):
  var_name= var_names[vv]
  data= data_all[var_name] 
  data_DJF= data.copy(); data_DJF= data_DJF[:,DJF,:,:]; 
  szD= data_DJF.shape; data_DJF= data_DJF.reshape((szD[0]*szD[1], szD[2], szD[3]))
  data_JJA= data.copy(); data_JJA= data_JJA[:,JJA,:,:]
  szJ= data_JJA.shape; data_JJA= data_JJA.reshape((szJ[0]*szJ[1], szJ[2], szJ[3]))
  sz= data.shape; data_ann= data.reshape((sz[0]*sz[1], sz[2], sz[3]))
  data_DJFm= np.nanmean(data_DJF,0)
  data_JJAm= np.nanmean(data_JJA,0)
  data_annm= np.nanmean(data_ann,0)
  data_DJFs= np.std(data_DJF,axis=0, ddof=1)
  data_JJAs= np.std(data_JJA,axis=0, ddof=1)
  data_anns= np.std(data_ann,axis=0, ddof=1)
  for season in ['DJF', 'JJA', 'ann']:
    if season == 'DJF':
      data=data_DJFm
      datas=data_DJFs
    elif season == 'JJA':
      data=data_JJAm
      datas=data_JJAs
    elif season == 'ann':
      data=data_annm
      datas=data_anns
    if var_name== 'FSNS' :
      data= -data
    elif  var_name== 'FLNT':
      data= -data

    Cent= np.round(areamean(data, lat, lon, [0, 360, -30, 30]))
#    cinter= 40
    cinter= 10
    cmax= Cent+ cinter*6
    cmin= Cent- cinter*5
    level= np.arange(cmin,cmax,cinter)
    data[(data<cmin)]=cmin
    data[(data>cmax-cinter)]=cmax- cinter
#    datas[(datas<scmin)]=scmin
#    datas[(datas>scmax-scinter)]=scmax- scinter

    fig= plt.figure(figsize=(8,2.5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    cs1 = plt.contourf(lon,lat,data,levels=level,cmap='RdYlBu_r')
    plt.clim((cmin,cmax))
#    cs2= plt.contour(lon,lat, datas, linewidths=0.5, levels= slevel, colors='grey')
    cax = fig.add_axes([0.1, 0.1, 0.8, 0.03])
    plt.colorbar(orientation='horizontal', cax=cax).set_ticks(level)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-30, 30)
    picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/Spatial/mean/'+var_name+'_'+season+'.png'
    plt.savefig(picsave)
    plt.show()

#    fig= plt.figure(figsize=(8,2.5))
#    ax = plt.axes(projection=ccrs.PlateCarree())
#    ax.coastlines()
#    cs1= plt.contourf(lon,lat,datas,levels=slevel,cmap='RdYlBu_r')
#    plt.clim((scmin,scmax))
#    cax = fig.add_axes([0.1, 0.1, 0.8, 0.03])
#    plt.colorbar(orientation='horizontal', cax=cax).set_ticks(slevel)
#    ax.set_xlim(-180, 180)
#    ax.set_ylim(-30, 30)
#    picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/Spatial/std/'+var_name+'_'+season+'.png'
#    plt.savefig(picsave)
#    plt.show()
#
