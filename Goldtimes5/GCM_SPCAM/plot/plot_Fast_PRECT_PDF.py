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

vars= ['PRECT']
data_filt={}
data_map_filt={}
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

pr    = data_filt['PRECT'][3:-4]
pr_map= data_map_filt['PRECT'][3:-4, :,:]
levels= np.arange()
Nt    = pr_map.shape[0]
for tt in range(Nt)
  tmp




plt.figure(figsize=(8,3))
plt.plot( np.arange(-5,6) ,np.zeros(11),'k',linewidth=0.5)
plt.plot( np.arange(-5,6) , lead_lag_RCE_Rad, 'b',linewidth=2.5, label= 'Rad. cooling')
plt.plot( np.arange(-5,6) , lead_lag_RCE_cv, 'r' ,linewidth=2.5, label= 'Conv. heating')
plt.plot( np.arange(-5,6) , lead_lag_RCE_pr, 'r--',linewidth=2.5, label= 'Prec. heating')
ax=plt.gca()
ax.set_xlim(-5, 5)
ax.set_ylim(-0.2, 1)
plt.legend(loc='lower right', fontsize=12)
picsave= filesaveT+ 'RCE_lead_lag_corr.png'
plt.savefig(picsave)
plt.show()


#  5 yr
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,6,5/365) ,np.zeros(RCEr.size),'k',linewidth=0.5)
plt.plot(np.arange(1,6,5/365) ,RCEr-np.nanmean(RCEr),'k',linewidth=1, label= 'RC imbalance')
plt.plot(np.arange(1,6,5/365) ,cvr-np.nanmean(cvr),'r',linewidth=1, label= 'Convection heating')
plt.plot(np.arange(1,6,5/365) ,prr-np.nanmean(prr),'r--',linewidth=.8, label= 'Precipitation')
plt.plot(np.arange(1,6,5/365) ,Radr-np.nanmean(Radr),'b',linewidth=1, label= 'Radiative cooling')
ax=plt.gca()
ax.set_xlim(1, 6)
ax.set_ylim(-20, 30)
plt.legend(loc='upper right', fontsize=8)
picsave= filesaveT+ 'RCE_decomp_5yr.png'
plt.savefig(picsave)
#plt.show()

# 1 yr
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,5*365+1,5) ,np.zeros(RCEr.size),'k',linewidth=0.5)
plt.plot(np.arange(1,5*365+1,5) ,RCEr-np.nanmean(RCEr),'k',linewidth=1, label= 'RC imbalance')
plt.plot(np.arange(1,5*365+1,5) ,cvr-np.nanmean(cvr),'r',linewidth=1, label= 'Convection heating')
plt.plot(np.arange(1,5*365+1,5) ,prr-np.nanmean(prr),'r--',linewidth=.8, label= 'Precipitation')
plt.plot(np.arange(1,5*365+1,5) ,Radr-np.nanmean(Radr),'b',linewidth=1, label= 'Radiative cooling')
ax=plt.gca()
ax.set_xlim(1, 366)
ax.set_ylim(-20, 30)
plt.legend(loc='upper right', fontsize=8)
picsave= filesaveT+ 'RCE_decomp_1yr.png'
plt.savefig(picsave)
#plt.show()

# Prec decompose
print(np.corrcoef(prr,CRMP)[0,1])
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,5*365+1,5) ,np.zeros(RCEr.size),'k',linewidth=0.5)
plt.plot(np.arange(1,5*365+1,5) ,prr-np.nanmean(prr),'r-',linewidth=1, label= 'Precipitation')
plt.plot(np.arange(1,5*365+1,5) ,CRMP-np.nanmean(CRMP),'k-',linewidth=1, label= 'Precipitation (CRM region)')
ax=plt.gca()
ax.set_xlim(1+365, 366+365)
ax.set_ylim(-20, 30)
plt.legend(loc='upper right', fontsize=8)
picsave= filesaveT+ 'PRECT_CRM_region.png'
plt.savefig(picsave)
#plt.show()

# Rad decompose
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,5*365+1,5) ,np.zeros(RCEr.size),'k',linewidth=0.5)
h1= plt.plot(np.arange(1,5*365+1,5) ,Radr-np.nanmean(Radr),'b',linewidth=2, label= 'Radiative cooling')
h2= plt.plot(np.arange(1,5*365+1,5) ,(SWtr-np.nanmean(SWtr))/5,'r',linewidth=1,label= 'SW (TOA)/5')
h3= plt.plot(np.arange(1,5*365+1,5) ,-(SWsr-np.nanmean(SWsr))/5,'r--',linewidth=1,label= 'SW (sur)/5')
h4= plt.plot(np.arange(1,5*365+1,5) ,-(LWtr-np.nanmean(LWtr)),'k',linewidth=1,label= 'LW (TOA)')
h5= plt.plot(np.arange(1,5*365+1,5) ,LWsr-np.nanmean(LWsr),'k--',linewidth=1,label= 'LW (sur)')
ax=plt.gca()
ax.set_xlim(1, 366)
ax.set_ylim(-10, 20)
plt.legend(loc='upper right', fontsize=8)
picsave= filesaveT+ 'Rad_decomp.png'
plt.savefig(picsave)
plt.show()

#plt.figure(figsize=(7,3))
#ax = plt.axes(projection=ccrs.PlateCarree())
#ax.coastlines()
#cs1= plt.contourf(lon,lat,cCv,levels=np.arange(-1,1.2,0.2),cmap='RdYlBu_r')


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
