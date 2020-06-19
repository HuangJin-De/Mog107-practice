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

data_file= data_path_P
with open(data_file, 'rb') as f:
    pEOF= pickle.load(f)
    plat= pEOF['lat']; plon= pEOF['lon']

data_file= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/CRM_region_mean/PRECT.pk'
with open(data_file, 'rb') as f:
    CRMP= pickle.load(f)

CRMP= CRMP.reshape((CRMP.size))
CRMP= panta(CRMP); CRMP_f= running_mean(CRMP,6)
CRMP= CRMP_f; #Slow
CRMP= CRMP*0.2583 # CRM region area / total tropical area

vars= ['RCE','PRECT','Rad_cool','SHFLX','FLNS','FSNS','FLNT','FSNTOA']
data_filt={}
data_map_filt={}
self_corr= {}
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
#  Y,X= tmp_map_f.shape[1:3]
#  cor= np.zeros(Y,X)
#  for yy in range(Y):
#    for xx in range(X):
#      tmp_map_f_loc= tmp_map_f[:, yy, xx]
#      cor [yy, xx]= np.corrcoef(tmp_map_f_loc, tmp_f)[0,1]
#  self_corr[vv] = cor

RCEr  = data_filt['RCE']
prr   = data_filt['PRECT']
cvr   = data_filt['SHFLX'] + prr
Radr  = data_filt['Rad_cool']
LWtr  = data_filt['FLNT']
LWsr  = data_filt['FLNS']
SWtr  = data_filt['FSNTOA']
SWsr  = data_filt['FSNS']

#=== lead lag correlation
lead_lag_RCE_Rad  = np.zeros(11)
lead_lag_RCE_cv   = np.zeros(11)
lead_lag_RCE_pr   = np.zeros(11)
Radr_s= Radr[3:-4]; RCEr_s= RCEr[3:-4]; cvr_s= cvr[3:-4]; prr_s= prr[3:-4]
lead_lag_RCE_Rad[5]= np.corrcoef(Radr_s, RCEr_s)[0,1]
lead_lag_RCE_cv[5] = np.corrcoef(cvr_s, RCEr_s)[0,1]
lead_lag_RCE_pr[5] = np.corrcoef(prr_s, RCEr_s)[0,1]
for i in range(1,6):
  # Rad cooling
  cor_lead= np.corrcoef(Radr_s[0:-1-i], RCEr_s[i:-1])[0,1]
  cor_lag= np.corrcoef(RCEr_s[0:-1-i], Radr_s[i:-1])[0,1]
  lead_lag_RCE_Rad[5-i]= cor_lead
  lead_lag_RCE_Rad[5+i]= cor_lag
  # convection heating
  cor_lead= np.corrcoef(cvr_s[0:-1-i], RCEr_s[i:-1])[0,1]
  cor_lag= np.corrcoef(RCEr_s[0:-1-i], cvr_s[i:-1])[0,1]
  lead_lag_RCE_cv[5-i]= cor_lead
  lead_lag_RCE_cv[5+i]= cor_lag
  # PRECT heating
  cor_lead= np.corrcoef(prr_s[0:-1-i], RCEr_s[i:-1])[0,1]
  cor_lag= np.corrcoef(RCEr_s[0:-1-i], prr_s[i:-1])[0,1]
  lead_lag_RCE_cv[5-i]= cor_lead
  lead_lag_RCE_cv[5+i]= cor_lag

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
ax.set_ylim(-10, 20)
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
ax.set_ylim(-10, 20)
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
