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
from PYfunctions.statistics.smooth import smooth121 as sm
from PYfunctions.statistics.smooth import panta_mean as panta

data_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/Trop_mean/'
data_path_P= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/EOF_PRECT.pk'

data_file= data_path+'RCE.pk'
with open(data_file, 'rb') as f:
    data= pickle.load(f)

data_file= data_path_P
with open(data_file, 'rb') as f:
    pEOF= pickle.load(f)

data_file= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/CRM_region_mean/PRECT.pk'
with open(data_file, 'rb') as f:
    CRMP= pickle.load(f)

CRMP= CRMP.reshape((CRMP.size))
CRMP= CRMP*0.2583 # CRM region area / total tropical area

plat= pEOF['lat']; plon= pEOF['lon']

RCE= data['RCE']
#RCEm= np.mean(RCE, axis=0); 
#RCEm = RCEm.reshape((1,RCEm.size)); RCEm= np.tile((RCE.shape[0],1))
RCEr= RCE.reshape(RCE.size)
cv= data['PRECT']+ data['SHFLX']; cvr= cv.reshape(cv.size)
Rad= data['Rad_cool']; Radr= Rad.reshape(Rad.size)
pr= data['PRECT']; prr= pr.reshape(pr.size)
LWs= data['FLNS']; LWsr= LWs.reshape(LWs.size)
SWs= data['FSNS']; SWsr= SWs.reshape(SWs.size)
LWt= data['FLNT']; LWtr= LWt.reshape(LWt.size)
SWt= data['FSNTOA']; SWtr= SWt.reshape(SWt.size)
prl= data['PRECT_largeVar']; prlr= prl.reshape(prl.size)

RCEr= panta(RCEr)
cvr= panta(cvr)
prr= panta(prr)
prlr= panta(prlr)
Radr= panta(Radr)
CRMP= panta(CRMP)
LWsr= panta(LWsr)
SWsr= panta(SWsr)
LWtr= panta(LWtr)
SWtr= panta(SWtr)

#  5 yr
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,6,5/365) ,np.zeros(RCEr.size),'k',linewidth=0.5)
plt.plot(np.arange(1,6,5/365) ,RCEr-np.mean(RCEr),'k',linewidth=1, label= 'RC imbalance')
plt.plot(np.arange(1,6,5/365) ,cvr-np.mean(cvr),'r',linewidth=1, label= 'Convection heating')
plt.plot(np.arange(1,6,5/365) ,prr-np.mean(prr),'r--',linewidth=.8, label= 'Precipitation')
plt.plot(np.arange(1,6,5/365) ,Radr-np.mean(Radr),'b',linewidth=1, label= 'Radiative cooling')
ax=plt.gca()
ax.set_xlim(1, 6)
ax.set_ylim(-20, 30)
plt.legend(loc='upper right', fontsize=8)
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/TimeSeries/panta/RCE_decomp_5yr.png'
#plt.savefig(picsave)
#plt.show()

# 1 yr
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,5*365+1,5) ,np.zeros(RCEr.size),'k',linewidth=0.5)
plt.plot(np.arange(1,5*365+1,5) ,RCEr-np.mean(RCEr),'k',linewidth=1, label= 'RC imbalance')
plt.plot(np.arange(1,5*365+1,5) ,cvr-np.mean(cvr),'r',linewidth=1, label= 'Convection heating')
plt.plot(np.arange(1,5*365+1,5) ,prr-np.mean(prr),'r--',linewidth=.8, label= 'Precipitation')
plt.plot(np.arange(1,5*365+1,5) ,Radr-np.mean(Radr),'b',linewidth=1, label= 'Radiative cooling')
ax=plt.gca()
ax.set_xlim(1, 366)
ax.set_ylim(-20, 30)
plt.legend(loc='upper right', fontsize=8)
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/TimeSeries/panta/RCE_decomp_1yr.png'
#plt.savefig(picsave)
#plt.show()

# Prec decompose
print(np.corrcoef(prr[0:365],CRMP[0:365])[0,1])
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,5*365+1,5) ,np.zeros(RCEr.size),'k',linewidth=0.5)
plt.plot(np.arange(1,5*365+1,5) ,prr-np.mean(prr),'r-',linewidth=1, label= 'Precipitation')
plt.plot(np.arange(1,5*365+1,5) ,CRMP-np.mean(CRMP),'k-',linewidth=1, label= 'Precipitation (CRM region)')
ax=plt.gca()
ax.set_xlim(1+365, 366+365)
ax.set_ylim(-20, 30)
plt.legend(loc='upper right', fontsize=8)
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/TimeSeries/panta/PRECT_CRM_region.png'
plt.savefig(picsave)
plt.show()

# Prec decompose
print(np.corrcoef(prr,CRMP)[0,1])
plt.figure(figsize=(8,3))
#plt.plot(np.arange(1,5*365+1,1) ,np.zeros(RCEr.size),'k',linewidth=0.5)
plt.plot(np.arange(1,5*365+1,5) ,prr,'r-',linewidth=1, label= 'Precipitation')
plt.plot(np.arange(1,5*365+1,5) ,CRMP,'k-',linewidth=1, label= 'Precipitation (CRM region)')
ax=plt.gca()
ax.set_xlim(1, 366)
#ax.set_ylim(-20, 30)
plt.legend(loc='upper right', fontsize=8)
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/TimeSeries/panta/PRECT_CRM_region_ori.png'
plt.savefig(picsave)
plt.show()


# Rad decompose
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,5*365+1,5) ,np.zeros(RCEr.size),'k',linewidth=0.5)
h1= plt.plot(np.arange(1,5*365+1,5) ,Radr-np.mean(Radr),'b',linewidth=2, label= 'Radiative cooling')
h2= plt.plot(np.arange(1,5*365+1,5) ,(SWtr-np.mean(SWtr))/5,'r',linewidth=1,label= 'SW (TOA)/5')
h3= plt.plot(np.arange(1,5*365+1,5) ,-(SWsr-np.mean(SWsr))/5,'r--',linewidth=1,label= 'SW (sur)/5')
h4= plt.plot(np.arange(1,5*365+1,5) ,-(LWtr-np.mean(LWtr)),'k',linewidth=1,label= 'LW (TOA)')
h5= plt.plot(np.arange(1,5*365+1,5) ,LWsr-np.mean(LWsr),'k--',linewidth=1,label= 'LW (sur)')
ax=plt.gca()
ax.set_xlim(1, 366)
ax.set_ylim(-10, 20)
plt.legend(loc='upper right', fontsize=8)
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/RCE/TimeSeries/panta/Rad_decomp.png'
plt.savefig(picsave)
plt.show()




# number index for DJF and JJA
day_num= np.zeros(365)
DJFm= day_num.copy()
JJAm= day_num.copy()
DJFm[0:59]=1
DJFm[334:]=1
JJAm[151:243]=1
DJF= (DJFm==1)
JJA= (JJAm==1)



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
