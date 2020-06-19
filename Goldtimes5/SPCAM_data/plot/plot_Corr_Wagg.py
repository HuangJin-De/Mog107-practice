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

data_path= '/data/cloud/Goldtimes5/data/SPCAM/CPL64/'

#var_name= 'MSE_sur'
var_name= 'TS'
#var_name= 'RHmid'
#var_name= 'OMEGA500'
#var_name= 'CAPE'


#=== load data
CCWm_file= data_path+'CCW_gridmean.pk'
CCWv_file= data_path+'Wvar.pk'
CCWa_file= data_path+'Wagg.pk'
data_file= data_path+var_name+'.pk'
with open(CCWm_file, 'rb') as f:
    CCWm= pickle.load(f)

with open(CCWv_file, 'rb') as f:
    CCWv= pickle.load(f)

with open(CCWa_file, 'rb') as f:
    CCWa= pickle.load(f)

with open(data_file, 'rb') as f:
    data= pickle.load(f)

CRM_temp= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h1.0001-01-01-00000.nc'
with Dataset(CRM_temp,'r') as D:
    lat=D.variables['LAT_15s_to_30n'][:]
    lon=D.variables['LON_60e_to_180e'][:]
    maxGY= lat.size
    maxGX= lon.size

DOFm= np.zeros((maxGY,maxGX))
DOFv= np.zeros((maxGY,maxGX))
DOFa= np.zeros((maxGY,maxGX))
DOFd= np.zeros((maxGY,maxGX))
N= data.shape[1]
for yy in np.arange(0,maxGY):
  for xx in np.arange(0,maxGX):
    Cm_temp= np.squeeze(CCWm[0, :, yy,xx])
    Cv_temp= np.squeeze(CCWv[0, :, yy,xx])
#    Ca_temp= np.squeeze(Ca[0, :, yy,xx])
    data_temp= np.squeeze(data[0, :, yy,xx])
    Rm= np.corrcoef(Cm_temp[0:-2],Cm_temp[1:-1])[0,1]
    Rv= np.corrcoef(Cv_temp[0:-2],Cv_temp[1:-1])[0,1]
    Rd= np.corrcoef(data_temp[0:-2],data_temp[1:-1])[0,1]
    DOFm[yy,xx]= (1-Rm) / (1+Rm)
    DOFv[yy,xx]= (1-Rv) / (1+Rv)
    DOFd[yy,xx]= (1-Rd) / (1+Rd)

DOFa= DOFv # assume rho(agg)= rho(var)

# number index for DJF and JJA
day_num= np.zeros(365)
DJFm= day_num.copy()
JJAm= day_num.copy()
DJFm[0:59]=1
DJFm[334:]=1
JJAm[151:243]=1
DJF= (DJFm==1)
JJA= (JJAm==1)

for ss in np.arange(0,2):
    if ss==0:
       smask= DJF
       season= 'DJF'
    elif ss==1:
       smask= JJA    
       season= 'JJA'
    Cm= np.squeeze(CCWm[:, smask,:,:])
    Cv= np.squeeze(CCWv[:, smask,:,:])
    Ca= np.squeeze(CCWa[:, smask,:,:])
    datam= np.squeeze(data[:, smask,:,:])
    cCm= np.zeros((maxGY,maxGX))
    cCv= np.zeros((maxGY,maxGX))
    cCa= np.zeros((maxGY,maxGX))
    cCmCa= np.zeros((maxGY,maxGX))
    cCmCv= np.zeros((maxGY,maxGX))
    cCaCv= np.zeros((maxGY,maxGX))
    cCm_sig= np.zeros((maxGY,maxGX))
    cCv_sig= np.zeros((maxGY,maxGX))
    cCa_sig= np.zeros((maxGY,maxGX))
    cCmCa_sig= np.zeros((maxGY,maxGX))
    cCmCv_sig= np.zeros((maxGY,maxGX))
    cCaCv_sig= np.zeros((maxGY,maxGX))
    # 1 year data for now, need additional loop if more than 1 year
    for yy in np.arange(0,maxGY):
      for xx in np.arange(0,maxGX):
        Cm_temp= np.squeeze(Cm[:, yy,xx])
        Cv_temp= np.squeeze(Cv[:, yy,xx])
        Ca_temp= np.squeeze(Ca[:, yy,xx])
        datam_temp= np.squeeze(datam[:, yy,xx])
        Nm= Cm.shape[0]*DOFm[yy,xx] # DOF
        Nv= Cv.shape[0]*DOFv[yy,xx]
        Na= Ca.shape[0]*DOFa[yy,xx]
        Nd= datam.shape[0]*DOFd[yy,xx]
        de_nan= (~np.isnan(Ca_temp))
        Nmd= Ca_temp[de_nan].shape[0]*DOFm[yy,xx]
        Nvd= Ca_temp[de_nan].shape[0]*DOFv[yy,xx]
        Nad= Ca_temp[de_nan].shape[0]*DOFa[yy,xx]
        Ndd= Ca_temp[de_nan].shape[0]*DOFd[yy,xx]
        # cCm
        corr= np.corrcoef(Cm_temp,datam_temp)[0,1]
        N= Nm+Nd-2; t_value= corr* (N)**0.5 / (1- corr**2)**0.5
        t_crit= t_dist.ppf(0.95,N);
        cCm[yy,xx]= corr
        cCm_sig[yy,xx]= abs(t_value) - t_crit # >0: significant
        # cCv
        corr= np.corrcoef(Cv_temp,datam_temp)[0,1]
        N= Nv+Nd-2; t_value= corr* (N)**0.5 / (1- corr**2)**0.5
        t_crit= t_dist.ppf(0.95,N);
        cCv[yy,xx]= corr
        cCv_sig[yy,xx]= abs(t_value) - t_crit # >0: significant
        # cCa
        corr= np.corrcoef(Ca_temp[de_nan],datam_temp[de_nan])[0,1]
        N= Nad+Ndd-2; t_value= corr* (N)**0.5 / (1- corr**2)**0.5
        t_crit= t_dist.ppf(0.95,N);
        cCa[yy,xx]= corr
        cCa_sig[yy,xx]= abs(t_value) - t_crit # >0: significant
        # cCmCa
        corr= np.corrcoef(Ca_temp[de_nan],Cm_temp[de_nan])[0,1]
        N= Nad+Nmd-2; t_value= corr* (N)**0.5 / (1- corr**2)**0.5
        t_crit= t_dist.ppf(0.95,N);
        cCmCa[yy,xx]= corr
        cCmCa_sig[yy,xx]= abs(t_value) - t_crit # >0: significant
        # cCaCv
        corr= np.corrcoef(Ca_temp[de_nan],Cv_temp[de_nan])[0,1]
        N= Nad+Nvd-2; t_value= corr* (N)**0.5 / (1- corr**2)**0.5
        t_crit= t_dist.ppf(0.95,N);
        cCaCv[yy,xx]= corr
        cCaCv_sig[yy,xx]= abs(t_value) - t_crit # >0: significant
    cCm[(cCm_sig<0)]=np.nan
    cCv[(cCv_sig<0)]=np.nan
    cCa[(cCa_sig<0)]=np.nan
    cCmCa[(cCmCa_sig<0)]=np.nan
    cCaCv[(cCaCv_sig<0)]=np.nan

    print('Calculation Done')
    # plot
    plt.figure(figsize=(7,3))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    cs1= plt.contourf(lon,lat,cCv,levels=np.arange(-1,1.2,0.2),cmap='RdYlBu_r')
    plt.clim(-1,1)
#    plt.colorbar().set_ticks((-1,-0.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1))
    ax.set_xlim(60, 180)
    ax.set_ylim(-15, 30)
    picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Corr/'+var_name+'_Wvar_'+season+'.png'
    plt.savefig(picsave)
    plt.show()

    plt.figure(figsize=(7,3))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    cs1= plt.contourf(lon,lat,cCa,levels=np.arange(-1,1.2,0.2),cmap='RdYlBu_r')
    plt.clim(-1,1)
#    plt.colorbar().set_ticks((-1,-0.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1))
    ax.set_xlim(60, 180)
    ax.set_ylim(-15, 30)
    picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Corr/'+var_name+'_Wagg_'+season+'.png'
    plt.savefig(picsave)
    plt.show()

#    plt.figure(figsize=(7,3))
#    ax = plt.axes(projection=ccrs.PlateCarree())
#    ax.coastlines()
#    cs1= plt.contourf(lon,lat,cCaCv,levels=np.arange(-1,1.2,0.2),cmap='RdYlBu_r')
#    plt.clim(-1,1)
##    plt.colorbar().set_ticks((-1,-0.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1))
#    ax.set_xlim(60, 180)
#    ax.set_ylim(-15, 30)
#    picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Corr/Wvar_Wagg_'+season+'.png'
#    plt.savefig(picsave)
#    plt.show()

#plt.figure(figsize=(7,3))
#cs1= plt.contourf(lon,lat,cCa,levels=np.arange(-1,1.2,0.2),cmap='RdYlBu_r')
#plt.clim(-1,1)
#plt.colorbar().set_ticks((-1,-0.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1))
#picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Corr/Colorbar.png'
#plt.savefig(picsave)
#plt.show()
