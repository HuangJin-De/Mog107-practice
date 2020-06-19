import csv
import numpy as np
import numpy.matlib
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy

data_path= '/data/cloud/Goldtimes5/data/SPCAM/CPL64/'
#=== load data
CCWm_file= data_path+'CCW_gridmean.pk'
CCWv_file= data_path+'CCW_gridvar.pk'
CCWa_file= data_path+'CCW_aggrgate.pk'
CAPE_file= data_path+'CAPE.pk'
RH_file= data_path+'RHmid.pk'
MSE_file= data_path+'MSE_sur.pk'

with open(CCWm_file, 'rb') as f:
    CCWm= pickle.load(f)

with open(CCWv_file, 'rb') as f:
    CCWv= pickle.load(f)

with open(CCWa_file, 'rb') as f:
    CCWa= pickle.load(f)

with open(CAPE_file, 'rb') as f:
    CAPE= pickle.load(f)

with open(RH_file, 'rb') as f:
    RH= pickle.load(f)

with open(MSE_file, 'rb') as f:
    MSE= pickle.load(f)

#CCWa[np.isnan(CCWa)]=50
#CaL= (CCWa==np.amax(CCWa))
#CaS= (CCWa==np.amin(CCWa))
#CAPE[~CaL]= np.mean(np.mean(np.mean(np.squeeze(CAPE),0),0),0)
#CAPE[~CaS]= np.mean(np.mean(np.mean(np.squeeze(CAPE),0),0),0)
#RH[~CaL]= np.mean(np.mean(np.mean(np.squeeze(RH),0),0),0)
#RH[~CaS]= np.mean(np.mean(np.mean(np.squeeze(RH),0),0),0)
#CAPE= np.squeeze(CAPE)
#RH= np.squeeze(RH)
#CAPESmax= np.argmax(CAPE)

CRM_path= '/data/dadm1/model_output/SPCAM/CPL64/'
CRM_temp= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h1.0001-01-01-00000.nc'
with Dataset(CRM_temp,'r') as D:
    lat=D.variables['LAT_15s_to_30n'][:]
    lon=D.variables['LON_60e_to_180e'][:]
    maxGY= lat.size
    maxGX= lon.size
    ilev = D.variables['ilev'][:]

dp= np.diff(ilev)
dp= dp[2:] # lowest 24 levels of GCM grid
dp2 = dp.reshape((dp.size,1))
dp2 = numpy.matlib.repmat(dp2,1,64)

# large CCWa
day1= 0; lat1= 10; lon1= 19;
# small CCWa
day2= 0; lat2= 1; lon2= 38;
print('S1: mean: '+str(CCWm[0,day1,lat1,lon1]))
print('S1: var: '+str(CCWv[0,day1,lat1,lon1]))
print('S1: agg: '+str(CCWa[0,day1,lat1,lon1]))
print('S1: CAPE: '+str(CAPE[0,day1,lat1,lon1]))
print('S1: RH: '+str(RH[0,day1,lat1,lon1]))
print('S1: MSE: '+str(MSE[0,day1,lat1,lon1]))

print('S2: mean: '+str(CCWm[0,day2,lat2,lon2]))
print('S2: var: '+str(CCWv[0,day2,lat2,lon2]))
print('S2: agg: '+str(CCWa[0,day2,lat2,lon2]))
print('S1: CAPE: '+str(CAPE[0,day2,lat2,lon2]))
print('S1: RH: '+str(RH[0,day2,lat2,lon2]))
print('S1: MSE: '+str(MSE[0,day2,lat2,lon2]))


agg_min= np.nanmin(CCWa)
agg_small= (CCWa== agg_min)
M= CCWm[agg_small]
V= CCWv[agg_small]
Mr= M.reshape((1,M.size))
Vr= V.reshape((1,V.size))
oMr= CCWm.reshape((1,CCWm.size))
oVr= CCWv.reshape((1,CCWv.size))

CRM_file1= CRM_path+'CPL64.cam.h1.0006-01-'+str(day1+1).zfill(2)+'-00000.nc'
CRM_file2= CRM_path+'CPL64.cam.h1.0006-01-'+str(day2+1).zfill(2)+'-00000.nc'
with Dataset(CRM_file1,'r') as D:
    qc= D.variables['CRM_QC_LON_60e_to_180e_LAT_15s_to_30n'][:,:,0,:,lat1,lon1]
    qc= np.squeeze(np.mean(qc,0))
    qc= qc * dp2 / np.sum(dp)
    hm1= np.sum(qc,0)
    ind1q= hm1.argsort()[-6:][::-1]
    W=  D.variables['CRM_W_LON_60e_to_180e_LAT_15s_to_30n'][:,16,0,:,lat1,lon2]
    W1= np.squeeze(np.mean(W,0))
    ind1W= W1.argsort()[-6:][::-1]

with Dataset(CRM_file2,'r') as D:
    qc= D.variables['CRM_QC_LON_60e_to_180e_LAT_15s_to_30n'][:,:,0,:,lat2,lon2]
    W=  D.variables['CRM_W_LON_60e_to_180e_LAT_15s_to_30n'][:,16,0,:,lat2,lon2]
    qc= np.squeeze(np.mean(qc,0))
    W2= np.squeeze(np.mean(W,0))
    qc= qc * dp2 / np.sum(dp)
    hm2= np.sum(qc,0)
    ind2q= hm2.argsort()[-6:][::-1]
    ind2W= W2.argsort()[-6:][::-1]

plt.figure(figsize=(7,3))
plt.plot(np.arange(1,65), hm1, color='r')
plt.plot(ind1q+1, hm1[ind1q], 'ro', linestyle='')
plt.plot(np.arange(1,65), hm2, color='b')
plt.plot(ind2q+1, hm2[ind2q], 'bo', linestyle='')
plt.xlim(1, 64)
picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Snapshot/Large_Small_agg.png'
plt.savefig(picsave)
plt.show()

plt.figure(figsize=(7,3))
plt.plot(np.arange(1,65), W1, color='r')
plt.plot(ind1W+1, W1[ind1W], 'ro', linestyle='')
plt.plot(np.arange(1,65), W2, color='b')
plt.plot(ind2W+1, W2[ind2W], 'bo', linestyle='')
plt.xlim(1, 64)
picsave= '/data/cloud/Goldtimes5/pic/SPCAM/Snapshot/Large_Small_agg_W.png'
plt.savefig(picsave)
plt.show()


