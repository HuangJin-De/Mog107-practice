import csv
import numpy as np
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt

years= np.arange(2006,2016)
TRMM_data= '/data/cloud/Goldtimes5/data/TRMM3b42size.pk'
ERA5_data= '/data/cloud/Goldtimes5/data/ERA5_MSE_interped.pk'
# lon: 100-240, lat: 0-30

#=== load lat lon TRMM
TRMM_temp= '/data/dadm1/obs/TRMM/TRMM3B42/3B42.2004.3hr.nc'
with Dataset(TRMM_temp,'r') as D:
    data= D.variables['latitude'][:]
    Tlat_full= np.array(data)
    data= D.variables['longitude'][:]
    Tlon_full= np.array(data)

Tlat_mask= ( (Tlat_full > 0) & (Tlat_full <29.8) )
Tlon_mask= ( (Tlon_full > 100) & (Tlon_full <239.77) )
Tlat= Tlat_full[Tlat_mask]
Tlon= Tlon_full[Tlon_mask]
with open(TRMM_data, 'rb') as f:
    TRMM= pickle.load(f)

with open(ERA5_data, 'rb') as f:
    ERA= pickle.load(f)

# remove 2008/12/31 and 2012/12/31 from ERA, because it is not considered in TRMM size object. Incredible.
rm2008_i= 365*8*3
rm2012_i= 365*8*7
ERA= np.delete(ERA, np.arange(rm2008_i, rm2008_i+8), 0)
ERA= np.delete(ERA, np.arange(rm2012_i, rm2012_i+8), 0)

ERA=ERA/1000 # J/kg => kJ/kg
ERAm= np.mean (ERA,axis=2)
ERAm= np.mean (ERAm,axis=0)
plt.figure(figsize=(8,4))
plt.plot(Tlat,np.zeros(Tlat.size),color=[.5, .5, .5],linewidth=1)
plt.plot(Tlat,ERAm,color="k",linewidth=2)
plt.xlim(0, 28)
plt.ylim(320, 345)
picsave= '/data/cloud/Goldtimes5/pic/TRMM_3b42/omega/Climatology.png'
plt.savefig(picsave)
plt.show()

critAs= np.array([0, 100, 200, 250, 300, 350])
for critA in critAs:
    cA= critA**2 *np.pi
    cC= (TRMM > cA)
    ratio= TRMM[cC].size / TRMM.size
    print(critA, ratio*100)
    MCS= ERA.copy()
    nMCS= ERA.copy()
    MCS [TRMM <= cA]=np.nan
    nMCS [TRMM > cA]=np.nan
    MCS_z= np.nanmean (MCS,axis=2); MCS_z= np.nanmean (MCS_z,axis=0)
    nMCS_z= np.nanmean (nMCS,axis=2); nMCS_z= np.nanmean (nMCS_z,axis=0)   
    plt.figure(figsize=(8,4))
    plt.plot(Tlat,np.zeros(Tlat.size),color=[.5, .5, .5],linewidth=1)
    plt.plot(Tlat,ERAm,color="k",linewidth=2)
    plt.plot(Tlat,MCS_z,color="r",linewidth=2)
    plt.plot(Tlat,nMCS_z,color="b", linewidth=2)
    plt.xlim(0, 28)
    plt.ylim(320, 345)
    picsave= '/data/cloud/Goldtimes5/pic/TRMM_3b42/MSE/CritA_'+str(critA)+'.png'
    plt.savefig(picsave)

