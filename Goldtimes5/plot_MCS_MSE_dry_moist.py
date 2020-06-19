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
#plt.figure(figsize=(8,4))
#plt.plot(Tlat,np.zeros(Tlat.size),color=[.5, .5, .5],linewidth=1)
#plt.plot(Tlat,ERAm,color="k",linewidth=2)
#plt.xlim(0, 28)
#plt.ylim(320, 345)
#picsave= '/data/cloud/Goldtimes5/pic/TRMM_3b42/omega/Climatology.png'
#plt.savefig(picsave)
#plt.show()

# number index for DJF and JJA
#day_num= np.zeros(365).reshape((1,365))
#DJFm= day_num.copy()
#JJAm= day_num.copy()
#DJFm[0, 0:58]=1
#DJFm[0, 334:]=1
#JJAm[0,151:242]=1
#DJF= (DJFm==1)
#JJA= (JJAm==1)

#DJF_ERA= ERA[DJF,:,:]
#DJF_TR = TRMM[DJF,:,:]
#JJA_ERA= ERA[JJA,:,:]
#JJA_TR = TRMM[JJA,:,:]


cA= 300**2 *np.pi
cN= (TRMM <= 0)
cL= (TRMM >= cA)
cS= ~cN & ~cL

ratio= TRMM[cS].size / TRMM.size
print(ratio*100)

#
Ns= ERA.copy()
Ls= ERA.copy()
Ss= ERA.copy()

Ns[~cN] = np.nan
Ls[~cL] = np.nan
Ss[~cS] = np.nan
N_z= np.nanmean (Ns,axis=2); N_z= np.nanmean (N_z,axis=0)
S_z= np.nanmean (Ss,axis=2); S_z= np.nanmean (S_z,axis=0)
L_z= np.nanmean (Ls,axis=2); L_z= np.nanmean (L_z,axis=0)
   
plt.figure(figsize=(8,4))
plt.plot(Tlat,np.zeros(Tlat.size),color=[.5, .5, .5],linewidth=1)
plt.plot(Tlat,ERAm,color="k",linewidth=2)
plt.plot(Tlat,N_z,color="b",linewidth=2)
plt.plot(Tlat,S_z,color="m",linewidth=2)
plt.plot(Tlat,L_z,color="r", linewidth=2)
plt.xlim(0, 28)
plt.ylim(320, 345)
picsave= '/data/cloud/Goldtimes5/pic/TRMM_3b42/MSE/Small_large_sys.png'
plt.savefig(picsave)


