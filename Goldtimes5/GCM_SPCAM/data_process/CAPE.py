import csv
import numpy as np
import numpy.matlib
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy
from scipy.interpolate import interp1d

data_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/'

var_name= 'CAPE'
# var_name= 'TS'
# var_name= 'OMEGA500'

#=== load data
Q_file= data_path+'Q.pk'
T_file= data_path+'T.pk'
Z_file= data_path+'Z3.pk'
PS_file= data_path+'PS.pk'
MSE_file= data_path+'MSE.pk'

with open(Q_file, 'rb') as f:
    Q_all= pickle.load(f)
    print(Q_all.shape)
    Q_all= Q_all[:,:,-1,:,:] # surface layer
    print(Q_all.shape)

with open(T_file, 'rb') as f:
    Te_all= pickle.load(f)
    sz= Te_all.shape

with open(MSE_file, 'rb') as f:
    MSE= pickle.load(f)
    MSEsur = MSE [:, :, -1, :, :]
    MSEsur = MSEsur.reshape((sz[0], sz[1], 1, sz[3], sz[4]))

with open(Z_file, 'rb') as f:
    Z3_all= pickle.load(f)

with open(PS_file, 'rb') as f:
    PS_all= pickle.load(f)

P= PS_all.reshape((sz[0], sz[1], 1, sz[3], sz[4]))

# parameters
Cp= 240 * 4.186 # J/kg/K
Lv= 2.477*10**6 # J/kg
g= 9.8

#=== load lev
GCM_temp= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(GCM_temp,'r') as D:
    lev=D.variables['lev'][:] # level in hPa
    ilev=D.variables['ilev'][:] # level in hPa

lev=lev / 1000 # ratio
ilev=ilev / 1000 # ratio
dp= np.diff(ilev)
maxZ= lev.size

#=== pressure field
P= np.tile(P, [1,1,lev.size,1,1]) # PS for now
dp3= dp.reshape((1,1,dp.size,1,1))
dp3= np.tile(dp3,[sz[0], sz[1], 1, sz[3], sz[4]])
dp3= P * dp3 # Pa

lev3= lev.reshape((1,1,lev.size,1,1))
lev3= np.tile(lev3,[sz[0], sz[1], 1, sz[3], sz[4]])
P=P * lev3 # Pa, PS to Pressure

#=== MSE*
es= (1.0007+(3.46e-6*(P/100)))*6.1121*np.exp(17.502*(Te_all-273.15)/(240.97+(Te_all-273.15))) # sat. vapor pres
ws  = .622*es/((P/100)-es)
qs= ws/(1+ws) # mixing ratio (kg/kg)

MSEsat= Cp*Te_all + Lv*qs + g*Z3_all
MSEsur= np.tile(MSEsur,[1,1,maxZ,1,1])

#=== difference
Diff= MSEsur - MSEsat
pos= (Diff > 0)
weight= np.zeros(sz)
weight[pos]=1
dpw= dp3*weight; 
Diffw= Diff*dpw
CAPE= np.sum(Diffw,2) / np.sum(dpw,2)

filesave= data_path+'CAPE.pk' 
with open(filesave, 'wb') as fs:
  pickle.dump(CAPE, fs, protocol = 4)






