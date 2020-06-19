import csv
import numpy as np
import numpy.matlib
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy
import sys
sys.path.append('/data/cloud/Goldtimes5/code/')
from PYfunctions.statistics.calc_areamean import areamean
from PYfunctions.statistics.EOFfun import EOF
from PYfunctions.statistics.smooth import panta_mean as panta

data_file= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/PRECT.pk'
datasave_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/'

# map
GCM_temp= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(GCM_temp,'r') as D:
  lat=D.variables['lat'][:] # level in hPa
  lon=D.variables['lon'][:] # level in hPa

maxY= lat.size
maxX= lon.size

Lv= 2.477*10**6 # J/kg
rho= 1000 # kg/m3

#=== load data and plot
with open(data_file, 'rb') as f:
  data= pickle.load(f)

#=== data process
Nyr, Nd=  data.shape[0:2]
data= data * Lv * rho # convert m/s to W/m2
data= data.reshape(( Nyr*Nd, maxY, maxX)) # Time, lat, lon
data= panta(data); Nt= data.shape[0] # panta mean
data_mean= np.mean(data, axis=0).reshape((1, maxY, maxX)); 
data_mean= np.tile(data_mean, [Nt, 1, 1])
data= data-data_mean # anomaly

#=== weighting matrix
weight= np.array( np.cos(lat* np.pi/180)**0.5, dtype= 'float' )
weight= weight.reshape(( 1, maxY, 1))
weight= np.tile(weight, [Nt, 1, maxX])
data= data* weight # weight by square cos(lat)

#=== mask
lat_mask= np.squeeze( np.array([(lat>= -30) & (lat<= 30)], dtype= 'bool'))
data_r= data[:, lat_mask, :]
Nt2, Ny, Nx = data_r.shape

#=== reshape
data_rs= data_r.reshape((Nt2, Ny*Nx))
EOFp= EOF(data_rs)
EOFp['lat']= lat[lat_mask]
EOFp['lon']= lon

filesave= datasave_path+ 'EOF_PRECT_panta.pk'
with open(filesave, 'wb') as fs:
  pickle.dump(EOFp, fs, protocol = 4)
