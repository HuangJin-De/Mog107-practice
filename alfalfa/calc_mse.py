import numpy as np
import xarray as xr
import math
import time
from mpi4py import MPI

print(time.asctime( time.localtime(time.time()) ))
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

dir1 = '/data/dadm1/reanalysis/ECMWF/ITM/daily/Q/'
dir2 = '/data/dadm1/reanalysis/ECMWF/ITM/daily/T/'
dir3 = '/data/dadm1/reanalysis/ECMWF/ITM/daily/Z/'
ds = xr.open_dataset(dir1 + 'daily_interim_Q_2000.nc')
lon = ds['longitude'].values
nlon = np.size(lon)
lat = ds['latitude'].values
nlat = np.size(lat)
lev = ds['level'].values
nlev = np.size(lev)
levtop = np.where(lev == 200)[0][0]
wlev = lev[:levtop+1]-lev[1:levtop+2]
wlev = wlev.reshape(1, levtop+1, 1, 1)

t_initial = 1998
nstep = 18
dt_bounds = nstep//size
tbounds = np.array([rank*dt_bounds, (rank+1)*dt_bounds])+t_initial

indt = np.arange(tbounds[0],tbounds[-1])
print(indt, rank)
lonrange = [45, 65, 85, 105, 125, 145, 165, 185, 205, 225]
latrange = [-30, -10, 10, 30]
#indt = np.arange(2000,2001)
crh_region = np.zeros((indt.size, 365, 3, 9))
Cp = 1004.      #[J K-1 kg-1][m2 s-2 K-1]
Lv = 2.5E6      #[m2 s-2]]
g = 9.8         #[m s-2]
rhow = 1000.
MSE_region = np.zeros((indt.size, 365, nlev, 3, 9))
for t in range(indt.size):
  ds1 = xr.open_dataset(dir1 + "daily_interim_Q_{:0>4d}.nc".format(indt[t]))
  ds2 = xr.open_dataset(dir2 + "daily_interim_T_{:0>4d}.nc".format(indt[t]))
  ds3 = xr.open_dataset(dir3 + "daily_interim_Z_{:0>4d}.nc".format(indt[t]))
  q_all = ds1['q'].values[:, :, :, :]
  t_all = ds2['t'].values[:, :, :, :]
  gz_all = ds3['z'].values[:, :, :, :]
  MSE = (Cp*t_all + gz_all + Lv*q_all)/Cp

  if ((indt[t]%4) != 0):
    for y in range(np.size(latrange)-1):
      for x in range(np.size(lonrange)-1):
        indy = np.where((lat >= latrange[y]) &
			(lat <= latrange[y+1]))[0]
        indx = np.where((lon >= lonrange[x]) &
			(lon <= lonrange[x+1]))[0]
        w = np.cos(lat[indy]*np.pi/180.)
        tmp = np.average(MSE[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1], axis = 2, weights = w)
        MSE_region[t, :, :, y, x] = np.mean(tmp, axis = 2)

  else:
    for y in range(np.size(latrange)-1):
      for x in range(np.size(lonrange)-1):
        indy = np.where((lat >= latrange[y]) &
                        (lat <= latrange[y+1]))[0]
        indx = np.where((lon >= lonrange[x]) &
                        (lon <= lonrange[x+1]))[0]
        w = np.cos(lat[indy]*np.pi/180.)
        tmp = np.average(MSE[:59, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1], axis = 2, weights = w)
        MSE_region[t, :59, :, y, x] = np.mean(tmp, axis = 2)
        tmp = np.average(MSE[60:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1], axis = 2, weights = w)
        MSE_region[t, 59:, :, y, x] = np.mean(tmp, axis = 2)

  print(t)
  print(time.asctime( time.localtime(time.time()) ))

if rank == 0:
  data = np.zeros((nstep, 365, nlev, 3, 9))
else:
  data = None

comm.Gather(MSE_region, data, root = 0)

if rank == 0:
  data.tofile('MSE_w20.dat')

print(time.asctime( time.localtime(time.time()) ))
