import numpy as np
import xarray as xr
import os

latrange = [-15, -5, 5, 15]
lonrange = [60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 181]

ds = xr.open_dataset('/data/dadm1/model_output/SPCAM/CTRL64/CTRL64.cam.h1.0001-04-01-00000.nc')
lat = ds['LAT_15s_to_30n'].values
nlat = np.size(lat)
lon = ds['LON_60e_to_180e'].values
nlon = np.size(lon)

gridx = [4,4,4,4,4,4,4,4,4,4,4,5]
gridy = [5,6,5]

# split to 10x10
cldnum = np.zeros([3650, 3, 12])
for t in range(10):
  for d in range(365):
    tmp = np.fromfile('/data/cloud/alfalfa/wk8/OBJ/SPCAM_CRM_{:0>2d}_{:0>3d}_cobj_daily.dat'.format(t+1, d+1), dtype = np.float32)
    size = os.path.getsize('/data/cloud/alfalfa/wk8/OBJ/SPCAM_CRM_{:0>2d}_{:0>3d}_cobj_daily.dat'.format(t+1, d+1))
    obj = tmp.reshape(size//4//8, 8)
    for y in range(3):
      for x in range(12):
        ind = np.where((obj[:,0] == d+1) &
		       (lon[obj[:,1].astype(int)] >= lonrange[x]) &
		       (lon[obj[:,1].astype(int)] <  lonrange[x+1]) &
		       (lat[obj[:,2].astype(int)] >= latrange[y]) &
		       (lat[obj[:,2].astype(int)] <  latrange[y+1]) &
		       (obj[:,5] <= 2.5E3) &
		       (obj[:,6] >= 1E4))[0]
        cldnum[t*365+d, y, x] = np.size(ind)#/(gridx[x]*gridy[y])	# per GCM grid
  print(t)
#cldnum.astype(np.float32).tofile('cvtcld_mean.dat')

'''# calculate every GCM grids
cldnum = np.zeros([3650, nlat, nlon])
for t in range(10):
  for d in range(365):
    tmp = np.fromfile('/data/cloud/alfalfa/wk8/OBJ/SPCAM_CRM_{:0>2d}_{:0>3d}_cobj_daily.dat'.format(t+1, d+1), dtype = np.float32)
    size = os.path.getsize('/data/cloud/alfalfa/wk8/OBJ/SPCAM_CRM_{:0>2d}_{:0>3d}_cobj_daily.dat'.format(t+1, d+1))
    obj = tmp.reshape(size//4//8, 8)
    for y in range(nlat):
      for x in range(nlon):
        ind = np.where((obj[:,1].astype(int) == x) &
                       (obj[:,2].astype(int) == y) &
                       (obj[:,5] <= 2.5E3) &
                       (obj[:,6] >= 1E4))[0]
        cldnum[t*365+d, y, x] = np.size(ind)
  print(t)
cldnum.astype(np.float32).tofile('cvtcld_all.dat')
'''
'''# no top and base constraint
cldnum = np.zeros([3650, nlat, nlon])
for t in range(10):
  for d in range(365):
    tmp = np.fromfile('/data/cloud/alfalfa/wk8/OBJ/SPCAM_CRM_{:0>2d}_{:0>3d}_cobj_daily.dat'.format(t+1, d+1), dtype = np.float32)
    size = os.path.getsize('/data/cloud/alfalfa/wk8/OBJ/SPCAM_CRM_{:0>2d}_{:0>3d}_cobj_daily.dat'.format(t+1, d+1))
    obj = tmp.reshape(size//4//8, 8)
    for y in range(nlat):
      for x in range(nlon):
        ind = np.where((obj[:,1].astype(int) == x) &
                       (obj[:,2].astype(int) == y))[0]
        cldnum[t*365+d, y, x] = np.size(ind)
  print(t)
cldnum.astype(np.float32).tofile('cld_all.dat')
'''

