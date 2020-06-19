import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import glob
import os
import sys
import time
from mpi4py import MPI

print(time.asctime( time.localtime(time.time()) ))
##### MPI settings
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

t_initial = 0
nstep = 3650
dt_bounds = nstep//size
tbounds = np.array([rank*dt_bounds, (rank+1)*dt_bounds])+t_initial

indt = np.arange(tbounds[0],tbounds[-1])
##### MPI settings #####
sys.setrecursionlimit(1000000)
##### define functions
## connect 4 directions
def stamp(data, mask, ny, nx, criteria, j, i, cldnum):
    if (j>0):
        if ((data[j-1,i]>=criteria) and (mask[j-1,i]==0)):
            mask[j-1,i] = cldnum
            mask = stamp(data, mask, ny, nx, criteria, j-1, i, cldnum)
    if (i>0):
        if ((data[j,i-1]>=criteria) and (mask[j,i-1]==0)):
            mask[j,i-1] = cldnum
            mask = stamp(data, mask, ny, nx, criteria, j, i-1, cldnum)
    if (j<ny-1):
        if ((data[j+1,i]>=criteria) and (mask[j+1,i]==0)):
            mask[j+1,i] = cldnum
            mask = stamp(data, mask, ny, nx, criteria, j+1, i, cldnum)
    if (i<nx-1):
        if ((data[j,i+1]>=criteria) and (mask[j,i+1]==0)):
            mask[j,i+1] = cldnum
            mask = stamp(data, mask, ny, nx, criteria, j, i+1, cldnum)
    return mask

## find object
def find_cloud(data, mask, ny, nx):
    criteria = 1E-5
    cldnum = 0
    for j in range(ny):
        for i in range(nx):
            if (data[j,i] < criteria):
                pass
            elif (mask[j,i] != 0):
                pass
            else:
                cldnum = cldnum+1
                mask[j,i] = cldnum
                mask = stamp(data, mask, ny, nx, criteria, j, i, cldnum)
    
    return mask, cldnum
##### define functions #####

path = '/data/dadm1/model_output/SPCAM/CTRL64/'
filename = glob.glob(path+'CTRL64.cam.h1.00*-*-*-00000.nc')
filename = sorted(filename)
filename2= glob.glob(path+'CTRL64.cam.h0.00*-*-*-00000.nc')
filename2= sorted(filename2)

#=== YJC load pressure
ds  = xr.open_dataset(filename2[100])
hyai= np.array(ds['hyai']); hybi= np.array(ds['hybi']);
PS  = np.array(ds['PS'])
PS  = np.mean(PS, axis=0)
PS  = PS.reshape((1, PS.shape[0], PS.shape[1])); PS= np.tile(PS, [hyai.size, 1, 1])
hyai= hyai.reshape((hyai.size, 1, 1)); hyai= np.tile(hyai,[1, PS.shape[1], PS.shape[2]])
hybi= hybi.reshape((hybi.size, 1, 1)); hybi= np.tile(hybi,[1, PS.shape[1], PS.shape[2]])
P   = hyai * 100000 + hybi * PS
dp  = P[2:-1,:,:] - P[1:-2,:,:]
#=== YJC end

ds = xr.open_dataset(filename[100])
lat = ds['LAT_15s_to_30n'].values
nlat = np.size(lat)
lon = ds['LON_60e_to_180e'].values
nlon = np.size(lon)

ds = xr.open_dataset(filename2[100])
latall = ds['lat'].values
lonall = ds['lon'].values
latindex = np.where((latall <= 30)&(latall >= -15))[0]
nlatindex = np.size(latindex)
lonindex = np.where((lonall <= 180)&(lonall >= 60))[0]
nlonindex = np.size(lonindex)

#objlist = np.empty([1,8])
for t in range(indt.size):
  ds = xr.open_dataset(filename[indt[t]])
  cld = ds['CRM_QC_LON_60e_to_180e_LAT_15s_to_30n'].values# [time, crm_nz, crm_ny, crm_nx, LAT_15s_to_30n, LON_60e_to_180e]
  tot = ds['CRM_QI_LON_60e_to_180e_LAT_15s_to_30n'].values
  cld_daily = np.nanmean(cld, axis = 0)
  cld_daily = np.where((cld_daily < 1E-5), 0.0, cld_daily)
  tot_daily = np.nanmean(tot, axis = 0)

  ds = xr.open_dataset(filename2[indt[t]])
  z3 = ds['Z3'].values# [time, lev, lat, lon]
  z_all = np.nanmean(z3, axis = 0)

  objlist = np.empty([1,10])
  for y in range(nlat):
    for x in range(nlon):
      data = cld_daily[:,0,:,y,x]
      totdata = tot_daily[:,0,:,y,x]
      shape = np.shape(data)
      mask = np.zeros(shape)
      mask_result, cldnum = find_cloud(data, mask, shape[0], shape[1])

      z = z_all[:, latindex[y], lonindex[x]]
      z = z[::-1]
      z = z[:shape[0]]
      z = z-z[0]
      #=== YJC dp
      dp_tmp = dp[:, latindex[y], lonindex[x]]
      dp_tmp= [::-1]; dp_tmp= [:shape[0]]
      
      dp_tmp = dp.reshape(())

      summary = np.empty([cldnum, 10])
      for n in range(cldnum):
        summary[n, 0] = (t%365)+1
        summary[n, 1] = x
        summary[n, 2] = y
        summary[n, 3] = n+1

        tmp = np.zeros(shape)
        #tmp[mask_result==n+1] = 1
        tmp = np.where((mask_result==n+1), 1, 0)

        summary[n, 4] = np.sum(tmp)
        for k in range(shape[0]):
          if any(tmp[k,:]):
            summary[n, 5] = z[k]
            break
        for k in reversed(range(shape[0])):
          if any(tmp[k,:]):
            summary[n, 6] = z[k]
            break
        summary[n, 7] = summary[n, 6] - summary[n, 5]

        tmp = np.where((mask_result==n+1), data, 0)
        summary[n, 8] = np.sum(tmp)
        tmp = np.where((mask_result==n+1), totdata, 0)
        summary[n, 9] = np.sum(tmp)
        
      objlist = np.append(objlist, summary, axis = 0)
  objlist = np.delete(objlist, 0, 0)

  objlist.astype(np.float32).tofile('/data/cloud/alfalfa/wk13/OBJ/SPCAM_CRM_{:0>2d}_{:0>3d}_cobj_daily.dat'.format((rank+1), ((t%365)+1)))
  print(rank, indt[t])
#objlist = np.delete(objlist, 0, 0)  

#objlist.astype(np.float32).tofile('SPCAM_CRM_{:0>2d}_cobj_daily.dat'.format((indt[t]//365)+1))
print(time.asctime( time.localtime(time.time()) ), ' rank ', rank, ' Done')
