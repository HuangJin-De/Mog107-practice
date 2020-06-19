import Ngl
import os
import xarray as xr
import numpy as np
#import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import glob
import time
from mpi4py import MPI
#import matplotlib.colors as colors
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#from matplotlib.axes import Axes
#from cartopy.mpl.geoaxes import GeoAxes
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib.ticker as mticker

print(time.asctime( time.localtime(time.time()) ))
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

t_initial = 0
nstep = 3650
dt_bounds = nstep//size
tbounds = np.array([rank*dt_bounds, (rank+1)*dt_bounds])+t_initial

indt = np.arange(tbounds[0],tbounds[-1])

path = '/data/dadm1/model_output/SPCAM/CTRL64/'
filename = glob.glob(path+'CTRL64.cam.h0.00*-*-*-00000.nc')
filename = sorted(filename)
# create new pressure level
pint = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775,
         750, 700, 650, 600, 550, 500, 450, 400, 350, 300]#,
         #250, 225, 200, 175, 150, 125, 100]
npint = np.size(pint)
pnew = []
for (item1, item2) in zip(pint[:npint-1], pint[1:npint]):
  pnew.append((item1+item2)/2)
dp = []
for (item1, item2) in zip(pint[:npint-1], pint[1:npint]):
  dp.append(item1-item2)
# focus on the region
ds = xr.open_dataset(filename[0])
lat = ds['lat']
nlat = np.size(lat)
lon = ds['lon']
nlon = np.size(lon)

latindex = np.where((lat <= 15)&(lat >= -15))[0]
nlatindex = np.size(latindex)
lonindex = np.where((lon <= 180)&(lon >= 60))[0]
nlonindex = np.size(lonindex)
lonrange = [60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 181]
latrange = [-15, -5, 5, 15]
#print(lon[lonindex])
#print(lat[latindex])

Cp = 1005.
Lv = 2.5E6
g = 9.81
#calculate crh in 60E-180E & 15S-15N (start)
#Td = np.zeros((indt.size, npint, np.size(latrange)-1, np.size(lonrange)-1))
#qd = np.zeros((indt.size, npint, np.size(latrange)-1, np.size(lonrange)-1))
#ud = np.zeros((indt.size, npint, np.size(latrange)-1, np.size(lonrange)-1))
#vd = np.zeros((indt.size, npint, np.size(latrange)-1, np.size(lonrange)-1))
omegad = np.zeros((indt.size, npint, np.size(latrange)-1, np.size(lonrange)-1))
for i in range(indt.size):
  ds = xr.open_dataset(filename[indt[i]])
  hyam = ds["hyam"]
  hybm = ds["hybm"]
#  T    = ds['T'][:, :, latindex, lonindex]
#  q    = ds['Q'][:, :, latindex, lonindex]
#  u    = ds['U'][:, :, latindex, lonindex]
#  v    = ds['V'][:, :, latindex, lonindex]
  omega = ds['OMEGA'][:, :, latindex, lonindex]
  psrf = ds["PS"][:, latindex, lonindex]
  P0mb =  0.01*ds["P0"]
  #lat = ds["lat"]
  #lon = ds["lon"]

  intyp = 1
  kxtrp = False

#  Tnew = Ngl.vinth2p(T,hyam,hybm,pint,psrf,intyp,P0mb,1,kxtrp)
#  Tnew[Tnew==1e30] = np.NaN
#  qnew = Ngl.vinth2p(q,hyam,hybm,pint,psrf,intyp,P0mb,1,kxtrp)
#  qnew[qnew==1e30] = np.NaN
#  unew = Ngl.vinth2p(u,hyam,hybm,pint,psrf,intyp,P0mb,1,kxtrp)
#  unew[unew==1e30] = np.NaN
#  vnew = Ngl.vinth2p(v,hyam,hybm,pint,psrf,intyp,P0mb,1,kxtrp)
#  vnew[vnew==1e30] = np.NaN
  omeganew = Ngl.vinth2p(omega,hyam,hybm,pint,psrf,intyp,P0mb,1,kxtrp)
  omeganew[omeganew==1e30] = np.NaN
  for y in range(np.size(latrange)-1):
    for x in range(np.size(lonrange)-1):
      indy = np.where((lat[latindex] >= latrange[y]) &
                      (lat[latindex] <  latrange[y+1]))[0]
      indx = np.where((lon[lonindex] >= lonrange[x]) &
                      (lon[lonindex] <  lonrange[x+1]))[0]
      w = np.cos(lat[latindex][indy]*np.pi/180.)

#      masked_T = np.ma.masked_array(Tnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1], np.isnan(Tnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1]))
#      masked_q = np.ma.masked_array(qnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1], np.isnan(qnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1]))
#      masked_u = np.ma.masked_array(unew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1], np.isnan(unew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1]))
#      masked_v = np.ma.masked_array(vnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1], np.isnan(vnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1]))
      masked_omega = np.ma.masked_array(omeganew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1], np.isnan(omeganew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1]))
#      T_w = np.ma.average(masked_T, axis = 2, weights = w)
#      q_w = np.ma.average(masked_q, axis = 2, weights = w)
#      u_w = np.ma.average(masked_u, axis = 2, weights = w)
#      v_w = np.ma.average(masked_v, axis = 2, weights = w)
      omega_w = np.ma.average(masked_omega, axis = 2, weights = w)

#      tmp = np.nanmean(T_w, axis = 2)
#      Td[i, :, y, x] = np.nanmean(tmp, axis = 0)
#      tmp = np.nanmean(q_w, axis = 2)
#      qd[i, :, y, x] = np.nanmean(tmp, axis = 0)
#      tmp = np.nanmean(u_w, axis = 2)
#      ud[i, :, y, x] = np.nanmean(tmp, axis = 0)
#      tmp = np.nanmean(v_w, axis = 2)
#      vd[i, :, y, x] = np.nanmean(tmp, axis = 0)
      tmp = np.nanmean(omega_w, axis = 2)
      omegad[i, :, y, x] = np.nanmean(tmp, axis = 0)

  print(rank, i)

R = 287.05
#rho = np.reshape(pint, [1, npint, 1, 1])/(R*Td)
print(time.asctime( time.localtime(time.time()) ))

if rank == 0:
#  data_rho = np.zeros((nstep, npint, np.size(latrange)-1, np.size(lonrange)-1))
#  data_qd = np.zeros((nstep, npint, np.size(latrange)-1, np.size(lonrange)-1))
#  data_ud = np.zeros((nstep, npint, np.size(latrange)-1, np.size(lonrange)-1))
#  data_vd = np.zeros((nstep, npint, np.size(latrange)-1, np.size(lonrange)-1))
  data_omegad = np.zeros((nstep, npint, np.size(latrange)-1, np.size(lonrange)-1))
else:
#  data_rho = None
#  data_qd = None
#  data_ud = None
#  data_vd = None
  data_omegad = None
#comm.Gather(rho, data_rho, root = 0)
#comm.Gather(qd, data_qd, root = 0)
#comm.Gather(ud, data_ud, root = 0)
#comm.Gather(vd, data_vd, root = 0)
comm.Gather(omegad, data_omegad, root = 0)
if rank == 0:
#  data_rho.astype(np.float32).tofile('rho_spcam_10x10.dat')
#  data_qd.astype(np.float32).tofile('q_spcam_10x10.dat')
#  data_ud.astype(np.float32).tofile('u_spcam_10x10.dat')
#  data_vd.astype(np.float32).tofile('v_spcam_10x10.dat')
  data_omegad.astype(np.float32).tofile('omega_spcam_10x10.dat')
print('Rank ', rank, ' all done!', time.asctime( time.localtime(time.time()) ))

