import Ngl
import os
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
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
         750, 700, 650, 600, 550, 500, 450, 400, 350, 300,
         250, 225, 200, 175, 150, 125, 100]
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

rh = np.zeros((indt.size, np.size(pint), np.size(latrange)-1, np.size(lonrange)-1))
for i in range(indt.size):
  ds = xr.open_dataset(filename[indt[i]])
  hyam = ds["hyam"]
  hybm = ds["hybm"]
  T    = ds["T"][:, :, latindex, lonindex]
#  z    = ds['Z3'][:, :, latindex, lonindex]
  q    = ds['Q'][:, :, latindex, lonindex]
  psrf = ds["PS"][:, latindex, lonindex]
  P0mb =  0.01*ds["P0"]
  #lat = ds["lat"]
  #lon = ds["lon"]

  intyp = 1
  kxtrp = False

  Tnew = Ngl.vinth2p(T,hyam,hybm,pint,psrf,intyp,P0mb,1,kxtrp)
  shape = np.shape(Tnew)
  qvs = np.empty(shape)
  qvs[Tnew==1e30] = np.NaN
  Tnew[Tnew==1e30] = np.NaN
  qnew = Ngl.vinth2p(q,hyam,hybm,pint,psrf,intyp,P0mb,1,kxtrp)
  qnew[qnew==1e30] = np.NaN

  tb0 = np.where(Tnew<273.15)
  ta0 = np.where(Tnew>=273.15)
  qvs[ta0[0], ta0[1], ta0[2], ta0[3]] = (380.*np.exp((17.27*(Tnew[ta0[0], ta0[1], ta0[2], ta0[3]]-273.15))
                                        /(Tnew[ta0[0], ta0[1], ta0[2], ta0[3]]-35.85)))/(np.array(pint)[ta0[1]]*100)
  qvs[tb0[0], tb0[1], tb0[2], tb0[3]] = (380.*np.exp((21.875*(Tnew[tb0[0], tb0[1], tb0[2], tb0[3]]-273.15))
                                        /(Tnew[tb0[0], tb0[1], tb0[2], tb0[3]]-7.65)))/(np.array(pint)[tb0[1]]*100)
  #qvs = (380.*np.exp((17.27*(Tnew-273.15))/(Tnew-35.85)))/(np.array(pint).reshape(1,npint,1,1)*100)
  for y in range(np.size(latrange)-1):
    for x in range(np.size(lonrange)-1):
      indy = np.where((lat[latindex] >= latrange[y]) &
                      (lat[latindex] <  latrange[y+1]))[0]
      indx = np.where((lon[lonindex] >= lonrange[x]) &
                      (lon[lonindex] <  lonrange[x+1]))[0]
      w = np.cos(lat[latindex][indy]*np.pi/180.)
      
#      ysize = np.size(indy)
#      xsize = np.size(indx)
#      shape = np.shape(Tnew)
#      qvs = np.empty((shape[0], shape[1], ysize, xsize))

#      Tfocus = Tnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1]
#      tb0 = np.where(Tfocus<273.15)
#      ta0 = np.where(Tfocus>=273.15)
#      qvs[ta0[0], ta0[1], ta0[2], ta0[3]] = (380.*np.exp((17.27*(Tfocus[ta0[0], ta0[1], ta0[2], ta0[3]]-273.15))
#					     /(Tfocus[ta0[0], ta0[1], ta0[2], ta0[3]]-35.85)))/(np.array(pint)[ta0[1]]*100)
#      qvs[tb0[0], tb0[1], tb0[2], tb0[3]] = (380.*np.exp((21.875*(Tfocus[tb0[0], tb0[1], tb0[2], tb0[3]]-273.15))
#					     /(Tfocus[tb0[0], tb0[1], tb0[2], tb0[3]]-7.65)))/(np.array(pint)[tb0[1]]*100)

      masked_qv = np.ma.masked_array(qnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1], np.isnan(qnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1]))
      masked_qvs = np.ma.masked_array(qvs[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1], np.isnan(qvs[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1]))
      qv_w = np.ma.average(masked_qv, axis = 2, weights = w)
      qvs_w = np.ma.average(masked_qvs, axis = 2, weights = w)
      tmp = 100*np.nanmean(qv_w, axis = 2)/np.nanmean(qvs_w, axis = 2)
       
      
      rh[i, :, y, x] = np.nanmean(tmp, axis = 0)
  
  print(i)
  #plt.plot(np.nanmean(tmp, axis = (0)), pint)
  #plt.ylim([1000,100])
  #plt.show()
  #print(np.shape(qvs))
  
print(rank, time.asctime( time.localtime(time.time()) ))

if rank == 0:
  data_rh = np.zeros((nstep, np.size(pint), np.size(latrange)-1, np.size(lonrange)-1))
else:
  data_rh = None
comm.Gather(rh, data_rh, root = 0)
if rank == 0:
  data_rh.astype(np.float32).tofile('rh_spcam_10x10.dat')

