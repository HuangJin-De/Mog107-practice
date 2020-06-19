import Ngl
import os
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import glob
import time
from mpi4py import MPI
import matplotlib.colors as colors
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker

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
         #250, 225, 200]#, 175, 150, 125, 100]
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

latindex = np.where((lat <= 25)&(lat >= -25))[0]
nlatindex = np.size(latindex)
lonindex = np.where((lon <= 190)&(lon >= 50))[0]
nlonindex = np.size(lonindex)
lonrange = [50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 181, 191]
latrange = [-25, -15, -5, 5, 15, 25]
#print(lon[lonindex])
#print(lat[latindex])

Cp = 1005.
Lv = 2.5E6
g = 9.81
#calculate crh in 60E-180E & 15S-15N (start)
qv = np.zeros((indt.size, npint, np.size(latrange)-1, np.size(lonrange)-1))
qvs = np.zeros((indt.size, npint, np.size(latrange)-1, np.size(lonrange)-1))
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
  #temp = Tnew
  #shape = np.shape(Tnew)
  #qvs = np.empty(shape)
  #qvs[temp==1e30] = np.NaN
  Tnew[Tnew==1e30] = np.NaN
  qnew = Ngl.vinth2p(q,hyam,hybm,pint,psrf,intyp,P0mb,1,kxtrp)
  qnew[qnew==1e30] = np.NaN

  for y in range(np.size(latrange)-1):
    for x in range(np.size(lonrange)-1):
      indy = np.where((lat[latindex] >= latrange[y]) &
                      (lat[latindex] <  latrange[y+1]))[0]
      indx = np.where((lon[lonindex] >= lonrange[x]) &
                      (lon[lonindex] <  lonrange[x+1]))[0]
      w = np.cos(lat[latindex][indy]*np.pi/180.)

      masked = np.ma.masked_array(Tnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1])
      tmp = np.ma.average(masked, axis = 2, weights = w)
      Tmean = np.nanmean(tmp, axis = 2)
      Tmean = np.nanmean(Tmean, axis = 0)
      qvs[i, :, y, x] = np.where((Tmean>=0), (380.*np.exp((17.27*(Tmean-273.15))/(Tmean-35.85)))/(np.array(pint)*100), \
                                             (380.*np.exp((21.875*(Tmean-273.15))/(Tmean-7.65)))/(np.array(pint)*100))

      masked = np.ma.masked_array(qnew[:, :, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1])
      tmp = np.ma.average(masked, axis = 2, weights = w)
      qmean = np.nanmean(tmp, axis = 2)
      qv[i, :, y, x] = np.nanmean(qmean, axis = 0)

  print(rank, i)


if rank == 0:
  data_qv = np.zeros((nstep, npint, np.size(latrange)-1, np.size(lonrange)-1))
  data_qvs= np.zeros((nstep, npint, np.size(latrange)-1, np.size(lonrange)-1))
else:
  data_qv = None
  data_qvs= None
comm.Gather(qv , data_qv , root = 0)
comm.Gather(qvs, data_qvs, root = 0)
if rank == 0:
  data_qv.astype(np.float32).tofile('qv_spcam_10x10_extend.dat')
  data_qvs.astype(np.float32).tofile('qvs_spcam_10x10_extend.dat')
print('Rank', rank, 'is done!', time.asctime( time.localtime(time.time()) ))

