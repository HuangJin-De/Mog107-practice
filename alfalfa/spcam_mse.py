import Ngl
import os
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import glob
import matplotlib.colors as colors
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker

# create filename list
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

latindex = np.where((lat <= 30)&(lat >= -30))[0]
nlatindex = np.size(latindex)
lonindex = np.where((lon <= 225)&(lon >= 45))[0]
nlonindex = np.size(lonindex)
lonrange = [45, 65, 85, 105, 125, 145, 165, 185, 205, 225]
latrange = [-30, -10, 10, 30]

Cp = 1005.
Lv = 2.5E6
g = 9.81
# calculate MSE
'''
mse = np.zeros((3650, 26, nlatindex, nlonindex))
for i in range(np.size(filename)):
  ds = xr.open_dataset(filename[i])
  hyam = ds["hyam"]
  hybm = ds["hybm"]
  T    = ds["T"][:, :, latindex, lonindex]
  z    = ds['Z3'][:, :, latindex, lonindex]
  q    = ds['Q'][:, :, latindex, lonindex]
  psrf = ds["PS"][:, latindex, lonindex]
  P0mb =  0.01*ds["P0"]
  #lat = ds["lat"]
  #lon = ds["lon"]

  intyp = 1
  kxtrp = False 

#  Tnew = Ngl.vinth2p(T,hyam,hybm,pnew,psrf,intyp,P0mb,1,kxtrp)
#  Tnew[Tnew==1e30] = np.NaN
#  znew = Ngl.vinth2p(z,hyam,hybm,pnew,psrf,intyp,P0mb,1,kxtrp)
#  znew[znew==1e30] = np.NaN
#  qnew = Ngl.vinth2p(q,hyam,hybm,pnew,psrf,intyp,P0mb,1,kxtrp)
#  qnew[qnew==1e30] = np.NaN
  tmp = Cp*T+g*z+Lv*q
#  mse[i, :, :, :] = np.nanmean((Cp*Tnew+g*znew+Lv*qnew)/Cp, axis = 0)
  mse_new = Ngl.vinth2p(tmp,hyam,hybm,pnew,psrf,intyp,P0mb,1,kxtrp)
  mse_new[mse_new==1e30] = np.NaN
  mse[i, :, :, :] = np.nanmean(mse_new, axis = 0)
  print(mse[i, :, 0, 36])
  print(i)
'''
#mse.astype(np.float32).tofile('msenew_spcam.dat')

mse = np.fromfile('msenew_spcam.dat', dtype = np.float32).reshape(3650, 26, nlatindex, nlonindex)
cimse = np.nansum(mse[:,:23,:,:]*1E2*np.reshape(dp[:23], (1, npint-1-3, 1, 1)), axis = 1)/g
print(np.shape(np.where(cimse>0)))
print(np.shape(cimse))

lat = lat[latindex]
lon = lon[lonindex]
# split region

cimse_d = np.zeros((3650,3,9))

for y in range(np.size(latrange)-1):
  for x in range(np.size(lonrange)-1):
    indy = np.where((lat >= latrange[y]) &
                    (lat <= latrange[y+1]))[0]
    indx = np.where((lon >= lonrange[x]) &
                    (lon <= lonrange[x+1]))[0]
    w = np.cos(lat[indy]*np.pi/180.)
    a = cimse[:, indy[0]:indy[-1]+1, indx[0]:indx[-1]+1]
    ma = np.ma.MaskedArray(a, mask=np.isnan(a))
    tmp = np.ma.average(ma, axis = 1, weights = w)
    cimse_d[:, y, x] = np.nanmean(tmp, axis = 1)

cimse_d.astype(np.float32).tofile('cimse_discrete.dat')

#cimse = np.reshape(cimse, (10, 365, nlatindex, nlonindex)

# time series
'''
fig, ax = plt.subplots(figsize = (12,4), constrained_layout = True)
ax.plot(np.arange(3650), np.nanmean(cimse, axis = (1, 2))/1E9)
ax.set_xlim([0,3650])
ax.set_xticks(np.arange(0,3650,365))
ax.set_xticklabels(np.arange(10))
ax.set_xlabel('Year')
ax.set_ylabel('CIMSE (10^9 J/m-2)')
ax.set_title('Column Integrated MSE')
plt.savefig('cimse_ts1.png')
plt.show()
'''

lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()

fig = plt.figure(figsize=(7,7))
proj = ccrs.PlateCarree(central_longitude=135.)
ax = fig.add_subplot(1,1,1, projection=proj)
cs = ax.contourf(lon, lat, np.nanmean(cimse/1E9, axis = 0), transform = ccrs.PlateCarree(),
                                                        cmap = plt.get_cmap('jet',10),
                                                        extend = 'max')
ax.set_extent([45, 225, -30, 30], crs=ccrs.PlateCarree())
ax.set_xticks([45, 90, 135, 180, 225], crs=ccrs.PlateCarree())
ax.set_yticks([-30, -10, 10, 30], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labelsize=12)
ax.coastlines()
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(cs, cax=cax, ticks = np.arange(5,51,5), label = '(J/m-2)')
ax.set_title('10 years Mean CIMSE', fontsize = 16)
plt.savefig('cimse_m.png')
plt.show()





