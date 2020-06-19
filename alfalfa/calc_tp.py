import numpy as np
import xarray as xr
import glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as colors
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import make_axes_locatable

path = '/data/dadm1/model_output/SPCAM/CTRL64/'
filename = glob.glob(path+'CTRL64.cam.h0.00*-*-*-00000.nc')
filename = sorted(filename)

pint = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, \
         750, 700, 650, 600, 550, 500, 450, 400, 350, 300]
npint = np.size(pint)
dp = np.array(pint[:-1])-np.array(pint[1:])

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
lonrange = [60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 181]
latrange = [-15, -5, 5, 15]
lonrange2= [60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
lonloc = (np.array(lonrange2[:-1])+np.array(lonrange2[1:]))/2
latloc = (np.array(latrange[:-1])+np.array(latrange[1:]))/2

g = 9.81
q = np.fromfile('../wk11/qv_spcam_10x10_extend2.dat', dtype = np.float32).reshape(3650, npint, 5, 14)
q = q[:,:,1:4,1:13]
qs = np.fromfile('../wk11/qvs_spcam_10x10_extend2.dat', dtype = np.float32).reshape(3650, npint, 5, 14)
qs = qs[:,:,1:4,1:13]
u = np.fromfile('../wk9/u_spcam_10x10.dat', dtype = np.float32).reshape(3650, npint, np.size(latrange)-1, np.size(lonrange)-1)
v = np.fromfile('../wk9/v_spcam_10x10.dat', dtype = np.float32).reshape(3650, npint, np.size(latrange)-1, np.size(lonrange)-1)
cn= np.fromfile('../wk8/cvtcld_mean.dat', dtype = np.float32).reshape(3650, np.size(latrange)-1, np.size(lonrange)-1)
#crh = np.fromfile('../wk7/crh_spcam_10x10.dat', dtype = np.float32).reshape(3650, 3, 12)
#tmp = np.fromfile('../wk10/crh_2lev_spcam_10x10.dat', dtype = np.float32).reshape(3650, 2, 3, 12)
crh = np.zeros([2, 3650, 3, 12])
crh[0, :, :, :] = 100*np.nansum(q[:,:11,:,:]*dp[:11].reshape(1, 11, 1, 1), axis = 1)/np.nansum(qs[:,:11,:,:]*dp[:11].reshape(1, 11, 1, 1), axis = 1)
crh[1, :, :, :] = 100*np.nansum(q[:,11:-1,:,:]*dp[11:].reshape(1, 8, 1, 1), axis = 1)/np.nansum(qs[:,11:-1,:,:]*dp[11:].reshape(1, 8, 1, 1), axis = 1)
#crh[0, :, :, :] = np.nansum(q[:,:11,:,:], axis = 1)
#crh[1, :, :, :] = np.nansum(q[:,11:,:,:], axis = 1)
interval = [0, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.8]

print(np.nansum(q[:,:11,:,:]*dp[:11].reshape(1, 11, 1, 1), axis = 1)[0,:,:])

uq = u*q
vq = v*q

utp = np.zeros([2, 3650, 3, 12])
vtp = np.zeros([2, 3650, 3, 12])
# 1000 to 300
#utp = np.nansum(((uq[:,:-1,:,:]+uq[:,1:,:,:])/2*dp.reshape(1,npint-1,1,1)*100), axis = 1)/g	#[kg/m/s] [3650,3,12]
#vtp = np.nansum(((vq[:,:-1,:,:]+vq[:,1:,:,:])/2*dp.reshape(1,npint-1,1,1)*100), axis = 1)/g	#[kg/m/s] 
# 1000 to 700
utp[0, :, :, :] = np.nansum(((uq[:,:11,:,:]+uq[:,1:12,:,:])/2*dp[:11].reshape(1,11,1,1)*100), axis = 1)/g  #[kg/m/s]
vtp[0, :, :, :] = np.nansum(((vq[:,:11,:,:]+vq[:,1:12,:,:])/2*dp[:11].reshape(1,11,1,1)*100), axis = 1)/g  #[kg/m/s]
# 700 to 300
utp[1, :, :, :] = np.nansum(((uq[:,11:-1,:,:]+uq[:,12:,:,:])/2*dp[11:].reshape(1,8,1,1)*100), axis = 1)/g  #[kg/m/s]
vtp[1, :, :, :] = np.nansum(((vq[:,11:-1,:,:]+vq[:,12:,:,:])/2*dp[11:].reshape(1,8,1,1)*100), axis = 1)/g  #[kg/m/s]

tp = np.zeros([2, 5, 3, 12])
utp_mean = np.zeros([2, 5, 3, 12])
vtp_mean = np.zeros([2, 5, 3, 12])
# all-time mean
utp_mean[0, 0, :, :] = np.mean(utp[0,:,:,:], axis = 0)
vtp_mean[0, 0, :, :] = np.mean(vtp[0,:,:,:], axis = 0)
utp_mean[1, 0, :, :] = np.mean(utp[1,:,:,:], axis = 0)
vtp_mean[1, 0, :, :] = np.mean(vtp[1,:,:,:], axis = 0)
tp[0, 0, :, :] = np.sqrt((utp_mean[0,0,:,:])**2+(vtp_mean[0,0,:,:])**2)
tp[1, 0, :, :] = np.sqrt((utp_mean[1,0,:,:])**2+(vtp_mean[1,0,:,:])**2)
#print(np.max(tp), np.min(tp))

# seasonal mean
for i in range(2):
  utps = utp[i,:,:,:].reshape(10, 365, 3, 12)
  vtps = vtp[i,:,:,:].reshape(10, 365, 3, 12)

  utp_mean[i, 1, :, :] = np.mean(utps[:,59:151,:,:], axis = (0, 1))
  vtp_mean[i, 1, :, :] = np.mean(vtps[:,59:151,:,:], axis = (0, 1))
  tp[i, 1, :, :] = np.sqrt((utp_mean[i, 1, :, :])**2+(vtp_mean[i, 1, :, :])**2)
  utp_mean[i, 2, :, :] = np.mean(utps[:,151:243,:,:], axis = (0, 1))
  vtp_mean[i, 2, :, :] = np.mean(vtps[:,151:243,:,:], axis = (0, 1))
  tp[i, 2, :, :] = np.sqrt((utp_mean[i, 2, :, :])**2+(vtp_mean[i, 2, :, :])**2)
  utp_mean[i, 3, :, :] = np.mean(utps[:,243:334,:,:], axis = (0, 1))
  vtp_mean[i, 3, :, :] = np.mean(vtps[:,243:334,:,:], axis = (0, 1))
  tp[i, 3, :, :] = np.sqrt((utp_mean[i, 3, :, :])**2+(vtp_mean[i, 3, :, :])**2)
  utp_s = np.concatenate((utps[:, 334:, :, :], utps[:, :59, :, :]), axis = 1)
  vtp_s = np.concatenate((vtps[:, 334:, :, :], vtps[:, :59, :, :]), axis = 1)
  utp_mean[i, 4, :, :] = np.mean(utp_s, axis = (0, 1))
  vtp_mean[i, 4, :, :] = np.mean(vtp_s, axis = (0, 1))
  tp[i, 4, :, :] = np.sqrt((utp_mean[i, 4, :, :])**2+(vtp_mean[i, 4, :, :])**2)
print(np.max(tp), np.min(tp))

crh_mean = np.zeros([2, 5, 3, 12])
crh_mean[0, 0, :, :] = np.nanmean(crh[0, :, :, :], axis = 0)
crh_mean[1, 0, :, :] = np.nanmean(crh[1, :, :, :], axis = 0)

# crh seasonal mean
for i in range(2):
  crhs = crh[i, :, :, :].reshape(10, 365, 3, 12)
  crh_s = crhs[:, 59:151, :, :]
  crh_mean[i, 1, :, :] = np.nanmean(crh_s, axis = (0, 1))
  crh_s = crhs[:, 151:243, :, :]
  crh_mean[i, 2, :, :] = np.nanmean(crh_s, axis = (0, 1))
  crh_s = crhs[:, 243:334, :, :]
  crh_mean[i, 3, :, :] = np.nanmean(crh_s, axis = (0, 1))
  crh_s = np.concatenate((crhs[:, 334:, :, :], crhs[:, :59, :, :]), axis = 1)
  crh_mean[i, 4, :, :] = np.nanmean(crh_s, axis = (0, 1))
  print(np.max(crh_mean[i,:,:,:], axis = (1,2)), np.min(crh_mean[i,:,:,:], axis = (1,2)))
'''
for i in range(np.size(interval)-1):
  tmp = np.where((cn>=interval[i])&(cn<interval[i+1]))
  print('From', interval[i], 'to', interval[i+1], ':',np.shape(tmp)[1])
print(np.shape(tmp))
'''
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()

# low[all, MAM, JJA, SON, DJF], mid[all, MAM, JJA, SON, DJF]

cmapnum = np.array([[15,20,20,20,20],[15,20,20,20,20]])
vmin = np.array([[70,65,65,65,65],[35,20,20,20,20]])
vmax = np.array([[85,85,85,85,85],[65,70,70,70,70]])
tpscale = np.array([10,1])
quiverkeylabel = np.array(['100 [kg m-1 s-1]','10 [kg m-1 s-1]'])
quiverkeypos = np.array([0.81, 0.825])
title1 = np.array(['Low','Mid'])
title2 = np.array(['(1000hPa-700hPa)', '(700hPa-300hPa)'])
title3 = np.array(['ANN', 'MAM', 'JJA', 'SON', 'DJF'])

for i in range(2):
  for j in range(5):
    fig = plt.figure(figsize=(12,4))
    proj = ccrs.PlateCarree(central_longitude=120.)
    ax = fig.add_subplot(1,1,1, projection=proj)
    cs = ax.pcolor(lonrange2, latrange, crh_mean[i,j,:,:], transform = ccrs.PlateCarree(), \
                   cmap = plt.cm.get_cmap('YlGnBu', cmapnum[i,j]), vmin = vmin[i,j], vmax = vmax[i,j], zorder=1)
    cs2= ax.quiver(lonloc, latloc, utp_mean[i,j,:,:]/tpscale[i], vtp_mean[i,j,:,:]/tpscale[i], transform = ccrs.PlateCarree(), \
                   color = [255/255, 75/255, 75/255], units = 'xy', width = 0.5, scale = 2.5, zorder=3)
    ax.quiverkey(cs2, quiverkeypos[i], 1.05, 10, label = quiverkeylabel[i],labelpos = 'E', fontproperties={'size':14})#, coordinates = 'figure')
    ax.set_extent([60, 180, -15, 15], crs=ccrs.PlateCarree())
    ax.set_xticks(lonrange2, crs=ccrs.PlateCarree())
    ax.set_yticks(latrange, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=16)
    ax.coastlines(zorder=2)
    ax.grid(ls=':')
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    cbar = plt.colorbar(cs, cax=cax, ticks=np.arange(vmin[i,j],vmax[i,j]+1,5))
    cbar.ax.get_yaxis().labelpad = 15
    cbar.set_label('(%)')
    ax.set_title('Water Vapor Transport & RH at '+title1[i]+'-Level '+title2[i]+'\n'+title3[j], fontsize = 18)
    plt.savefig('tpcrh_'+title1[i]+'_'+title3[j]+'.png')


