import csv
import numpy as np
import pickle
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import matplotlib as mpl

with open('/data/cloud/Goldtimes5/data/MCS_TRMM_3b42.pk', 'rb') as f:
    data= pickle.load(f)

Time= data[:,0]
Time= Time.astype(int)
L1= (Time <= 59*8)
L2= (Time > 334*8)
DJF= data[L1 | L2, :]

L1= (Time > 59*8)
L2= (Time <= 151*8)
MAM= data[L1 & L2, :]

L1= (Time > 151*8)
L2= (Time <= 243*8)
JJA= data[L1 & L2, :]

L1= (Time > 243*8)
L2= (Time <= 334*8)
SON= data[L1 & L2, :]

# Data x,y
Dlon= DJF[:,3]
Dlat= DJF[:,4]
Dlon= Dlon.astype(float)
Dlat= Dlat.astype(float)
#Dlon[Dlon<30]=Dlon[Dlon<30]+360
Dlon[Dlon>180]=Dlon[Dlon>180]-360

Mlon= MAM[:,3]
Mlat= MAM[:,4]
Mlon= Mlon.astype(float)
Mlat= Mlat.astype(float)
#Mlon[Mlon<30]=Mlon[Mlon<30]+360
Mlon[Mlon>180]=Mlon[Mlon>180]-360

Jlon= JJA[:,3]
Jlat= JJA[:,4]
Jlon= Jlon.astype(float)
Jlat= Jlat.astype(float)
#Jlon[Jlon<30]=Jlon[Jlon<30]+360
Jlon[Jlon>180]=Jlon[Jlon>180]-360

Slon= SON[:,3]
Slat= SON[:,4]
Slon= Slon.astype(float)
Slat= Slat.astype(float)
#Slon[Slon<30]=Slon[Slon<30]+360
Slon[Slon>180]=Slon[Slon>180]-360

Tlon= data[:,3]
Tlat= data[:,4]
Tlon= Tlon.astype(float)
Tlat= Tlat.astype(float)
#Tlon[Tlon<30]=Tlon[Tlon<30]+360
Tlon[Tlon>180]=Tlon[Tlon>180]-360

TR= data[:,5]
TR= TR.astype(float)

# plot

for i in ['DJF', 'MAM', 'JJA', 'SON', 'Total']:
    if i=='DJF':
        cb=np.array([56,108,176])/256
        lon= Dlon
        lat= Dlat
    elif i== 'MAM':
        cb=np.array([190,174,212])/256
        lon= Mlon
        lat= Mlat
    elif i== 'JJA':
        cb=np.array([240,2,127])/256
        lon= Jlon
        lat= Jlat
    elif i== 'SON':
        cb=np.array([253,192,134])/256
        lon= Slon
        lat= Slat
    elif i== 'Total':
        cb='b'
        lon= Tlon
        lat= Tlat
    plt.figure(figsize=(12,6))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude= 0))
    ax.coastlines()
#    print(np.amin(lon),' ',np.amax(lon),np.amin(lat),np.amax(lat))
    plt.scatter(lon,lat,1,marker='o',color=cb,transform=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -45, 45], crs=ccrs.PlateCarree())
#    ax.set_xlim(30, 389)
#    ax.set_ylim(-49, 49)
    picsave= '/data/cloud/Goldtimes5/pic/TRMM_3b42/Season/'+i+'.png'
    plt.savefig(picsave)
    plt.show()

larger= (TR>5.)
smaller= (TR<=5.)
TRl= TR[larger]
TRs= TR[smaller]
lonl= Tlon[larger]
lons= Tlon[smaller]
latl= Tlat[larger]
lats= Tlat[smaller]


plt.figure(figsize=(12,6))
cmap = plt.cm.jet  # define the colormap
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude= 0))
ax.coastlines()
cs1= plt.scatter(lons,lats,1,c=TRs,cmap=cmap,marker='o',transform=ccrs.PlateCarree())
cs2= plt.scatter(lonl,latl,1,c=TRl,cmap=cmap,marker='o',transform=ccrs.PlateCarree())
plt.colorbar()
plt.clim(2,7)
ax.set_extent([-180, 180, -49, 49], crs=ccrs.PlateCarree())
#ax.set_xlim(30, 389)
#ax.set_ylim(-49, 49)
picsave= '/data/cloud/Goldtimes5/pic/TRMM_3b42/Season/Rainfall.png'
plt.savefig(picsave)
plt.show()

