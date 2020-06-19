import csv
import numpy as np
import pickle
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import matplotlib as mpl

with open('/data/cloud/Goldtimes5/data/MCS_TRMM_PR.pk', 'rb') as f:
    data= pickle.load(f)

# remove lat-biased point
lat= data[:,6]
lat= lat.astype(float)
lat= lat*0.1-40
NormalData= (lat<90)
#data= np.delete(data, np.where(lat>90), 0)
data= data[NormalData, :]

Julian= data[0:,8]
Julian= Julian.astype(int)
L1= (Julian  < 60)
L2= (Julian >= 335)
DJF= data[L1 | L2, :]

L1= (Julian >= 60)
L2= (Julian  < 152)
MAM= data[L1 & L2, :]

L1= (Julian >= 152)
L2= (Julian  < 244)
JJA= data[L1 & L2, :]

L1= (Julian >= 244)
L2= (Julian  < 335)
SON= data[L1 & L2, :]

# Data x,y
Dlon= DJF[:,5]
Dlat= DJF[:,6]
Dlon= Dlon.astype(float)
Dlat= Dlat.astype(float)
Dlon= Dlon*0.1
Dlat= Dlat*0.1+40-90
Dlon[Dlon<30]=Dlon[Dlon<30]+360

Mlon= MAM[:,5]
Mlat= MAM[:,6]
Mlon= Mlon.astype(float)
Mlat= Mlat.astype(float)
Mlon= Mlon*0.1
Mlat= Mlat*0.1+40-90
Mlon[Mlon<30]=Mlon[Mlon<30]+360

Jlon= JJA[:,5]
Jlat= JJA[:,6]
Jlon= Jlon.astype(float)
Jlat= Jlat.astype(float)
Jlon= Jlon*0.1
Jlat= Jlat*0.1+40-90
Jlon[Jlon<30]=Jlon[Jlon<30]+360

Slon= SON[:,5]
Slat= SON[:,6]
Slon= Slon.astype(float)
Slat= Slat.astype(float)
Slon= Slon*0.1
Slat= Slat*0.1+40-90
Slon[Slon<30]=Slon[Slon<30]+360

Tlon= data[:,5]
Tlat= data[:,6]
Tlon= Tlon.astype(float)
Tlat= Tlat.astype(float)
Tlon= Tlon*0.1
Tlat= Tlat*0.1+40-90
Tlon[Tlon<30]=Tlon[Tlon<30]+360

TR= data[:,1]
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
    plt.figure(figsize=(4.8,4.8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
#plt.scatter(Tlon,Tlat,1,marker='o',color='b')
    plt.scatter(lon,lat,1,marker='o',color=cb)
    ax.set_xlim(80, 150)
    ax.set_ylim(-45, 26)
    picsave= '/data/cloud/Goldtimes5/pic/TRMM_PR/'+i+'.png'
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


plt.figure(figsize=(4.8,4.8))
cmap = plt.cm.jet  # define the colormap
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
cs1= plt.scatter(lons,lats,1,c=TRs,cmap=cmap,marker='o')
cs2= plt.scatter(lonl,latl,1,c=TRl,cmap=cmap,marker='o')
plt.colorbar()
plt.clim(2,7)
ax.set_xlim(80, 150)
ax.set_ylim(-45, 26)
picsave= '/data/cloud/Goldtimes5/pic/TRMM_PR/Rainfall.png'
plt.savefig(picsave)
plt.show()



#for year in years:
#    print('year: ',year)
#    filename= file+'Rain_object_'+str(year)+'.csv'
#    with open(filename, newline='') as csvfile:
#        data = list(csv.reader(csvfile))
#        data= np.array(data)
#    lat= data[:,6]
#    lat= lat.astype(float)
#    ma= np.amax(lat)
#    mi= np.amin(lat)    
#    ma= ma*0.1+40-90
#    mi= mi*0.1+40-90
#    print('year: ',year, ', max= ', ma,', min= ',mi)
#

