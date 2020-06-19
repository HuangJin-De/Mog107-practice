import csv
import numpy as np
import numpy.matlib
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt
from scipy.stats import t as t_dist
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy
import sys
sys.path.append('/data/cloud/Goldtimes5/code/')
from PYfunctions.statistics.smooth import smooth121 as sm
from PYfunctions.statistics.smooth import running_mean as rm

#=== map
mapfile= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(mapfile,'r') as D:
    lat=D.variables['lat'][:]
    lon=D.variables['lon'][:]
    maxY= lat.size
    maxX= lon.size

data_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/'

data_file= data_path+'VMSE_cor_inte.pk'
with open(data_file, 'rb') as f:
    VT= pickle.load(f)

data_file= data_path+'EFE_cor.pk'
with open(data_file, 'rb') as f:
    EFE= pickle.load(f)

VT1= VT[0,0,:]
VT2= VT[0,182,:]
VTm= np.mean(np.mean(VT,axis=1),axis=0)

print(EFE[0,0])
print(EFE[0,182])

#  5 yr
plt.figure(figsize=(8,3))
plt.plot(lat ,np.zeros(lat.size),'k',linewidth=0.5)
plt.plot(lat , VT1/10**15 ,'b',linewidth=1, label= '1/1')
plt.plot(lat , VT2/10**15 ,'k',linewidth=1, label= '7/1')
plt.plot(lat , VTm/10**15 ,'r',linewidth=1, label= 'mean')
plt.plot( EFE[0,0], 0, 'bo')
plt.plot( EFE[0,182], 0, 'ko')
ax=plt.gca()
ax.set_xlim(-90, 90)
ax.set_ylim(-6, 7)
plt.legend(loc='upper right', fontsize=8)
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/EFE_PC/VMSE.png'
plt.savefig(picsave)
plt.show()


