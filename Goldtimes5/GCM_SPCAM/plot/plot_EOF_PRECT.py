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
from PYfunctions.statistics.calc_areamean import areamean


data_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/Trop_mean/'
data_path_P= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/EOF_PRECT.pk'

data_file= data_path+'RCE.pk'
with open(data_file, 'rb') as f:
    data= pickle.load(f)

data_file= data_path_P
with open(data_file, 'rb') as f:
    pEOF= pickle.load(f)

plat= pEOF['lat']; plon= pEOF['lon']
Ny= len(plat); Nx= len(plon)

PC= np.zeros((100, 365*5))
Em= np.zeros(100)
PC_ori= pEOF['PC']
PC_std= np.std( PC_ori, axis= 1)
for i in np.arange(0,100):
  EOF= pEOF['EOF'][i,:].reshape((Ny,Nx))
  EOFm = areamean(EOF, plat, plon, [0,360,-30,30])
  Em [i] = EOFm
  PC[i,:]= PC_ori[i,:] * EOFm

RCE= data['RCE']
#RCEm= np.mean(RCE, axis=0); 
#RCEm = RCEm.reshape((1,RCEm.size)); RCEm= np.tile((RCE.shape[0],1))
RCEr= RCE.reshape(RCE.size)
cv= data['PRECT']+ data['SHFLX']; cvr= cv.reshape(cv.size)
Rad= data['Rad_cool']; Radr= Rad.reshape(Rad.size)
pr= data['PRECT']; prr= pr.reshape(pr.size)
LWs= data['FLNS']; LWsr= LWs.reshape(LWs.size)
SWs= data['FSNS']; SWsr= SWs.reshape(SWs.size)
LWt= data['FLNT']; LWtr= LWt.reshape(LWt.size)
SWt= data['FSNTOA']; SWtr= SWt.reshape(SWt.size)
prl= data['PRECT_largeVar']; prlr= prl.reshape(prl.size)

# Prec decompose
plt.figure(figsize=(8,3))
plt.plot(np.arange(1,5*365+1,1) ,np.zeros(RCEr.size),'k',linewidth=0.5)
plt.plot(np.arange(1,5*365+1,1) ,prr-np.mean(prr),'k-',linewidth=2, label= 'Precipitation')
plt.plot(np.arange(1,5*365+1,1) ,np.sum(PC[0:3,:],axis=0),'b-',linewidth=1, label= 'Sum PC1-3')
plt.plot(np.arange(1,5*365+1,1) ,np.sum(PC[:,:],axis=0),'r-',linewidth=1, label= 'Sum PC1-100')
ax=plt.gca()
ax.set_xlim(1, 366)
ax.set_ylim(-20, 30)
plt.legend(loc='upper right', fontsize=8)
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/PRECT_EOF/PRECT_PC_explain.png'
plt.savefig(picsave)
#plt.show()
plt.close()

plt.figure(figsize=(8,8))
plt.plot(np.arange(1,101) , abs(Em),'ro', linestyle='none')
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/PRECT_EOF/PRECT_EOF_tropmean.png'
plt.savefig(picsave)
#plt.show()
plt.close()

for i in np.arange(0, 10): 
  EOF= pEOF['EOF'][i,:].reshape((Ny,Nx)); PCp= PC_ori[i,:]; var= np.around(pEOF['var'][i], decimals=1)
  cmax= 200; cmin= -cmax; cinter= cmax/5; clevel= np.arange(cmin,cmax+cinter,cinter)
  EOF[(EOF> cmax)]=cmax; EOF[(EOF< cmin)]=cmin;
  plt.figure(figsize=(8,3))
  ax = plt.axes(projection=ccrs.PlateCarree())
  ax.coastlines()
  cs1= plt.contourf(plon,plat,EOF,levels=clevel ,cmap='RdYlBu_r')
  plt.text(130, 35, 'Var: '+str(var), fontsize= 12)
  plt.clim(cmin,cmax)
#  plt.colorbar().set_ticks(clevel)
  ax.set_xlim(-180, 180)
  ax.set_ylim(-30, 30)
  picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/PRECT_EOF/EOF'+str(i+1).zfill(2)+'.png'
  plt.savefig(picsave)
#  plt.show()
  plt.close()
  plt.figure(figsize=(8,3))
  plt.plot(np.arange(1, 6, 1/365), np.zeros(RCEr.size),'k',linewidth=0.5)
  plt.plot(np.arange(1, 6, 1/365) , PCp ,'b-',linewidth=1)
  ax=plt.gca(); ax.set_xlim(1, 6-1/365); ax.set_ylim(-3, 3)
  picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/PRECT_EOF/PC'+str(i+1).zfill(2)+'.png'
  plt.savefig(picsave)
#  plt.show()
  plt.close()

EOF= pEOF['EOF'][0,:].reshape((Ny,Nx)); 
cmax= 200; cmin= -cmax; cinter= cmax/5; clevel= np.arange(cmin,cmax+cinter,cinter)
EOF[(EOF> cmax)]=cmax; EOF[(EOF< cmin)]=cmin;
plt.figure(figsize=(8,3))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
cs1= plt.contourf(plon,plat,EOF,levels=clevel ,cmap='RdYlBu_r')
plt.clim(cmin,cmax)
plt.colorbar().set_ticks(clevel)
ax.set_xlim(-180, 180); ax.set_ylim(-30, 30)
picsave= '/data/cloud/Goldtimes5/pic/GCM_SPCAM/PRECT_EOF/EOF_colorbar.png'
plt.savefig(picsave)

