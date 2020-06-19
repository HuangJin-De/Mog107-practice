import csv
import numpy as np
import numpy.matlib
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
#import cartopy.crs as ccrs
#import cartopy
#import sys
#sys.path.append('/data/cloud/Goldtimes5/code/')
#from PYfunctions.statistics.calc_areamean import areamean

data_path= '/data/cloud/Goldtimes5/data/GCM_SPCAM/CPL64/'
var_names= np.array(['VMSE_cor_inte','PRECT'])

mapfile= '/data/dadm1/model_output/SPCAM/CPL64/CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(mapfile,'r') as D:
    lat=D.variables['lat'][:]
    lon=D.variables['lon'][:]
    maxGY= lat.size
    maxGX= lon.size

weight= np.cos(lat*np.pi/180)
A= 510100000 * 10**6
A_w= A * weight / np.sum(weight)

#=== load data and plot
data_all={}
for vv in np.arange(0,var_names.shape[0]):
  var_name= var_names[vv]
  print(var_name)
  data_file= data_path+var_name+'.pk'
  with open(data_file, 'rb') as f:
    data= pickle.load(f)
    data_all[var_name]= data

#=== 
VTall= data_all['VMSE_cor_inte']
Prec= data_all['PRECT']
maxYr, maxD, maxY, maxX = Prec.shape

#=== Energy flux equator
EFE= np.zeros([maxYr, maxD])
for yy in np.arange(0, maxYr):
  for dd in np.arange(0, maxD):
    temp= VTall[yy,dd, 32:63] # 30S-30N
    VTmin= np.argmin(np.abs(temp)) + 32 # minimun point
    tempr= VTall[yy,dd, VTmin-3 : VTmin+4]
    latr= lat[ VTmin-3 : VTmin+4]
    EFE [yy,dd]= np.polyfit(tempr, latr, 1)[1]

filesave= data_path+ 'EFE_cor.pk'
with open(filesave, 'wb') as fs:
  pickle.dump(EFE, fs, protocol = 4)

#=== Prec center
Prec= np.mean(Prec,3)
Cent= np.zeros([maxYr, maxD])
for yy in np.arange(0, maxYr):
  for dd in np.arange(0, maxD):
    temp= Prec[yy,dd, 32:63] # 30S-30N
    temp= temp**10
    latr= lat[32:63]
    Cent[yy,dd] = np.sum(latr * temp)/np.sum(temp)

filesave= data_path+ 'Prec_Center.pk'
with open(filesave, 'wb') as fs:
  pickle.dump(Cent, fs, protocol = 4)
