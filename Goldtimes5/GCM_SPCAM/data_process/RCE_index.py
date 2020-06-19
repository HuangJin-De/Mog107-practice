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
var_names= np.array(['PRECT','FLNT','FLNS','FSNTOA','FSDTOA','FSUTOA','FSNS','SHFLX'])

#Lv= 2.266*10**6 # J/kg
Lv= 2.477*10**6 # J/kg
rho= 1000 # kg/m3

#=== load data and plot
data_all={}
for vv in np.arange(0,var_names.shape[0]):
  var_name= var_names[vv]
  print(var_name)
  data_file= data_path+var_name+'.pk'
  with open(data_file, 'rb') as f:
    data= pickle.load(f)
    if var_name== 'PRECT':
      print('    Unit converted!!')
      data= data * rho * Lv # m/s * kg/m3 * J/kg = W/m2
    data_all[var_name]= data

#=== RCE index
Rad= data_all['FSNTOA'] - data_all['FSNS'] - data_all['FLNT'] + data_all['FLNS']
SW=  data_all['FSNTOA'] - data_all['FSNS']
LW=   - data_all['FLNT'] + data_all['FLNS']
RCE= data_all['PRECT'] + Rad + data_all['SHFLX']
data_all['Rad_cool']= Rad
data_all['Rad_cool_SW']= SW
data_all['Rad_cool_LW']= LW
data_all['RCE']= RCE

filesave= data_path+ 'RCE.pk'
with open(filesave, 'wb') as fs:
  pickle.dump(data_all, fs, protocol = 4)
