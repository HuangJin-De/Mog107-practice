# load processed SPCAM data
# SPCAM(varname, typename)
# Input variable name, type name ('ann','DJF','JJA')


import numpy as np
import pickle
from netCDF4 import Dataset

def SPCAM(varname, typename):
  data_path= '/data/cloud/Goldtimes5/data/SPCAM/CPL64/'

  data_file= data_path+varname+'.pk'
  with open(data_file, 'rb') as f:
    data= pickle.load(f)

  if typename== 'ann':
    return data
  elif typename== 'DJF':
    mask= np.zeros(365)
    mask[0:59]=1;    mask[334:]=1;
    mask= (mask==1)
    data= data[:,mask,:,:]
    return data
  elif typename== 'JJA':
    mask= np.zeros(365)
    mask[151:243]=1
    mask= (mask==1)
    data= data[:,mask,:,:]
    return data

