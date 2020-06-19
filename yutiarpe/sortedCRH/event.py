import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


#regime=[32.5,41.97778,45.03034,47.664032,50.807507,62.737694]
regime=np.array([[32.5,41.97778],[50.807507,62.737694]])
aggindex=np.fromfile('aggindex2.dat',np.float32)

def findymd(days):
  monthdays=      [31,28,31,30 ,31 ,30 ,31 ,31 ,30 ,31 ,30 ,31]
  mnday1=np.array([0 ,31,59,90 ,120,151,181,212,243,273,304,334])
  mnday2=np.array([31,59,90,120,151,181,212,243,273,304,334,365])+1
  year=(days-1)//365+1
  #print(days, days-365*(year-1))
  mon=np.where((days-365*(year-1)>=mnday1)*(days-365*(year-1)<mnday2))[0][0]+1
  day=days-365*(year-1)-mnday1[mon-1]
  return year, mon, day

fil=open('output.dat','wb')
for iregime in range(regime.shape[0]):
  data=[]
  regimeday=np.where((regime[iregime,0]<=aggindex)*(aggindex<=regime[iregime,1]))[0]+1
  for day in regimeday:
    a,b,c=findymd(day)
    #print(iregime, day, a, b, c)
    nc=Dataset('./CPL64_g/CPL64_mf_CRH_regime-%04d-%02d-%02d.nc'%(a,b,c))
    tmp=np.nanmean(nc.variables['mf'][:,:,:],axis=0)
    data.append(tmp)
  data=np.array(data)
  data.mean(axis=0).astype(np.float32).tofile(fil)
  data.std(axis=0).astype(np.float32).tofile(fil)
  print(iregime, regimeday.size)

    
    



