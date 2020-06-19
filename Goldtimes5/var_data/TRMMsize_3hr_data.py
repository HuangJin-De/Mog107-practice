import csv
import numpy as np
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt

years= np.arange(2006,2016)
TRMM_path= '/data/dadm1/obs/TRMM/TRMM3B42size/' #TRMMsize_3hrs_2001.nc
# lon: 0-360, lat: -49-49

#=== load lat lon TRMM
TRMM_temp= '/data/dadm1/obs/TRMM/TRMM3B42/3B42.2004.3hr.nc'
with Dataset(TRMM_temp,'r') as D:
    data= D.variables['latitude'][:]
    Tlat_full= np.array(data)
    data= D.variables['longitude'][:]
    Tlon_full= np.array(data)
Tlat_mask= ( (Tlat_full > 0) & (Tlat_full <29.8) )
Tlon_mask= ( (Tlon_full > 100) & (Tlon_full <239.77) )
Tlat= Tlat_full[Tlat_mask]
Tlon= Tlon_full[Tlon_mask]

for year in years:
    print('Year ',year)
    # TRMM file
    TRMM_file= TRMM_path+ 'TRMMsize_3hrs_'+str(year)+'.nc'
    with Dataset(TRMM_file,'r') as D:
        TRMM= D.variables['objsize'][:,Tlat_mask,Tlon_mask] # correspond to t, Tlat, Tlon
    if year==years[0]:
        TRMM_all= TRMM
    else:
        TRMM_all = np.append(TRMM_all,TRMM,axis=0)

data= TRMM_all
with open('/data/cloud/Goldtimes5/data/TRMM3b42size.pk', 'wb') as fs:
     pickle.dump(data, fs, protocol = 4)
