import csv
import numpy as np
from netCDF4 import Dataset
import pickle
import matplotlib.pyplot as plt

years= np.arange(2006,2016)
TRMM_path= '/data/dadm1/obs/TRMM/TRMM3B42size/' #TRMMsize_3hrs_2001.nc
# lon: 0-360, lat: -49-49
ERA5_path= '/data/dadm1/reanalysis/ECMWF/ERA5/0.28125/hourly/PRS/' #/yr/ERA5_hourly_p28125_deg_20150411.nc
# lon: 100-240, lat: 0-30

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


#=== load lat lon ERA5
ERA5_temp= '/data/dadm1/reanalysis/ECMWF/ERA5/0.28125/hourly/PRS/2015/ERA5_hourly_p28125_deg_20150411.nc'
with Dataset(ERA5_temp,'r') as D:
    data= D.variables['latitude'][:]
    Elat= np.array(data)
    data= D.variables['longitude'][:]
    Elon= np.array(data)
    data= D.variables['level'][:]
    Elevel= data
lev= (Elevel==500)

for year in years
    print('Year ',year)
    # numbers of days in each months
    if year % 4 !=0:
        day_num= [31 28 31 30 31 30 31 31 30 31 30 31]
    elif year % 100 == 0:
        day_num= [31 28 31 30 31 30 31 31 30 31 30 31]
    else:
        day_num= [31 29 31 30 31 30 31 31 30 31 30 31]

    # TRMM file
    TRMM_file= TRMM_path+ 'TRMMsize_3hrs_'+str(year)+'.nc'
    with Dataset(TRMM_file,'r') as D:
        TRMM= D.variables['objsize'][:,Tlat_mask,Tlon_mask] # correspond to t, Tlat, Tlon
    Time= np.arange(1,TRMM.shape[0])
    

    for tt in Time
        #  calculate current time
        Jday= np.ceil(tt,8)

    for mon in np.arange(1,13,1):
        


        for day in 























    if year==years[0]:
        TRMM_all= TRMM
    else:
        TRMM_all = np.append(TRMM_all,TRMM,axis=0)


