import csv
from datetime import datetime
import numpy as np
from netCDF4 import Dataset
import pickle
from scipy import interpolate
import matplotlib.pyplot as plt

# Timer
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Start  =", current_time)

years= np.arange(2006,2016)
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
ERA5_temp= '/data/dadm1/reanalysis/ECMWF/ERA5/0.28125/hourly/PRS/2006/ERA5_hourly_p28125_deg_20060101.nc'
with Dataset(ERA5_temp,'r') as D:
    data= D.variables['latitude'][:]
    Elat= np.array(data)
    data= D.variables['longitude'][:]
    Elon= np.array(data)
    data= D.variables['level'][:]
    Elevel= data

lev= (Elevel==975)
Cp= 240 * 4.186 # J/kg/K
Lv= 2.266*10**6 # J/kg

first=1
for year in years:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print('Year ',str(year),", Current Time =", current_time)
    # numbers of days in each months
    if year % 4 !=0:
        day_num= [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    elif year % 100 == 0:
        day_num= [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
        day_num= [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for mon in np.arange(0,12,1):
        Nday= day_num[mon]
        for day in np.arange(0,Nday,1):
            cur_date= str(year)+str(mon+1).zfill(2)+str(day+1).zfill(2)
#            print(cur_date)
            ERA5_file= ERA5_path+str(year)+'/ERA5_hourly_p28125_deg_'+cur_date+'.nc'
            with Dataset(ERA5_file,'r') as D:
                T= D.variables['t'][0:-1:3, lev, :, :]
                T= np.squeeze(T)
                Q= D.variables['q'][0:-1:3, lev, :, :]
                Q= np.squeeze(Q)
                Z= D.variables['z'][0:-1:3, lev, :, :]
                Z= np.squeeze(Z)

            data= Cp*T + Z # DSE
            data2= Lv * Q # LHT
            # interpolate
            data_i= np.zeros((8, np.size(Tlat), np.size(Tlon) ))
            data_i2= np.zeros((8, np.size(Tlat), np.size(Tlon) ))
            for i in np.arange(0,8): # 8 times per day
                temp= data[i,:,:]
                f= interpolate.interp2d(Elon, Elat, temp, kind='cubic')
                temp_i= f(Tlon,Tlat)
                data_i [i,:,:]= temp_i
                temp2= data2[i,:,:]
                f2= interpolate.interp2d(Elon, Elat, temp2, kind='cubic')
                temp_i2= f2(Tlon,Tlat)
                data_i2 [i,:,:]= temp_i2
            if first ==1:
                data_all = data_i
                data_all2 = data_i2
            else:
                data_all= np.append(data_all,data_i,axis=0)
                data_all2= np.append(data_all2,data_i,axis=0)
            first=0

data=data_all
with open('/data/cloud/Goldtimes5/data/ERA5_MSEdry_interped.pk', 'wb') as fs:
     pickle.dump(data, fs, protocol = 4)
data=data_all2
with open('/data/cloud/Goldtimes5/data/ERA5_MSEmoist_interped.pk', 'wb') as fs:
     pickle.dump(data, fs, protocol = 4)

