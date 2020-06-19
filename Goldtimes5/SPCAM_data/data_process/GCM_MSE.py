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

years= np.arange(6,7)
year_num= years[-1]- years[0] +1

SPCAM_path= '/data/dadm1/model_output/SPCAM/CPL64/'
# CPL64.cam.h0.0009-02-20-00000.nc , lat lon 96 144
# CPL64.cam.h1.0007-04-18-00000.nc
# LAT_15s_to_30n: 24, -14.21 : 1.89: 29.37
# LON_60e_to_180e: 49, 60:2.5:180

#=== load data
GCM_temp= SPCAM_path + 'CPL64.cam.h0.0001-01-01-00000.nc'
with Dataset(GCM_temp,'r') as D:
    lat=D.variables['lat'][:]
    lon=D.variables['lon'][:]
    lev=D.variables['lev'][:]
latr= lat[40:64]
lonr= lon[24:73]
maxGY= latr.size
maxGX= lonr.size

Cp= 240 * 4.186 # J/kg/K
Lv= 2.266*10**6 # J/kg
g= 9.8

data_all= np.zeros((year_num, 365, maxGY, maxGX)) # grid mean column condensed water (CCW)
for year in years:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print('Year ',str(year),", Current Time =", current_time)
    # numbers of days in each months
    day_num= [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for mon in np.arange(0,12,1):
        Nday= day_num[mon]
        for day in np.arange(0,Nday,1):
            cur_date= str(year).zfill(4)+'-'+str(mon+1).zfill(2)+'-'+str(day+1).zfill(2)
            cur_day= np.sum(day_num[0:mon])+day # Julian day, 0:364
            CRM_file= SPCAM_path+'CPL64.cam.h0.'+cur_date+'-00000.nc'
            with Dataset(CRM_file,'r') as D:
                T= D.variables['T'][:, -1, 40:64, 24:73]
                T= np.squeeze(np.mean(T,0))
                Q= D.variables['Q'][:, -1, 40:64, 24:73]
                Q= np.squeeze(np.mean(Q,0))
                Geo= D.variables['Z3'][:, -1, 40:64, 24:73]
                Geo= np.squeeze(np.mean(Geo,0))

            data= Cp*T + Lv*Q + g*Geo # near-surface MSE (975 hPa)
#            data2= Cp*T + Z # DSE
#            data3= Lv * Q # LHT
            data_all[year-years[0], int(cur_day), : ,:]= data

with open('/data/cloud/Goldtimes5/data/SPCAM/CPL64/MSE_sur.pk', 'wb') as fs:
     pickle.dump(data_all, fs, protocol = 4)

