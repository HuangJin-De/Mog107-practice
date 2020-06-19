import csv
from datetime import datetime
import numpy as np
import numpy.matlib
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
CRM_temp= SPCAM_path + 'CPL64.cam.h1.0001-01-01-00000.nc'
with Dataset(GCM_temp,'r') as D:
    ilev= D.variables['ilev'][:]

dp= np.diff(ilev)
dp= dp[2:] # lowest 24 levels of GCM grid
dp2 = dp.reshape((dp.size,1))
dp2 = numpy.matlib.repmat(dp2,1,64)

with Dataset(CRM_temp,'r') as D:
    lat=D.variables['LAT_15s_to_30n'][:]
    lon=D.variables['LON_60e_to_180e'][:]
    maxGY= lat.size
    maxGX= lon.size

# parameters
g= 9.8
rho= 1000

#data_all= np.zeros((year_num, 365, maxGY, maxGX)) # grid mean column condensed water (CCW)
#data_all2= np.zeros((year_num, 365, maxGY, maxGX)) # grid varaiation of CCW
data_all3= np.zeros((year_num, 365, maxGY, maxGX)) # mean distance of top 10% CCW columns
for year in years:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print('Year ',str(year),", Current Time =", current_time)
    # numbers of days in each months
    day_num= [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for mon in np.arange(0,12):
        Nday= day_num[mon]
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print('     ',str(mon+1),", Current Time =", current_time)
        for day in np.arange(0,Nday):
            cur_date= str(year).zfill(4)+'-'+str(mon+1).zfill(2)+'-'+str(day+1).zfill(2)
            cur_day= np.sum(day_num[0:mon])+day # Julian day, 0:364
            CRM_file= SPCAM_path+'CPL64.cam.h1.'+cur_date+'-00000.nc'
            with Dataset(CRM_file,'r') as D:
                data= D.variables['CRM_QC_LON_60e_to_180e_LAT_15s_to_30n']
                data= np.squeeze(np.mean(data,0))
                for gy in np.arange(0,maxGY):
                    for gx in np.arange(0,maxGX):
                        temp= np.squeeze(data[:,:,gy,gx])
                        #print(data.shape)
                        temp= temp * dp2 / np.sum(dp)
                        hm= np.sum(temp,0)
                        ind= hm.argsort()[-6:][::-1]
                        hm_min6= np.amin(np.sort(hm)[-6:])
                        dist= np.zeros(15)
                        ni=0
                        if hm_min6==0:
                          data_all3[year-years[0], int(cur_day), gy, gx]=np.nan
                        else:
                          for n in np.arange(0,6):
                            for d in np.arange(n+1,6):
                                dist[ni] = np.amin([abs(ind[n]-ind[d]), 64-abs((ind[n]-ind[d]))])
                                ni += 1 
#                        data_all[year-years[0], int(cur_day), gy, gx]= np.mean(hm)
#                        data_all2[year-years[0], int(cur_day), gy, gx]= np.var(hm,ddof=1)    
                          data_all3[year-years[0], int(cur_day), gy, gx]= np.mean(dist)*4

#with open('/data/cloud/Goldtimes5/data/SPCAM/CPL64/CCW_gridmean.pk', 'wb') as fs:
#     pickle.dump(data_all, fs, protocol = 4)

#with open('/data/cloud/Goldtimes5/data/SPCAM/CPL64/CCW_gridvar.pk', 'wb') as fs:
#     pickle.dump(data_all2, fs, protocol = 4)

with open('/data/cloud/Goldtimes5/data/SPCAM/CPL64/CCW_aggrgate.pk', 'wb') as fs:
     pickle.dump(data_all3, fs, protocol = 4)

