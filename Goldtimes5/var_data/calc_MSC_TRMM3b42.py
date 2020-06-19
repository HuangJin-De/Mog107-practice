import csv
import numpy as np
import pickle
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

years= np.arange(1998,2014)

path= '/data/cloud/yuyu/hw2/TRMMmask/TRMM_3hr_mask/'
# filename style= 'TRMM_3hr_1998_prop.txt'

# indt, id, radius[km,sqrt(area/pi)], lon[deg], lat[deg], average_pcp[mm/hr]
# 3hr, cloud id, r(km), lon, lat, mean rainfall (mm/hr)



for year in years:
    print('year: ',year)
    filename= path+'TRMM_3hr_'+str(year)+'_prop.txt'
    with open(filename, newline='') as f:
        data = pd.read_csv(f, delim_whitespace=True, header=None)
        data= np.array(data)
    if year == years[0]:
         data_all= data
    else:
         data_all= np.append(data_all,data,axis=0)
print('Data loaded')
print(data_all.shape)
data=data_all

print('Datasize: ',data.shape)
sz= data[0:,2]
sz= sz.astype(float)
critS= np.array(350) # 250 km
Sz= (sz>critS)

mr= data[0:,5]
mr= mr.astype(float)
mr_s= np.sort(mr)
topR= int(np.round(np.size(mr_s)*0.8))
critR= mr_s[topR]
Rain= (mr>critR)

data= data[Sz & Rain,:]

print('Datasize: ',data.shape)
with open('/data/cloud/Goldtimes5/data/MCS_TRMM_3b42.pk', 'wb') as fs:
     pickle.dump(data, fs)
