import csv
import numpy as np
import pickle
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import matplotlib as mpl

years= np.arange(1998,2014)
file= '/data/dadm1/obs/TRMM_PR/PR_object/'

for year in years:
    print('year: ',year)
    filename= file+'Rain_object_'+str(year)+'.csv'
    with open(filename, newline='') as csvfile:
        data = list(csv.reader(csvfile))
        data= np.array(data)
#        sz= data[0:,0]
#        mr= data[0:,1]
#    print(sz.shape)
#    print(mr.shape)
#    sz= sz.reshape((np.size(sz), 1)
#    mr= mr.reshape((np.size(mr), 1))
    if year == years[0]:
 #       sz_all =sz
 #       mr_all =mr
         data_all= data
    else:
#        sz_all= np.append(sz_all,sz,axis=0)
#        mr_all= np.append(mr_all,mr,axis=0)
         data_all= np.append(data_all,data,axis=0)
print('Data loaded')
print(data_all.shape)

#==== for loop for each time
years= np.arange(1998,2014)
days= np.arange(1,366)
hours= np.arange(1,25)
data_yr= data_all[:,9].astype(int)
data_day= data_all[:,8].astype(int)
data_hr= data_all[:,7].astype(int)

#for year in years[2]
#  yr_mask= (data_yr == year)
#  for day in days[0]
#    day_mask= (data_day == day)
#    for hour in hours[0]
#      hr_mask= (data_hr == hour)
#      cur_time= yr_mask & day_mask & hr_mask
#      data= data_all [cur_time,:]
   
      

year=years[2]
yr_mask= (data_yr == year)
day = days[0]
day_mask= (data_day == day)
hour = hours[0]
hr_mask= (data_hr == hour)
print('Yr: ',year,', day: ',day, 'hr:', hour)
cur_time= yr_mask & day_mask & hr_mask
data= data_all [cur_time,:]

sub_data= data
#=== N
N= sub_data.shape[0]

#=== distance






















