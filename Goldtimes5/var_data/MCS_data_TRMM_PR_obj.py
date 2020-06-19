import csv
import numpy as np
import pickle
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

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

# find top 25% objects
mr= data_all[0:,1]
mr= mr.astype(float)
mr_s= np.sort(mr)
topR= int(np.round(np.size(mr_s)*0.8))
critR= mr_s[topR]
Rain= (mr>critR)

sz= data_all[0:,0]
sz= sz.astype(int)
critS= np.array(400) # 200 km
Sz= (sz>critS)


MCS= Rain & Sz
data= data_all[MCS, :]
# data= np.delete(data_all, np.where(notMCS), 0)

with open('../data/MCS_TRMM_PR2.pk', 'wb') as f:
     pickle.dump(data, f)
