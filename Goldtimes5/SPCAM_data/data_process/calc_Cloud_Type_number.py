import csv
from datetime import datetime
import numpy as np
from netCDF4 import Dataset
import pickle
from scipy import interpolate
import matplotlib.pyplot as plt
#import xarray as xr
#import os

# Timer
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Start  =", current_time)

years= np.arange(6,11)
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

cL= 2500; cH= 6000; # critical height of low/mid/high clouds

num_all= np.zeros((year_num, 365, maxGY, maxGX, 6)) # ocurrance of each cloud type 
rain_all= np.zeros((year_num, 365, maxGY, maxGX, 6)) # cloud rain of each cloud type
CW_all= np.zeros((year_num, 365, maxGY, maxGX, 6)) # condensed water (CCW)
for year in years:
  now = datetime.now()
  current_time = now.strftime("%H:%M:%S")
  print('Year ',str(year),", Current Time =", current_time)
  for day in np.arange(0,365,1):
    data_file= '/data/cloud/alfalfa/wk13/OBJ/SPCAM_CRM_'+str(year).zfill(2)+'_'+str(day+1).zfill(3)+'_cobj_daily.dat'
    tmp = np.fromfile(data_file, dtype = np.float32)
    tmp= tmp.reshape((int(len(tmp)/10),10))
    c_top= tmp[:,6]; c_base= tmp[:,5]
    #=== masks
    mLow= (c_top < cL)
    mMid= (c_top < cH) & (c_base > cL)
    mHigh=(c_base > cH) 
    mSC = ((c_top < cH)&(c_top > cL)) & (c_base < cL)
    mDC = (c_top > cH) & (c_base < cL)
    mMH = (c_top > cH) & ((c_base < cH) & (c_base > cL))
    Low = tmp[mLow,:]; Mid = tmp[mMid,:]; High = tmp[mHigh,:];
    SC  = tmp[mSC, :]; DC  = tmp[mDC, :]; MH   = tmp[mMH,:];
    
    #=== Each grid point
    tmp_all= np.zeros((maxGY, maxGX, 6))
    tmp_rain_all = np.zeros((maxGY, maxGX, 6))
    tmp_CW_all = np.zeros((maxGY, maxGX, 6))
    for yy in np.arange(0,maxGY):
      for xx in np.arange(0,maxGX):
        tmp_y = tmp[:,2]; tmp_x = tmp[:,1];
        tmp_data= tmp[((tmp_y==yy) & (tmp_x==xx) ) ,:]
        c_top= tmp_data[:,6]; c_base= tmp_data[:,5]
        tmp_CW= tmp_data[:,8]; 
        tmp_rain= tmp_data[:,9] - tmp_data[:,8];
        #=== type mask
        mLow= (c_top < cL)
        mMid= (c_top < cH) & (c_base > cL)
        mHigh=(c_base > cH)
        mSC = ((c_top < cH)&(c_top > cL)) & (c_base < cL)
        mDC = (c_top > cH) & (c_base < cL)
        mMH = (c_top > cH) & ((c_base < cH) & (c_base > cL))
        #=== occurance
        Low = len(c_top[mLow]); Mid=len(c_top[mMid]); High= len(c_top[mHigh])
        SC  = len(c_top[mSC]);  DC =len(c_top[mDC]);  MH=   len(c_top[mMH])
        tmp_all[yy,xx,0]= Low; tmp_all[yy,xx,1]=Mid; tmp_all[yy,xx,2]=High
        tmp_all[yy,xx,3]= SC;  tmp_all[yy,xx,4]=DC;  tmp_all[yy,xx,5]=MH
        #=== CW
        Low_CW  = np.sum(tmp_CW [mLow]);
        Mid_CW  = np.sum(tmp_CW [mMid]);
        High_CW = np.sum(tmp_CW [mHigh]);
        SC_CW   = np.sum(tmp_CW [mSC]);
        DC_CW   = np.sum(tmp_CW [mDC]);
        MH_CW   = np.sum(tmp_CW [mMH]);
        tmp_CW_all[yy,xx,0]= Low_CW; tmp_CW_all[yy,xx,1]=Mid_CW; tmp_CW_all[yy,xx,2]=High_CW
        tmp_CW_all[yy,xx,3]= SC_CW;  tmp_CW_all[yy,xx,4]=DC_CW;  tmp_CW_all[yy,xx,5]=MH_CW
        #=== rain
        Low_rain  = np.sum(tmp_rain [mLow]);
        Mid_rain  = np.sum(tmp_rain [mMid]);
        High_rain = np.sum(tmp_rain [mHigh]);
        SC_rain   = np.sum(tmp_rain [mSC]);
        DC_rain   = np.sum(tmp_rain [mDC]);
        MH_rain   = np.sum(tmp_rain [mMH]);
        tmp_rain_all[yy,xx,0]= Low_rain; tmp_rain_all[yy,xx,1]=Mid_rain; tmp_rain_all[yy,xx,2]=High_rain
        tmp_rain_all[yy,xx,3]= SC_rain;  tmp_rain_all[yy,xx,4]=DC_rain;  tmp_rain_all[yy,xx,5]=MH_rain
    num_all  [year-years[0], day, : ,:, :]= tmp_all
    CW_all   [year-years[0], day, : ,:, :]= tmp_CW_all
    rain_all [year-years[0], day, : ,:, :]= tmp_rain_all

data_all={}
data_all['num'] = num_all
data_all['CW']  = CW_all
data_all['rain']= rain_all

filesave= '/data/cloud/Goldtimes5/data/SPCAM/CPL64/Cloud_Type.pk'
with open(filesave, 'wb') as fs:
     pickle.dump(data_all, fs, protocol = 4)

