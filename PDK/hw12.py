nc = Dataset('/data/dadm1/obs/CERES_SYN1deg-Day_Terra-Aqua-MODIS/CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_2000.nc')
print(nc.variables)
top_flux = nc.variables['toa_net_all_daily'][:]

sfc_flux = nc.variables['ini_sfc_sw_down_all_daily'][:] + nc.variables['ini_sfc_sw_up_all_daily'][:]\
        + nc.variables['ini_sfc_lw_down_all_daily'][:] + nc.variables['ini_sfc_lw_up_all_daily'][:]

ATM_rad = flux - sfc_flux
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 20:11:09 2020

@author: Kevin
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


#enso = np.loadtxt('/data/cloud/PDK/ENSO/index.txt')
nc = Dataset('/data/dadm1/obs/CERES_SYN1deg-Day_Terra-Aqua-MODIS/CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_2000.nc')
print(nc.variables)
top_flux = nc.variables['toa_net_all_daily'][:]

sfc_flux = nc.variables['ini_sfc_sw_down_all_daily'][:] + nc.variables['ini_sfc_sw_up_all_daily'][:]\
	+ nc.variables['ini_sfc_lw_down_all_daily'][:] + nc.variables['ini_sfc_lw_up_all_daily'][:]

ATM_rad = flux - sfc_flux





sy = 1979
sy_idx = sy - 1950
enso_fl = enso[sy_idx:sy_idx+10,1:].flatten()
ym = np.arange(sy,sy+10,1/12)
#####################################
'''
j = int(0)
EVENT = -12*enso_fl
for ef in enso_fl:
    if ef >=0.5:
        EVENT[j] = 2
    elif ef <=0.5:
        EVENT[j] = 0
    else:
        EVENT[j] = 1
    j = int(j+1)
'''
'''
plt.figure(figsize=(4, 12))
plt.plot(enso_fl,ym,'k')
#plt.plot(enso_fl[enso_fl>=0.5],ym[enso_fl>=0.5],'r.')
#plt.plot(enso_fl[enso_fl<=-0.5],ym[enso_fl<=-0.5],'b.')
plt.plot(ym*0+0.5,ym,'r--')
plt.plot(ym*0-0.5,ym,'b--')
'''
years = np.arange(sy,sy+11,1)
'''
plt.yticks(years)
plt.ylim(years[0],years[-1])
plt.xlim(-2.5,2.5)
'''
#ax.set_ylim([y_ticks[0] - 0.5, y_ticks[-1] + 0.5])

#plt.savefig(str(sy) + '.png')
#np.savetxt('ENSO_index.txt', enso[29:,:], delimiter=' ')   
#####################################################
crh_pdf_mean = np.zeros((12,100))
crh_pdf_enso_mean = np.zeros((12,100))
crh_pdf_NINO_mean = np.zeros((12,100))

for j in range(1,int(crh_pdf.shape[0]/12)-2):
    #print((enso[j,1]+enso[j-1,-1] + enso[j,-2]))
    if (enso[j,1]+enso[j-1,-1] + enso[j-1,-2]) > 3:
        crh_pdf_enso_mean = crh_pdf_enso_mean + np.mean(crh_pdf[j*12-2:j*12+1,:],axis=0)
        print(j+1979,'ENSO')
    else:
        if (enso[j,1]+enso[j-1,-1] + enso[j-1,-2]) < -3:
            crh_pdf_NINO_mean = crh_pdf_NINO_mean + np.mean(crh_pdf[j*12-2:j*12+1,:],axis=0)
            print(j+1979,'NINO')
        else:
            crh_pdf_mean = crh_pdf_mean + np.mean(crh_pdf[j*12-2:j*12+1,:],axis=0)

crh_pdf_mean = crh_pdf_mean/np.mean(crh_pdf_mean,axis=1)[0]
crh_pdf_enso_mean = crh_pdf_enso_mean/np.mean(crh_pdf_enso_mean,axis=1)[0]
crh_pdf_NINO_mean = crh_pdf_NINO_mean/np.mean(crh_pdf_NINO_mean,axis=1)[0]
######################################
# plot

'''
plt.figure(figsize=(5.5,5.5))
plt.title('all')
plt.contourf(crh_pdf_mean,cmap = plt.cm.jet)
plt.colorbar()


plt.figure(figsize=(5.5,5.5))
plt.title('ONI>1.5')
plt.contourf(crh_pdf_enso_mean,cmap = plt.cm.jet)
plt.colorbar()
'''
plt.figure(figsize=(8,6))
plt.title('CRH frequency in DJF')
plt.plot(np.arange(1,101,1),crh_pdf_mean[0,:],'k',label='all',linewidth=1)
plt.plot(np.arange(1,101,1),crh_pdf_enso_mean[0,:],'r',label='ONI>1.0',linewidth=1)
plt.plot(np.arange(1,101,1),crh_pdf_NINO_mean[0,:],'b',label='ONI<-1.0',linewidth=1)
plt.legend()
plt.show()

'''
ST_enso_fl = enso_fl[np.where(enso_fl[1]>1.5)]
enso_fl[enso_fl>1.5]
'''
