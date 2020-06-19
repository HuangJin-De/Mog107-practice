import numpy as np
import math
import time
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as colors
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import make_axes_locatable

# create new pressure level
pint = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775,
         750, 700, 650, 600, 550, 500, 450, 400, 350, 300]#,
         #250, 225, 200, 175, 150, 125, 100]
npint = np.size(pint)
pnew = []
for (item1, item2) in zip(pint[:npint-1], pint[1:npint]):
  pnew.append((item1+item2)/2)
dp = []
for (item1, item2) in zip(pint[:npint-1], pint[1:npint]):
  dp.append(item1-item2)
# focus on the region
lonrange = [60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 181]
latrange = [-15, -5, 5, 15]
lonrange2= [60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]

cn = np.fromfile('../wk8/cvtcld_mean.dat', dtype = np.float32).reshape(3650, 3, 12)
omega = np.fromfile('../wk9/omega_spcam_10x10.dat', dtype = np.float32).reshape(3650, npint, 3, 12)
#crh = np.fromfile('../wk7/crh_spcam_10x10.dat', dtype = np.float32).reshape(3650, 3, 12)
#rh = np.fromfile('../wk7/rh_spcam_10x10.dat', dtype = np.float32).reshape(3650, npint2, 3, 12)
interval = [0, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.8]
'''
crh_cn = np.zeros(np.size(interval)-1)
sam = np.zeros(np.size(interval)-1)
for i in range(np.size(interval)-1):
  n = 0
  crhsum = 0
  for t in range(3650):
    #tmp = np.where((cn[t,:,:] >= interval[i])&(cn[t,:,:] < interval[i+1]), crh[t,:,:], np.nan)
    one = np.where((cn[t,:,:] >= interval[i])&(cn[t,:,:] < interval[i+1]), 1, np.nan)
    n = n+np.nansum(one)
    crhsum = crhsum+np.nansum(tmp)

  crh_cn[i] = crhsum/n
  sam[i] = n
  print('from',interval[i],'to',interval[i+1],'have',n)
'''
g = 9.81
vmf = -omega/g
cnlist = np.zeros([3650*3*12, 4])
n=0
for t in range(3650):
  for y in range(3):
    for x in range(12):
      cnlist[n, 0] = cn[t, y, x]
      cnlist[n, 1] = t
      cnlist[n, 2] = y
      cnlist[n, 3] = x
      n+=1
np.set_printoptions(threshold=np.inf)
cnlist = cnlist.tolist()
cnlist.sort(key = lambda s: s[0])
cnlist = np.array(cnlist, dtype = np.float32)

nnum = np.zeros(100)
sortnvmf = np.zeros([npint, 100])
for n in range(100):
  vmf_c = vmf[int(cnlist[n*1314, 1]), :, int(cnlist[n*1314, 2]), int(cnlist[n*1314, 3])].reshape(20,1)
  nnum[n] = cnlist[n*1314, 0]
  for i in range(1, 1314):
    index = cnlist[n*1314+i, :]
    vmf_c = np.concatenate((vmf_c, vmf[int(index[1]), :, int(index[2]), int(index[3])].reshape(20,1)), axis=1)
  sortnvmf[:, n] = np.nanmean(vmf_c, axis = 1)

#print(np.max(sortnvmf), np.min(sortnvmf))
nnum = np.round(nnum, 2)
fig, ax = plt.subplots()
cs = ax.pcolor(np.arange(100), pint, sortnvmf, cmap = plt.cm.get_cmap('bwr', 12), vmin = -0.015, vmax = 0.015)
ax.set_xticks(np.arange(0, 100, 10))
ax.set_xticklabels(nnum[::10])
ax.set_ylim([1000, 300])
ax.set_xlabel('Number of Convective Objects', fontsize = 14)
ax.set_ylabel('Pressure (hPa)', fontsize = 14)
ax.set_title('Vertical Mass Flux', fontsize = 16)
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
cbar = plt.colorbar(cs, cax=cax, ticks=np.arange(-0.015, 0.016, 0.005))
cbar.set_label('[kg m-2 s-1]')
plt.savefig('Vert_MF.png', bbox_inches = 'tight')
plt.show()

