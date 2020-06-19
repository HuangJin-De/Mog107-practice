import numpy as np
import math
import time
from mpi4py import MPI

print(time.asctime( time.localtime(time.time()) ))
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

t_initial = 1998
nstep = 18
dt_bounds = nstep//size
tbounds = np.array([rank*dt_bounds, (rank+1)*dt_bounds])+t_initial

indt = np.arange(tbounds[0],tbounds[-1])

lonrange = [45, 65, 85, 105, 125, 145, 165, 185, 205, 225]
latrange = [-30, -10, 10, 30]

def calc_GCD(lat1, lat2, lon1, lon2):
  """Great-Circle Distance in km"""
  r = 6371.
  lats = math.radians(lat1)
  latf = math.radians(lat2)
  lons = math.radians(lon1)
  lonf = math.radians(lon2)
  sigma = 2*math.asin(math.sqrt( math.sin(abs(latf-lats)/2)**2 + math.cos(lats)*math.cos(latf)*math.sin(abs(lonf-lons)/2)**2) )
  return r*sigma
print(calc_GCD(10, 30, 120, 140))
def calc_IP(r1, r2, d):
  """interaction potential"""
  if ((r1+r2) >= d):
    return 1
  else:
    return (r1+r2)/d
indt = np.arange(1998, 1999)
COP_all = np.zeros((indt.size, 365, 3, 9))
preci_mean = np.zeros((indt.size, 365, 3, 9))
for year in range(indt.size):
  f = open(r'/data/cloud/yuyu/hw2/TRMMmask/TRMM_daily_mask/TRMM_daily_{}_prop.txt'.format(indt[year]))
  text = []
  for line in f:
    line = line.strip().split()
    text.append(line)
  data = np.array(text, dtype=float)

  if ((indt[year]%4) != 0):
    for d in range(0,1):
      for y in range(np.size(latrange)-1):
        for x in range(np.size(lonrange)-1):
          ind = np.where((data[:,0] == d+1) & 
			 (data[:,4] >= latrange[y]) & 
			 (data[:,4] <  latrange[y+1]) & 
			 (data[:,3] >= lonrange[x]) & 
			 (data[:,3] <  lonrange[x+1]))
          cn = np.size(ind[0])
          if cn == 1:
            COP_all[year, d, y, x] = 0
            preci_mean[year, d, y, x] = data[ind[0][0], 5]
          elif cn != 0:
            COP = 0
            preci = data[ind[0][cn-1], 5]
            for i in range(cn-1):
              preci = preci + data[ind[0][i], 5]
              for j in range(i+1, cn):
                COP = COP + calc_IP(data[ind[0][i],2], data[ind[0][j],2], calc_GCD(data[ind[0][i],4], data[ind[0][j],4], data[ind[0][i],3], data[ind[0][j],3]))
            n = cn*(cn-1)/2
            COP_all[year, d, y, x] = COP/n
            preci_mean[year, d, y, x] = preci/cn
          else:
            COP_all[year, d, y, x] = 0.
            preci_mean[year, d, y, x] = 0.
  else:
    for d in range(0,59):
      for y in range(np.size(latrange)-1):
        for x in range(np.size(lonrange)-1):
          ind = np.where((data[:,0] == d+1) &
                         (data[:,4] >= latrange[y]) &
                         (data[:,4] <  latrange[y+1]) &
                         (data[:,3] >= lonrange[x]) &
                         (data[:,3] <  lonrange[x+1]))
          cn = np.size(ind[0])
          if cn == 1:
            COP_all[year, d, y, x] = 0
            preci_mean[year, d, y, x] = data[ind[0][0], 5]
          elif cn != 0:
            COP = 0
            preci = data[ind[0][cn-1], 5]
            for i in range(cn-1):
              preci = preci + data[ind[0][i], 5]
              for j in range(i+1, cn):
                COP = COP + calc_IP(data[ind[0][i],2], data[ind[0][j],2], calc_GCD(data[ind[0][i],4], data[ind[0][j],4], data[ind[0][i],3], data[ind[0][j],3]))
            n = cn*(cn-1)/2
            COP_all[year, d, y, x] = COP/n
            preci_mean[year, d, y, x] = preci/cn
          else:
            COP_all[year, d, y, x] = 0.
            preci_mean[year, d, y, x] = 0.

    for d in range(60,366):
      for y in range(np.size(latrange)-1):
        for x in range(np.size(lonrange)-1):
          ind = np.where((data[:,0] == d+1) &
                         (data[:,4] >= latrange[y]) &
                         (data[:,4] <  latrange[y+1]) &
                         (data[:,3] >= lonrange[x]) &
                         (data[:,3] <  lonrange[x+1]))
          cn = np.size(ind[0])
          if cn == 1:
            COP_all[year, d-1, y, x] = 0
            preci_mean[year, d-1, y, x] = data[ind[0][0], 5]
          elif cn != 0:
            COP = 0
            preci = data[ind[0][cn-1], 5]
            for i in range(cn-1):
              preci = preci + data[ind[0][i], 5]
              for j in range(i+1, cn):
                COP = COP + calc_IP(data[ind[0][i],2], data[ind[0][j],2], calc_GCD(data[ind[0][i],4], data[ind[0][j],4], data[ind[0][i],3], data[ind[0][j],3]))
            n = cn*(cn-1)/2
            COP_all[year, d-1, y, x] = COP/n
            preci_mean[year, d-1, y, x] = preci/cn
          else:
            COP_all[year, d-1, y, x] = 0.
            preci_mean[year, d-1, y, x] = 0.
  print(indt[year])
  print(time.asctime( time.localtime(time.time()) ))
'''
if rank == 0:
  data_COP = np.zeros((nstep, 365, 3, 9))
  data_preci = np.zeros((nstep, 365, 3, 9))
else:
  data_COP = None
  data_preci = None

comm.Gather(COP_all, data_COP, root=0)
comm.Gather(preci_mean, data_preci, root=0)

if rank == 0:
  data_COP.tofile('COP_w20.dat')
  data_preci.tofile('preci_w20.dat')
'''
