# Calculation running smooth
# data [time, :]
# N (iteration steps)

import numpy as np

def smooth121(data,N):
  data= np.array(data,dtype= 'float')
  sz= data.shape
  if np.array(sz).size > 1 :
    if data.shape[1] == data.size:  # 1*N matrix
      print(data.size)
      data= data.reshape((data.size),1)
    else:
      data= data.reshape((sz[0], int(data.size/sz[0])))  
  else:
    data= data.reshape((data.size,1))
  
  for i in np.arange(0,N):
    temp= data.copy()
    temp[1:-2,:] = (2*data[1:-2,:] + data[0:-3,:] + data[2:-1,:])/4
    temp[0,:] = (2*data[0,:] + data[1,:] )/3
    temp[-1,:] = (2*data[-1,:] + data[-2,:] )/3
    data= temp  
  data= data.reshape(sz)
  return data

def running_mean(data,N):
  N= np.floor(N/2)
  N= int(N)
  data= np.array(data,dtype= 'float')
  sz= data.shape
  if np.array(sz).size > 1 :
    if data.shape[1] == data.size:  # 1*N matrix
      print(data.size)
      data= data.reshape((data.size),1)
    else:
      data= data.reshape((sz[0], int(data.size/sz[0])))
  else:
    data= data.reshape((data.size,1))

  temp= np.zeros(data.shape)
  for i in np.arange(N, sz[0]-N):
    temp[i,:] = np.nanmean(data[i-N: i+N+1,:], axis=0 )
  temp[0:N,:]=np.nan
  temp[-N:,:]= np.nan
  data= temp.reshape(sz)
  return data

def panta_mean(data):
  data= np.array(data,dtype= 'float')
  sz= data.shape
  szl= np.array(sz) # final shape
  if szl.size > 1 :
    if data.shape[1] == data.size:  # 1*N matrix
      data= data.reshape((data.size),1)
      szl[1]= szl[1]/5
    else:
      data= data.reshape((sz[0], int(data.size/sz[0])))
      szl[0]= szl[0]/5
  else:
    data= data.reshape((data.size,1))
    szl[0]= szl[0]/5

  temp= np.zeros((int(data.shape[0]/5), data.shape[1]))
  for i in np.arange(0, data.shape[0]/5, 1):
    i= int(i)
    temp[i,:] = np.nanmean(data[5*i: 5*(i+1),:], axis=0 )

  data= temp.reshape(szl)
  return data





