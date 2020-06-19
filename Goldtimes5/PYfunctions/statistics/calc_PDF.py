# Calculation areamean
# areamnea (data, LAT, LON, region )
# LAT -90~90
# LON 0~360
# region[ minLON, maxLON, minLAT, maxLAT ]

import numpy as np

def PDF(data, level):
  data=    np.array(data)
  level=   np.array(level)
  data=    data.reshape(data.size)
  N=       level.size +1
  profile= np.zeros(N)
  #=== 1st
  mask=    (data< level[0])
  profile[0]= len(data[mask])
  #=== inter
  for i in range(1,N-1):
    upper=  level[i]; lower= level[i-1]
    mask=   (data>= lower) & (data< upper)
    profile[i] = len(data[mask])
  #=== last
  mask=         (data>= level[-1])
  profile[N-1]= len(data[mask])

  profile= profile/ data.size
  return profile
