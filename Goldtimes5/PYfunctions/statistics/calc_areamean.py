# Calculation areamean
# areamnea (data, LAT, LON, region )
# LAT -90~90
# LON 0~360
# region[ minLON, maxLON, minLAT, maxLAT ]

import numpy as np

def areamean(data,LAT,LON,region):
  [mLON,mLAT]=np.meshgrid(LON,LAT)
  mask=np.ones(mLON.shape,dtype=bool);
  mask[(mLON< region[0])]=0; 
  mask[(mLON> region[1])]=0; 
  mask[(mLAT< region[2])]=0; 
  mask[(mLAT> region[3])]=0;
  mask[np.isnan(data)]=0;
  weight=np.cos(mLAT* np.pi/180);
  data=data*weight;
  area=np.nansum(np.nansum(data[mask]))/np.nansum(np.nansum(weight[mask]));
  return area
