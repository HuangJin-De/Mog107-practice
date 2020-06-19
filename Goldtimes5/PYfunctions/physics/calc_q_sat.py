# Output saturated water vapor mixing ratio (kg/kg) under T, P
# Input: T(K), P(hPa)

import numpy as np

def q_sat(T,P):
  es= (1.0007+(3.46e-6*P))*6.1121*np.exp(17.502*(T-273.15)/(240.97+(T-273.15)))
  ws= 0.622*es/(P-es)
  return ws
