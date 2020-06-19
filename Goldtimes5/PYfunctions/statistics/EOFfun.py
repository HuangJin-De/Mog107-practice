# Input 2D data (time, space), area-weighted
# Output data_out {'EOF' / 'PC' / 'var'}

import numpy as np
import scipy.linalg as la

def EOF(data):
  N, K = data.shape
  max_mode = 100
  EOFs= np.zeros((max_mode, K))
  PCs= np.zeros((max_mode, N))
  var_exp= np.zeros((max_mode))

  data_trans= np.transpose(data)  
  COVdata= data_trans.dot(data) / (N-1)
  Evals, Evecs = la.eig(COVdata)
  Evals= Evals.real
  Evecs= Evecs.real

  var_e = Evals *100 / np.sum(Evals)

  idx = Evals.argsort()[::-1]
  Evals= Evals[idx]
  Evecs= Evecs[:,idx]

  for i in np.arange(0,max_mode):
    tempVec = Evecs[:, i].reshape((K,1))
    tempVal = Evals[i]
    PC      = data.dot(tempVec) / tempVal**0.5
    PCs[i,:]= PC.reshape((1,N))
    EOFs [i,:] = tempVec.reshape((1,K)) * tempVal**0.5
    var_exp[i] = var_e [i]

  out= {}
  out['EOF']= EOFs
  out['PC' ]= PCs
  out['var']= var_exp
  return out
