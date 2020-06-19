
# =================================================================== #
#   import packages
import numpy as np
import xarray as xr
import Ngl
import matplotlib.pyplot as plt
import matplotlib.colors as mcls
from global_land_mask import globe
from matplotlib import cm
import sys
from mpi4py import MPI

# =================================================================== #
#   call MPI
def getMPI():
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()
  return rank, size
mpirank,mpisize=getMPI()

# =================================================================== #
#   start program

homedir = '/data/dadm1/model_output/SPCAM/'
casename = 'CPL64'
fhead = casename+'.cam.h0.00'
ftail = '-00000.nc'

year = np.arange(6)+1
mon = np.arange(12)+1
day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])

#   parameter setting
cp = 1005.
lv = 2.5e6
g = 9.8
a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i=\
 6.11147274, 0.503160820, 0.188439774e-1,\
 0.420895665e-3, 0.615021634e-5,0.602588177e-7,\
 0.385852041e-9, 0.146898966e-11, 0.252751365e-14
a0,a1,a2,a3,a4,a5,a6,a7,a8=\
 6.11239921, 0.443987641, 0.142986287e-1, \
 0.264847430e-3, 0.302950461e-5, 0.206739458e-7,\
 0.640689451e-10,-0.952447341e-13,-0.976195544e-15


#   example file
fname_ex = homedir+casename+'/'+fhead+'01-01-01'+ftail
fex = xr.open_dataset(fname_ex)
lato = fex['lat']
#mask_dgr = (lato>=-45.)*(lato<=45.)
mask_dgr = (lato>=-30.)*(lato<=30.)
lat = lato[mask_dgr]
lon = fex['lon']
hyam = fex['hyam']
hybm = fex['hybm']
p0 = fex['P0']*1e-2

plevo = np.arange(1e2,1e3+10.,50.)
#plevo = np.array([900,800,700,600,500,400,300,200,100,50])
plevo = xr.DataArray(plevo, coords=[('lev', plevo)])
plevo.attrs['long_name'] = 'pressure level'
plevo.attrs['units'] = 'hPa'
dplev = plevo[1:].data-plevo[:-1].data
clat = np.cos(lat*np.pi/180)
clat_4d = np.tile(clat.data.reshape(1,1,clat.size,1), (24,plevo.size,1,lon.size))

lon2 = lon.data.copy()
lon2[lon2>180.] = lon2[lon2>180.]-360.
lon2d, lat2d = np.meshgrid(lon2, lat.data)
#ocean = ~(globe.is_land(lat2d, lon2d))
#ocean2 = ocean.copy()
#ocean2[1:,:] = ocean2[1:,:]*ocean[:-1,:]
#ocean2[:-1,:] = ocean2[:-1,:]*ocean[1:,:]
#ocean2[:,1:] = ocean2[:,1:]*ocean[:,:-1]
#ocean2[:,:-1] = ocean2[:,:-1]*ocean[:,1:]
#ocean2[1:,1:] = ocean2[1:,1:]*ocean[:-1,:-1]
#ocean2[1:,:-1] = ocean2[1:,:-1]*ocean[:-1,1:]
#ocean2[:-1,:-1] = ocean2[:-1,:-1]*ocean[1:,1:]
#ocean2[:-1,1:] = ocean2[:-1,1:]*ocean[1:,:-1]
ocean2=np.ones(lon2d.shape)*True

reglow = 0
reghgh = 100
regnum = 100
dreg = (reghgh-reglow)/regnum	#   dreg = 1.6e7

regime = np.arange(reglow,reghgh,dreg)
regime = xr.DataArray(regime, coords=[('regime', regime)])
regime.attrs['long_name'] = 'CRH regime'
regime.attrs['units'] = '%'
'''
reglow = 0.
reghgh = 70.
regnum = 50.
dreg = (reghgh-reglow)/regnum	#   dreg = 1.4
'''

#for indf in np.arange(mpirank*filelist.size/mpisize,(mpirank+1)*filelist.size/mpisize,dtype=np.int):
#for yri in range(year.size):
#for yri in [0]:
yri=mpirank
#yri=0
#print()
for moni in range(mon.size):
#for moni in [0]:
    #moni=mpirank
    #print('mon: %02d'%mon[moni])
    for dayi in range(day[moni]):
    #for dayi in [1]:
        print(mpirank,moni,dayi+1)
        #print(dayi+1,mpirank)
        fname = homedir+casename+'/'+fhead+'%02d-'%year[yri]+'%02d-'%mon[moni]+'%02d'%(dayi+1)+ftail
        f = xr.open_dataset(fname)
    
        time = f['time']
        t = f['T'][:,:,mask_dgr,:]     #K
        q = f['Q'][:,:,mask_dgr,:]     #kg/kg
        p = f['PRES'][:,:,mask_dgr,:]  #Pa
        ps = f['PS'][:,mask_dgr,:]     #hPa
        ps4d = np.tile(ps.data.reshape(time.size,1,lat.size,lon.size),(1,t.shape[1],1,1))
        t = np.where(p<=ps4d, t, np.nan)
        q = np.where(p<=ps4d, q, np.nan)
    
        #   column-integrated MSE
        #Flatau formulation:
        dt       = np.where(t-273.16>-80,t-273.16,-80)
        polysvp1 = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*\
                   (a4i+dt*(a5i+dt*(a6i+dt*(a7i+\
                   a8i*dt)))))))
        polysvp1 = polysvp1*100.
        tmp1      = 0.622*polysvp1/np.where(p-polysvp1>1e-3,p-polysvp1,1e-3)
        #       Flatau formulation:, vapor
        polysvp1 = a0 + dt*(a1+dt*(a2+dt*(a3+dt*\
                   (a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
        polysvp1 = polysvp1*100.
        tmp2      = 0.622*polysvp1/np.where(p-polysvp1>1e-3,p-polysvp1,1e-3)
        qs       = np.where(t>273.15,tmp2,tmp1)  #kg/kg


        

        p100 = np.where(p>=1e4,p,np.nan)
        cq = np.nansum((q[:,1:,:,:]+q[:,:-1,:,:])/2*(p100[:,1:,:,:]-p100[:,:-1,:,:])/g,axis=1)
        cqs = np.nansum((qs[:,1:,:,:]+qs[:,:-1,:,:])/2*(p100[:,1:,:,:]-p100[:,:-1,:,:])/g,axis=1)
        crh = cq/cqs*100
    #    cimse = np.trapz(mse,p,axis=1)/g
        del t, q, qs, cq, cqs, dt, tmp1, tmp2, polysvp1

        #   veritcal intp on pressures. OMEGA
        omega = f['OMEGA'][:,:,mask_dgr,:]
        omega = Ngl.vinth2p(omega,hyam,hybm,plevo,ps,1,p0,1,False)
        omega = np.where(omega>=1e30,np.nan,omega)
        omega = xr.DataArray(omega, coords=[('time',time),('lev',plevo),('lat',lat),('lon',lon)])
        
        #   sorting
        omega_reg = np.zeros((24,plevo.size,int(regnum)))
        reg_wgtbase = omega_reg.copy()
        for regi in range(int(regnum)):
            ocean_3d = np.tile(ocean2.reshape(1,lat.size,lon.size), (time.size,1,1))
            mask_reg = (crh>reglow+(regi*dreg))*(crh<=reglow+((regi+1)*dreg))*ocean_3d
            mask_reg4d = np.tile(mask_reg.reshape(time.size,1,lat.size,lon.size), (1,plevo.size,1,1))
            tmp = np.nansum(np.where(mask_reg4d,omega*clat_4d,np.nan),axis=(2,3))
            tmp_wgtbase = np.nansum(np.where(mask_reg4d,clat_4d,np.nan),axis=(2,3))
            omega_reg[:,:,regi] += tmp
            reg_wgtbase[:,:,regi] += tmp_wgtbase
    
        f.close(); del f
      
        omega_reg = -1*omega_reg/g
    
        omega_reg = xr.DataArray(omega_reg, coords=[('time',np.arange(24)),('lev', plevo), ('regime', regime)])
        omega_reg.attrs['long_name'] = 'Vertical mass flux sorted by CRH regime'
        omega_reg.attrs['units'] = 'kg/m2/s'
        reg_wgtbase = xr.DataArray(reg_wgtbase, coords=[('time',np.arange(24)),('lev', plevo), ('regime', regime)])
        reg_wgtbase.attrs['long_name'] = 'Weighted base of CRH regime'
        reg_wgtbase.attrs['units'] = '1'
        
        outdir = '/data/cloud/yuyu/hw6/dataCRH/'+casename+'_g'
        fout = outdir+'/'+casename+'_mf_CRH_regime-00%02d'%year[yri]+'-%02d-%02d.nc'%(mon[moni],dayi+1)
        ds = xr.Dataset({'mf': omega_reg, 'reg_wgtbase': reg_wgtbase, 'lev': plevo, 'regime': regime,'time':np.arange(24)})
        ds.to_netcdf(path=fout, mode='w')
        del ds
    
            
            



