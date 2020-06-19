"reinit"

path="/data/dadm1/reanalysis/ECMWF/ITM/daily/"

"sdfopen "path"landsea_mask.nc"
"open /data/cloud/der/hw07/.code/grid.ctl"

"set lon 0 357.5"
"set lat -88.75 88.75"

"set undef -99999."
"set fwrite /data/cloud/der/hw07/for_you/daily_interim_lsm_2p5d.dat"
"set gxout fwrite"

"d lterp(lsm.1,grid.2(z=1,t=1),aave)"

"disable fwrite"

"quit"
