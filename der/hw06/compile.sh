#!/usr/bin/bash

for run in CTRL64 PRSDAY64 #CTRL64
do
  cp .code/read3.F read_${run}.f
  sed -i "s/sed_run_name/$run/g" read_${run}.f
  mpifort -O3 -free read_${run}.f -I/data/cloud/.local/include -L/data/cloud/.local/lib -lnetcdff -lnetcdf
  time mpirun -n 1 ./a.out 
  rm -f read_${run}.f
  #cp .code/isen.ctl ${run}_isen.ctl
  #sed -i "s/run/$run/g" ${run}_isen.ctl
  #cp .code/map.ctl ${run}_map.ctl
  #sed -i "s/run/$run/g" ${run}_map.ctl
  cp .code/cloud.ctl ${run}_cloud.ctl
  sed -i "s/run/$run/g" ${run}_cloud.ctl
done
