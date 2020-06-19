#!/usr/bin/bash

for run in ec_era #CTRL64 #PRSDAY64 #CTRL64
do
  cp .code/jointpdf.F read_${run}.f
  sed -i "s/sed_run_name/$run/g" read_${run}.f
  mpifort -O3 -free read_${run}.f -I/data/cloud/.local/include -L/data/cloud/.local/lib -lnetcdff -lnetcdf
  time mpirun -n 1 ./a.out 
  rm -f read_${run}.f
done
