#!/usr/bin/bash

yr_s=1979
yr_e=2015
for model in interim
do
  cp .code/jointpdf.F  crh_${model}.f
  sed -i "s/yr_s/${yr_s}/g" crh_${model}.f
  sed -i "s/yr_e/${yr_e}/g" crh_${model}.f
  mpifort -O3 -free crh_${model}.f -I/data/cloud/.local/include -L/data/cloud/.local/lib -lnetcdff -lnetcdf
  time mpirun -n 1 ./a.out 
done
