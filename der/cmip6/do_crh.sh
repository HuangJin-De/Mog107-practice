#!/usr/bin/bash

yr_s=1980
yr_e=2014
for model in $(ls cmip_2p5d/)
do
  #echo ${model}
  cp CODE/crh_freq.F  crh_${model}.f
  sed -i "s/model_name/${model}/g" crh_${model}.f
  sed -i "s/yr_s/${yr_s}/g" crh_${model}.f
  sed -i "s/yr_e/${yr_e}/g" crh_${model}.f
  mpifort -O3 -free crh_${model}.f 
  mpirun -n 1 ./a.out 
  rm -f crh_${model}.f a.out
  cp ./template/freq_template.ctl ./gs_ctl_files/${model}_freq.ctl
  sed -i "s/model/${model}/g" ./gs_ctl_files/${model}_freq.ctl

  cp CODE/crh_pdf.F  crh_${model}.f
  sed -i "s/model_name/${model}/g" crh_${model}.f
  sed -i "s/yr_s/${yr_s}/g" crh_${model}.f
  sed -i "s/yr_e/${yr_e}/g" crh_${model}.f
  mpifort -O3 -free crh_${model}.f
  mpirun -n 1 ./a.out
  rm -f crh_${model}.f a.out
  cp ./template/pdf_template.ctl ./gs_ctl_files/${model}_pdf.ctl
  sed -i "s/model/${model}/g" ./gs_ctl_files/${model}_pdf.ctl
done
