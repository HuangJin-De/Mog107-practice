#!/usr/bin/bash

path=/data/der0318/CMIP6/
target=/data/der0318/work/cmip6/cmip_ln/
template=/data/der0318/work/cmip6/template/
reoutput=/data/der0318/work/cmip6/cmip_re/

cd ${target}
#rm -rf *

for model in KACE-1-0-G  #$(ls ${path})
do
if [ "$(echo ${model} | cut -c 1-1)" == "F" ]
then
  echo ${model}
elif [ "$(echo ${model} | cut -c 1-1)" == "A" ]
then
  echo ${model}
else
  echo ${model} "true"
  mkdir ${model}
  cd ${model}
  for var in ta hus zg wap
  do
    # link and make ctl
    ctlname=${model}_${var}.ctl
    cp ${template}ctl_template.ctl ${ctlname}
    
    days=0
    for filename in $(ls ${path}${model}/${var}_*)
    do
      precut=$(echo ${filename} | rev | cut -f1 -d'_' | rev )

      year1=$(echo ${precut} | cut -c 1-4)
      if [ "$days" == "0" ]
      then
        sed -i "s/year/${year1}/g" ${ctlname}
      fi
      year2=$(echo ${precut} | cut -c 10-13)

      lnname=${var}_${model}_${year1}_${year2}.nc

      ln -s ${filename} ${lnname}
      
      # for ctl
      prenccut=$(ncdump -h ${lnname} | grep time)
 
      n=$(echo ${prenccut} | cut -f2 -d'(' | cut -c 1-5)
      if [ "${model}" == "KACE-1-0-G" ]
      then
        n=$(echo ${prenccut} | cut -f2 -d'=' | cut -c 2-6)
      elif [ "${model}" == "MIROC6" ]
      then
        n=$(echo ${prenccut} | cut -f2 -d'(' | cut -c 1-3)
      fi

      echo "CHSUB $((${days}+1)) $((${days}+${n})) ${lnname}" >> ${model}_${var}.ctl
      days=$((${days}+${n}))
    done
    
    sed -i "s/totalday/${days}/g" ${ctlname}   
    # 365_day_calendar or not
    day365=$(($days%365))
    echo $days $day365
    if [ "$day365" == "0" ]
    then
      sed -i "s/template/template 365_day_calendar/g" ${ctlname}
    fi

  done
  
  #cd ${reoutput}
  #rm -rf ${model}
  #mkdir ${model}
  #cd ${model}

  ## regrid 
  #cp ${template}grid_2p5d.dat .
  #cp ${template}grid_2p5d.ctl .
  #if [ "$day365" == "0" ]
  #then
  #  sed -i "s/template/template 365_day_calendar/g" grid_2p5d.ctl
  #else 
  #  sed -i "s/template/zrev/g" grid_2p5d.ctl
  #fi
  #cp ${template}regrid_template.gs ${model}_regrid.gs
  #sed -i "s/modelname/${model}/g" ${model}_regrid.gs
  #grads -blc ${model}_regrid.gs
  #cp ${template}ctl_2p5_template.ctl ${model}.ctl 
  #sed -i "s/modelname/${model}/g" ${model}.ctl
  #cd ${target}
fi
done


