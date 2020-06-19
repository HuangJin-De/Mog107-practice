"reinit"
"set display color white"
"c"

"set lwid 13 8.5"
"set lwid 14 14.0"
"set annot 1 13"
"set strsiz 0.19"
"set xlopts 1 13 0.19"
"set ylopts 1 13 0.19"
"set clopts 1 13 0.19"
"set rgb 200 100 100 100 220"
"set grid on 3 200 6"
"set map 1 1 12"

model.1="ACCESS-CM2"
model.2="AWI-ESM-1-1-LR"
model.3="BCC-CSM2-MR"
model.4="BCC-ESM1"
model.5="CanESM5"
model.6="CESM2"
model.7="CESM2-FV2"
model.8="CESM2-WACCM"
model.9="CESM2-WACCM-FV2"
model.10="EC-Earth3"
model.11="GFDL-CM4"
model.12="GFDL-ESM4"
model.13="INM-CM4-8"
model.14="INM-CM5-0"
model.15="IPSL-CM6A-LR"
model.16="MIROC6"
model.17="MPI-ESM-1-2-H"
model.18="MPI-ESM1-2-HR"
model.19="MPI-ESM1-2-LR"
model.20="MRI-ESM2-0"
model.21="NorESM2-LM"
model.22="NorESM2-MM"
model.23="TaiESM1"
model.24="ERA-I"

i=1
while(i<=24)
"open ../gs_ctl_files/"model.i"_freq.ctl"
i=i+1
endwhile

"set lon 0 360"
"set lat -30 30"

i=1
while (i<=24)
"set grads off"
"set parea 1 10 4.3 7.3"
"set ylint 30"
"color 0 1 0.1 -kind white->white->burlywood->sienna"
"set gxout grfill"
"d f."i""
"xcbar 10.1 10.2 2 6.5 -fs 2 -ft 13 -fw 0.15 -fh 0.15 -line on"
"set string 1 l 13 0"
"set strsiz 0.18"
"draw string 10.3 2 %"
"draw title "model.i""

ii=i+1
"set grads off"
"set parea 1 10 1.2 4.2"
"color 0 1 0.1 -kind white->white->burlywood->sienna"
"set gxout grfill"
"d f."ii""
"draw title "model.ii""

iii=ii/2
"printim dry_"iii".png x2048 y1536"
"c"

i=i+2
endwhile

"quit"

