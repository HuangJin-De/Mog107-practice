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

cc.1=2
cc.2=4
cc.3=8
cc.4=9
cc.5=11
cc.6=3

i=1
while(i<=24)
"open ../gs_ctl_files/"model.i"_pdf.ctl"
i=i+1
endwhile

j=1
while(j<=4)

"set mproj off"

"set z 1"
"set x 1 100"
"set grads off"
"set parea 1.2 10 1.2 7.5"
"set xlabs 0|10|20|30|40|50|60|70|80|90|100"
"set vrange 0 0.03"
"set ylint 0.005"

"set cthick 14"
"set cmark 0"
"d f.24"

"set frame off"
"set xlab off"
"set ylab off"
"set grid off"

*"set line 0 1 13"
*"draw recf 1.3 4.4 3.9 7.4"
*"set line 1 1 13"
*"draw rec 1.3 4.4 3.9 7.4"

x=1.5
y=7.2
i=1
while (i<=6)

ii=6*(j-1)+i

if (ii<24)
"set cthick 13"
"set ccolor "cc.i""
"set cmark 0"
"set cstyle 2"
"d f."ii""

"set line "cc.i" 2 14"
"draw line "x" "y" "x+0.5" "y""
"set strsiz 0.16"
"set string 1 l 13 0"
"draw string "x+0.6" "y" "model.ii""

y=y-0.5
endif
i=i+1
say ii
endwhile


"draw title PDF"
"draw xlab column relative humidity (%)"
"printim cmip5_pdf_"j".png x2048 y1536"
"c"

"set frame on"
"set xlab on"
"set ylab on"
"set grid on"

j=j+1
endwhile

"quit"

