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
cc.2=2
cc.3=2
cc.4=2
cc.5=4
cc.6=4
cc.7=4
cc.8=4
cc.9=4
cc.10=4
cc.11=9
cc.12=9
cc.13=9
cc.14=9
cc.15=9
cc.16=9
cc.17=3
cc.18=3
cc.19=3
cc.20=3
cc.21=3
cc.22=3
cc.23=2
cc.24=1

cm.1=2
cm.2=4
cm.3=7
cm.4=8
cm.5=1
cm.6=2
cm.7=4
cm.8=7
cm.9=8
cm.10=10
cm.11=1
cm.12=2
cm.13=4
cm.14=7
cm.15=8
cm.16=10
cm.17=1
cm.18=2
cm.19=4
cm.20=7
cm.21=8
cm.22=10
cm.23=1
cm.24=3

"open ../gs_ctl_files/multimodel.ctl"

"set grads off"
"set parea 1.2 7.5 1.2 7.5"

x=8
y=8

i=2
while (i<=24)
ii=i-1
"set t "i" "i""
"set x 1 2"
"set gxout scatter"
"set vrange 22 42"
"set vrange2 -10 10"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d dry;wet*1000"
if (i=2)
"set frame off"
"set xlab off"
"set ylab off"
"set grid off"
endif

"set line "cc.ii" 1 14"
"draw mark "cm.ii" "x" "y" 0.25"
"set strsiz 0.125"
"set string 1 l 13 0"
"draw string "x+0.2" "y" "model.ii""

y=y-0.35

i=i+1
endwhile

ii=24
"set t 0 1"
"set gxout scatter"
"set vrange 22 42"
"set vrange2 -10 10"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d dry;wet*1000"

"draw xlab delta CRH % (wet-dry)"
"draw ylab delta peak 0.1% (dry-wet)"

"printim multimodel_dcrh_dpeak.png x2048 y1536"
"c"

"set frame on"
"set xlab on"
"set ylab on"
"set grid on"
"set grads off"
"set parea 1.2 7.5 1.2 7.5"

x=8
y=8

i=2
while (i<=24)
ii=i-1
"set t "i" "i""
"set x 1 2"
"set gxout scatter"
"set vrange 22 42"
"set vrange2 -20 50"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d dry;adry*100-(100-adry*100)"
if (i=2)
"set frame off"
"set xlab off"
"set ylab off"
"set grid off"
endif

"set line "cc.ii" 1 14"
"draw mark "cm.ii" "x" "y" 0.25"
"set strsiz 0.125"
"set string 1 l 13 0"
"draw string "x+0.2" "y" "model.ii""

y=y-0.35

i=i+1
endwhile

ii=24
"set t 0 1"
"set gxout scatter"
"set vrange 22 42"
"set vrange2 -20 50"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d dry;adry*100-(100-adry*100)"

"draw xlab delta CRH % (wet-dry)"
"draw ylab delta area % (dry-wet)"

"printim multimodel_dcrh_darea.png x2048 y1536"
"c"


"set frame on"
"set xlab on"
"set ylab on"
"set grid on"
"set grads off"
"set parea 1.2 7.5 1.2 7.5"

x=8
y=8

i=2
while (i<=24)
ii=i-1
"set t "i" "i""
"set x 1 2"
"set gxout scatter"
"set vrange 50 70"
"set vrange2 -20 50"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d inter;adry*100-(100-adry*100)"
if (i=2)
"set frame off"
"set xlab off"
"set ylab off"
"set grid off"
endif

"set line "cc.ii" 1 14"
"draw mark "cm.ii" "x" "y" 0.25"
"set strsiz 0.125"
"set string 1 l 13 0"
"draw string "x+0.2" "y" "model.ii""

y=y-0.35

i=i+1
endwhile

ii=24
"set t 0 1"
"set gxout scatter"
"set vrange 50 70"
"set vrange2 -20 50"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d inter;adry*100-(100-adry*100)"

"draw xlab intersection CRH %"
"draw ylab delta area % (dry-wet)"

"printim multimodel_inter_darea.png x2048 y1536"
"c"


"set frame on"
"set xlab on"
"set ylab on"
"set grid on"
"set grads off"
"set parea 1.2 7.5 1.2 7.5"

x=8
y=8

i=2
while (i<=24)
ii=i-1
"set t "i" "i""
"set x 1 2"
"set gxout scatter"
"set vrange 22 42"
"set vrange2 22 30"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d dry;(1-addry-awwet)*100"
if (i=2)
"set frame off"
"set xlab off"
"set ylab off"
"set grid off"
endif

"set line "cc.ii" 1 14"
"draw mark "cm.ii" "x" "y" 0.25"
"set strsiz 0.125"
"set string 1 l 13 0"
"draw string "x+0.2" "y" "model.ii""

y=y-0.35

i=i+1
endwhile

ii=24
"set t 0 1"
"set gxout scatter"
"set vrange 22 42"
"set vrange2 22 30"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d dry;(1-addry-awwet)*100"

"draw xlab delta CRH % (wet-dry)"
"draw ylab transition zone % (dry+wet)"

"printim multimodel_dcrh_tarea.png x2048 y1536"
"c"


"set frame on"
"set xlab on"
"set ylab on"
"set grid on"
"set grads off"
"set parea 1.2 7.5 1.2 7.5"

x=8
y=8

i=2
while (i<=24)
ii=i-1
"set t "i" "i""
"set x 1 2"
"set gxout scatter"
"set vrange 22 42"
"set vrange2 10 20"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d dry;(adry-addry)*100"
if (i=2)
"set frame off"
"set xlab off"
"set ylab off"
"set grid off"
endif

"set line "cc.ii" 1 14"
"draw mark "cm.ii" "x" "y" 0.25"
"set strsiz 0.125"
"set string 1 l 13 0"
"draw string "x+0.2" "y" "model.ii""

y=y-0.35

i=i+1
endwhile

ii=24
"set t 0 1"
"set gxout scatter"
"set vrange 22 42"
"set vrange2 10 20"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d dry;(adry-addry)*100"

"draw xlab delta CRH % (wet-dry)"
"draw ylab transition zone % (dry)"

"printim multimodel_dcrh_dtarea.png x2048 y1536"
"c"

"set frame on"
"set xlab on"
"set ylab on"
"set grid on"
"set grads off"
"set parea 1.2 7.5 1.2 7.5"

x=8
y=8

i=2
while (i<=24)
ii=i-1
"set t "i" "i""
"set x 1 2"
"set gxout scatter"
"set vrange 22 42"
"set vrange2 7 17"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d dry;(1-adry-awwet)*100"
if (i=2)
"set frame off"
"set xlab off"
"set ylab off"
"set grid off"
endif

"set line "cc.ii" 1 14"
"draw mark "cm.ii" "x" "y" 0.25"
"set strsiz 0.125"
"set string 1 l 13 0"
"draw string "x+0.2" "y" "model.ii""

y=y-0.35

i=i+1
endwhile

ii=24
"set t 0 1"
"set gxout scatter"
"set vrange 22 42"
"set vrange2 7 17"
"set digsiz "0.3+ii/80.""
"set ccolor "cc.ii""
"set cthick 13"
"set cmark "cm.ii""
"set line 1 1 13"
"d dry;(1-adry-awwet)*100"

"draw xlab delta CRH % (wet-dry)"
"draw ylab transition zone % (wet)"

"printim multimodel_dcrh_wtarea.png x2048 y1536"
"c"


"quit"
