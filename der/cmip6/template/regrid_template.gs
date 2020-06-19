"reinit"

path="/data/der0318/work/cmip6/cmip_ln/modelname/"
outpath="/data/der0318/work/cmip6/cmip_2p5d/modelname/"

moda.1=31
moda.2=28
moda.3=31
moda.4=30
moda.5=31
moda.6=30
moda.7=31
moda.8=31
moda.9=30
moda.10=31
moda.11=30
moda.12=31

mon.1="JAN"
mon.2="FEB"
mon.3="MAR"
mon.4="APR"
mon.5="MAY"
mon.6="JUN"
mon.7="JUL"
mon.8="AUG"
mon.9="SEP"
mon.10="OCT"
mon.11="NOV"
mon.12="DEC"

z.1="1000"
z.2="850"
z.3="700"
z.4="500"
z.5="250"
z.6="100"

var.1="zg"
var.2="ta"
var.3="hus"
var.4="wap"

"set cachesf 10000000"

i=1
while(i<=4)
"xdfopen "path"modelname_"var.i".ctl"
i=i+1
endwhile
"open grid_2p5d.ctl"

"set lon 0 357.5"
"set lat -88.75 88.75"

yr=1980
while(yr<=2014)

"set undef -999."
"set fwrite "outpath"modelname_"math_format('%4g',yr)"_2p5d.dat"
"set gxout fwrite"

mo=1
while(mo<=12)

da=1
while(da<=moda.mo)

"set time "da""mon.mo""yr""

i=1
while(i<=4)

k=1
while(k<=6)
"set lev "z.k""

"d lterp("var.i"."i",grid.5(z=1,t=1),aave,50)"

k=k+1
endwhile

i=i+1
endwhile


da=da+1
endwhile

out=yr"-"mo
say out

mo=mo+1
endwhile

"disable fwrite"

yr=yr+1
endwhile


"quit"
