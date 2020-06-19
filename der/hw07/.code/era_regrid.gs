"reinit"

path="/data/dadm1/reanalysis/ECMWF/ITM/daily/"

var.1="Z"
var.2="T"
var.3="Q"
var.4="W"
var.5="SP"

t=1979
while (t<=2017)
i=1
while (i<=5)
"sdfopen "path""var.i"/daily_interim_"var.i"_"t".nc"
i=i+1
endwhile

"open /data/cloud/der/hw07/.code/grid.ctl"

"q file"
temp=sublin(result,5)
tend=subwrd(temp,12)
zend=subwrd(temp,9)

"set lon 0 357.5"
"set lat -88.75 88.75"

"set undef -99999."
"set fwrite /data/cloud/der/hw07/for_you/daily_interim_"t"_2p5d.dat"
"set gxout fwrite"

tt=1
while (tt<=tend)
"set t "tt""

k=1
while (k<=zend)
"set z "k""
"d lterp(z.1,grid.6(z=1,t=1),aave)"
k=k+1
endwhile
k=1
while (k<=zend)
"set z "k""
"d lterp(t.2,grid.6(z=1,t=1),aave)"
k=k+1
endwhile
k=1
while (k<=zend)
"set z "k""
"d lterp(q.3,grid.6(z=1,t=1),aave)"
k=k+1
endwhile
k=1
while (k<=zend)
"set z "k""
"d lterp(w.4,grid.6(z=1,t=1),aave)"
k=k+1
endwhile
"set z 1"
"d lterp(sp.5(z=1),grid.6(z=1,t=1),aave)"

tt=tt+1
endwhile

i=1
while (i<=6)
"close "7-i""
i=i+1
endwhile

"disable fwrite"

t=t+1
endwhile


"quit"
