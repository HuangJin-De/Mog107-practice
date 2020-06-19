"reinit"

path="/data/dadm1/obs/OLR_NOAA/"

yr=1979
while (yr<=2017)

"sdfopen "path"olr-daily_v01r02_"yr".nc"
"open /data/cloud/der/hw07/.code/grid.ctl"

"set undef -99999."
"set fwrite /data/cloud/der/hw07/for_you/daily_olr_"yr"_2p5d.dat"
"set gxout fwrite"

"q file"
temp=sublin(result,5)
tend=subwrd(temp,12)

"set lon 0 357.5"
"set lat -88.75 88.75"

t=1
while (t<=tend)
"set t "t""

"d lterp(olr,grid.2(z=1,t=1),aave)"

t=t+1
endwhile

"disable fwrite"

"close 2"
"close 1"
say yr
yr=yr+1
endwhile

"quit"
