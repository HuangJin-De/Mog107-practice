"reinit"
"set display color white"
"c"

"set lwid 13 6"
"set lwid 14 8.0"
"set lwid 15 12.0"
"set annot 1 14"
"set strsiz 0.16"
"set xlopts 1 14 0.16"
"set ylopts 1 14 0.16"
"set clopts 1 13 0.08"
"set rgb 200 100 100 100 220"
"set grid on 3 200 4"
"set map 1 1 12"

run="interim"
"open cloud.ctl"
"open isen.ctl"

"set x 1 100"

"set grads off"
"set parea 1.2 10 1 2.6"
"set xlabs 0|10|20|30|40|50|60|70|80|90|100"
"set vrange 0 0.07"
"set ylint 0.02"
"set cthick 14"
"set ccolor 1"
"set cmark 0"
"set z 1"
"d cs.2(z=1)/sum(cs.2(z=1),x=1,x=100)"
"draw ylab PDF"
"draw xlab column relative humidity (%)"
"set lev 0 20000"
"set parea 1.2 10 2.6 7.6"
"set ylint 2000"
"set xlabs 0|10|20|30|40|50|60|70|80|90|100"
"color -7 -2 0.5  -kind white->p31->p32->p33->p34->p35->p36->p37->p38->p39->p310->p311"
"set gxout grfill"
"set xlabs ||||||||||"
"d log10(cst/sum(sum(cst,x=1,x=100),z=1,z=20))"
"xcbar 10.1 10.2 2.8 7.4 -fs 2 -fw 0.14 -fh 0.14 -ft 14"
"draw title Covective cloud top height"
"draw ylab height (m)"
"printim cst_"run".png x2048 y1536"
"c"

"set grads off"
"set parea 1.2 10 1 2.6"
"set xlabs 0|10|20|30|40|50|60|70|80|90|100"
"set vrange 0 0.07"
"set ylint 0.02"
"set cthick 14"
"set ccolor 1"
"set cmark 0"
"set z 1"
"d cs.2(z=1)/sum(cs.2(z=1),x=1,x=100)"
"draw ylab PDF"
"draw xlab column relative humidity (%)"
"set lev 0 20000"
"set parea 1.2 10 2.6 7.6"
"set ylint 2000"
"set xlabs 0|10|20|30|40|50|60|70|80|90|100"
"color -7 -2 0.5  -kind white->p31->p32->p33->p34->p35->p36->p37->p38->p39->p310->p311"
"set gxout grfill"
"set xlabs ||||||||||"
"d log10(csb/sum(sum(csb,x=1,x=100),z=1,z=20))"
"xcbar 10.1 10.2 2.8 7.4 -fs 2 -fw 0.14 -fh 0.14 -ft 14"
"draw title Covective cloud base"
"draw ylab height (m)"
"printim csb_"run".png x2048 y1536"
"c"

"set grads off"
"set parea 1.2 10 1 2.6"
"set xlabs 0|10|20|30|40|50|60|70|80|90|100"
"set vrange 0 0.03"
"set ylint 0.02"
"set cthick 14"
"set ccolor 1"
"set cmark 0"
"set z 1"
"d ncs.2(z=1)/sum(ncs.2(z=1),x=1,x=100)"
"draw ylab PDF"
"draw xlab column relative humidity (%)"
"set lev 0 20000"
"set parea 1.2 10 2.6 7.6"
"set ylint 2000"
"set xlabs 0|10|20|30|40|50|60|70|80|90|100"
"color -7 -2 0.5  -kind white->p31->p32->p33->p34->p35->p36->p37->p38->p39->p310->p311"
"set gxout grfill"
"set xlabs ||||||||||"
"d log10(ncst/sum(sum(ncst,x=1,x=100),z=1,z=20))"
"xcbar 10.1 10.2 2.8 7.4 -fs 2 -fw 0.14 -fh 0.14 -ft 14"
"draw title Non-convective cloud top height"
"draw ylab height (m)"
"printim ncst_"run".png x2048 y1536"
"c"

"set grads off"
"set parea 1.2 10 1 2.6"
"set xlabs 0|10|20|30|40|50|60|70|80|90|100"
"set vrange 0 0.03"
"set ylint 0.02"
"set cthick 14"
"set ccolor 1"
"set cmark 0"
"set z 1"
"d ncs.2(z=1)/sum(ncs.2(z=1),x=1,x=100)"
"draw ylab PDF"
"draw xlab column relative humidity (%)"
"set lev 0 20000"
"set parea 1.2 10 2.6 7.6"
"set ylint 2000"
"set xlabs 0|10|20|30|40|50|60|70|80|90|100"
"color -7 -2 0.5  -kind white->p31->p32->p33->p34->p35->p36->p37->p38->p39->p310->p311"
"set gxout grfill"
"set xlabs ||||||||||"
"d log10(ncsb/sum(sum(ncsb,x=1,x=100),z=1,z=20))"
"xcbar 10.1 10.2 2.8 7.4 -fs 2 -fw 0.14 -fh 0.14 -ft 14"
"draw title Non-convective cloud base"
"draw ylab height (m)"
"printim ncsb_"run".png x2048 y1536"
"c"


"quit"
