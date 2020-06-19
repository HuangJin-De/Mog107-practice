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

"open isen.ctl"

"set mproj off"

"set z 1"
"set x 1 100"
"set grads off"
"set parea 1.2 10 1.2 7.5"
"set xlabs 0|10|20|30|40|50|60|70|80|90|100"
"set vrange 0 0.025"
"set ylint 0.005"

"set cthick 14"
"set cmark 0"
"d f(t=1)"

"set frame off"
"set xlab off"
"set ylab off"
"set grid off"

"set cthick 13"
"set ccolor 2"
"set cmark 0"
"d f(t=2)"

"set cthick 13"
"set ccolor 4"
"set cmark 0"
"d f(t=3)"

"draw title PDF"
"draw xlab column relative humidity (%)"
"printim oni.png x2048 y1536"
"c"

"quit"
