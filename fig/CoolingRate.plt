set term pdfcairo enhanced dl 2. color font 'Times, 16' size 5, 4

s2d=24*3600
s2h=1*3600
llw=3 
mat='iro'
siz='01'
nam1='500K'
nam2='300K'
nam3='100K'
set output 'CoolingRate_'.mat.siz.'nm.pdf'

f1='../out/Cool-grain-'.mat.siz.'nm-Tmax'.nam1.'.out'
f2='../out/Cool-grain-'.mat.siz.'nm-Tmax'.nam2.'.out'
f3='../out/Cool-grain-'.mat.siz.'nm-Tmax'.nam3.'.out'

set xrange[1e-4:20]
set yrange[10:*]
unset grid 
set mxtics 10
set mytics 10
set logscale x
unset logscale y
#set format '10^{%L}'
set xlabel 'Time after absorption [hours]'

set ylabel 'Grain temperature T [K]'

plot f1 u ($1/s2h):3 w l lc 20 lw llw  t 'T_{ini}='.nam1\
,    f2 u ($1/s2h):3 w l lc 15 dt 2 lw llw  t 'T_{ini}='.nam2\
,    f3 u ($1/s2h):3 every 10 w l lc 22 dt 3 lw llw  t 'T_{ini}='.nam3\

set logscale y
set format '10^{%L}'

set yrange[*:*]
set ylabel 'Cooling rate -dT/dt [K/s]'
plot f1 u ($1/s2h):(-$4) w l lc 20 lw llw t ''\
,    f2 u ($1/s2h):(-$4) w l lc 15 dt 2 lw llw t ''\
,    f3 u ($1/s2h):(-$4) every 10 w l lc 22 dt 3 lw llw t ''\


set ylabel '-T/ (dT/dt) [s]'
plot f1 u ($1/s2h):(-$3/$4) w l lw llw lc 20 t ''\
,    f2 u ($1/s2h):(-$3/$4) w l lw llw lc 15 dt 2 t ''\
,    f3 u ($1/s2h):(-$3/$4) w l lw llw lc 22 dt 3 t ''\
