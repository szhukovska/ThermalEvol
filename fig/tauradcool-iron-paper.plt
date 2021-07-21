set term pdfcairo enhanced dl 1. ps 0.7 color font 'Times, 16' size 5, 4.5

set mxtics 10
set mytics 10
set logscale x
set logscale y
set format '10^{%L}'
set key left bottom reverse L
set ylabel 'Cooling time {/=18  {/Symbol t}_{rad} [s]'
set xlabel 'Vibrational energy E/hc [cm^{-1}]'


set output 'tauradcool-iron-sil.pdf'
set yrange[2e-1:3e7]

m="#0B5fa5"
k="#00C618"

p='../out/'
f0=p.'tauradcool-iro100A.out'
f10=p.'tauradcool-iro50A.out'
f1=p.'tauradcool-iro30A.out'
f21=p.'tauradcool-iro15A.out'
f2=p.'tauradcool-iro10A.out'
ff0=p.'tauradcool-sil50A.out'
ff1=p.'tauradcool-sil30A.out'
ff2=p.'tauradcool-sil20A.out'
ff3=p.'tauradcool-sil10A.out'
m="#0B5fa5"
s='{/Times-Italic {/=18 a}}='
################################

mm=5; ll=7; nn=9
pp=30
El=109690.8e99

plot    f10 u ($1):($1< El ? $2:1/0) t 'iron, '.s.'50 nm' w lp pt ll pi pp lw 4 lc 20 dt 1\
,    f1 u ($1):($1< El ? $2:1/0) t 'iron, '.s.'3 nm' w lp pt mm pi pp  lw 4 lc 7 dt 1\
,    f2 u ($1):($1< El ? $2:1/0) t 'iron, '.s.'1 nm' w lp pt nn pi pp  lw 4 lc 3 dt 1\
,    ff0 u ($1):($1< El ? $2:1/0) t 'silicate, '.s.'50 nm' w lp pt ll pi pp  lw 3.5 lc 20 dt 2\
,    ff1 u ($1):($1< El ? $2:1/0) t 'silicate, '.s.'3 nm' w lp pt mm pi pp  lw 3.5 lc 7 dt 2\
,    ff3 u ($1):($1< El ? $2:1/0) t 'silicate, '.s.'1 nm' w lp pt nn pi pp  lw 3.5 lc 3 dt 2\
