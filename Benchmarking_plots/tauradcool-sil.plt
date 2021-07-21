set term pdfcairo enhanced dl 1. color font 'Arial, 16' size 5, 5

set yrange[1e-4:1e7]
#set yrange[1:4e2]
set grid 
set mxtics 10
set mytics 10
set logscale x
set logscale y
set format '10^{%L}'
set key left bottom
set ylabel '{/Symbol t}_{rad}, s'
set xlabel 'E/hc [cm^{-1}]'


p='../out/'

set yrange[2e-3:1e6]
set xrange[10:2e5]
f1=p.'tauradcool-sil30A.out'
f2=p.'tauradcool-sil20A.out'
f3=p.'tauradcool-sil10A.out'

set output 'tauradcool-sil.pdf'

plot  f1 t '30A' w l lw 3\
,    f2 t '20A' w l lw 3\
,    f3 t '10A' w l lw 3\

#,    f4 t '6A' w l lw 3 \ no data for Qabs
#,    f2 t '3A' w l lw 3 \
