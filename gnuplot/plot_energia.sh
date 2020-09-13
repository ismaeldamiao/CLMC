#!/usr/bin/gnuplot -c

set term png size 1920, 1080 font ",18"
set output sprintf("%s.png", ARG1)

set key off

set palette rgbformulae 23,28,3
set ticslevel 0
set surface
set hidden3d
set grid

set xlabel "t"
set ylabel "Massa n"
set zlabel "Energia" offset -1,1,0 rotate by 90
set view 50.0, 75.0

splot ARG1 with pm3d palette
