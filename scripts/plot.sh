#!/usr/bin/bash

source info.txt

if [ -e "energia_*.dat.dat" ]; then
   ENERGIA="energia_*.dat.dat"
else
   ENERGIA="${1}"
fi

if [ -e "dispersao_*.dat.dat" ]; then
   DISPERSAO="dispersao_*.dat.dat"
else
   DISPERSAO="${2}"
fi

gnuplot <<EOF
set term png size 1080, 1080 font "serif,26" enhanced

set key off
set grid
set title "N=${N}, α=${ALPHA}, v_0=${V0}\nη_2=${ETA2},η_3=${ETA3},η_4=${ETA4}"

# ###
# Sigma
# ###

set output "sigma.png"
set xlabel "t"
set ylabel "σ(t)"
plot "${DISPERSAO}" using 1:2 w l lt rgb "black"

# ###
# Z
# ###

set output "z.png"
set xlabel "t"
set ylabel "Z(t)"
plot "${DISPERSAO}" using 1:3 w l lt rgb "black"

# ###
# Sigma
# ###

set term png size 1920, 1080 font "serif,26" enhanced

set palette rgbformulae 23,28,3
set ticslevel 0
set surface
set hidden3d

set output "energia.png"
set xlabel "t"
set ylabel "Massa n"
set zlabel "Energia" offset -1,1,0 rotate by 90
set view 50.0, 75.0

splot "${ENERGIA}" with pm3d palette
exit
EOF

exit 0