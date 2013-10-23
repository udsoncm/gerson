reset
unset key
set pm3d map interpolate 0,0

set xr [-200:200]

set xlabel "x [nm]"
set ylabel "E [meV]"
set title "LDOS along the wire showing the Van Hove singularity (red)"

splot [][0:70] "ldos.dat" u 1:2:3