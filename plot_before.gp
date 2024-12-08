set terminal png
set output "before.png"

set title "before positions"
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set xrange [0:10]
set yrange [0:10]
set zrange [0:10]

set border linewidth 1
set tics out

splot "Before.dat" with points pointtype 7 pointsize 0.5