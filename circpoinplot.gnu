reset
set terminal x11
set format xy "%g"
unset key
set xlabel "theta"
set ylabel "alpha"
set xrange [0:2]
set xtics 0,0.5,2
plot 'gnupoincareCircle.txt' lc rgb "red" with points
unset multiplot
unset terminal
