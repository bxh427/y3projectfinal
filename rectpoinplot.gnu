reset
set terminal x11
set format xy "%g"
unset key
set xlabel "theta"
set ylabel "alpha"
set xrange [0:2]
plot "gnupoincareRectangle.txt" lc rgb "red"
unset terminal 
