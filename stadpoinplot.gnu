reset
set terminal x11
unset key
set format xy "%g"
set xlabel "theta"
set ylabel "alpha"
set xrange [0:2]
set yrange [-0.5:0.5]
set xtics 0,0.5,2
set ytics -0.5,0.1,0.5
plot 'gnupoincareStadium.txt' with dots lc rgb "red"
unset multiplot
unset output
