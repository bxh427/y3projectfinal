reset
set terminal x11
set format xy "%g"
set xlabel "s"
set ylabel "p"
set yrange [-1.2:1.2]
unset key
plot 'gnupoincareEllipse.txt' with dots lc rgb "red"
