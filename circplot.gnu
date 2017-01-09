reset
set terminal x11
set format xy ""
set size square	
set parametric
unset key
set xrange [-1.1:1.1]
set yrange [-1.1:1.1]
plot [0:2*pi] cos(t),sin(t) lc rgb "black",\
"gnuCircle.txt" with lines lc rgb "red" lt 1
unset multiplot
