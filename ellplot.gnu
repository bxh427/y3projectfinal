reset
set terminal x11
set format xy ""
unset key
a = 2
b = 1
frac = b/a
set size ratio frac
set xrange [-0.5-a:a+0.5]
set yrange [-0.5-b:b+0.5]
set object 1 ellipse center 0,0 size 2*a,2*b
plot "gnuEllipse.txt" with lines lc rgb "red" lt 1
