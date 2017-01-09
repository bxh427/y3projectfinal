set terminal epslatex
set output "StadDivTraj.tex"
set format xy ""
set parametric
unset key
a = 1
b = 2
frac = a/b
set size ratio frac 
set xrange [(-0.5*a)-0.1:b+(0.5*a)+0.1]
set yrange [-0.1:a+0.1]
set arrow 1 from 0,0 to b,0 nohead lw 0.1
set arrow 2 from 0,a to b,a nohead lw 0.1
plot [-pi/2:pi/2] -0.5*a*cos(t),(0.5*a)+0.5*a*sin(t) lc rgb "black" lt 1,\
b+0.5*a*cos(t),0.5*a+0.5*a*sin(t) lc rgb "black" lt 1,\
"gnuStadium.txt" with lines lc rgb "red" lt 1
unset output
unset terminal

set terminal x11
set format xy ""
set parametric
a = 1
b = 2
frac = a/b
set size ratio frac
set xrange [(-0.5*a)-0.1:b+(0.5*a)+0.1]
set yrange [-0.1:a+0.1]
set arrow 1 from 0,0 to b,0 nohead lw 0.1
set arrow 2 from 0,a to b,a nohead lw 0.1
plot [-pi/2:pi/2] -0.5*a*cos(t),(0.5*a)+0.5*a*sin(t) lc rgb "black" lt 1,\
b+0.5*a*cos(t),0.5*a+0.5*a*sin(t) lc rgb "black" lt 1,\
"gnuStadium.txt" with lines lc rgb "red" lt 1
unset terminal
