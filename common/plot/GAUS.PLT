#==========================================================================
# Gaussian 
# gnuplot command file
# daniel greenhoe
# in gnuplot, enter
#    >load "<this file>"
#==========================================================================
reset
set title
set data style lines
a = -4
set xrange [ -a:a ]
set yrange [ -a:a ]
set zrange [ : 1./(2.*pi) ]
set hidden3d
set contour base
set nokey
set hidden3d
set isosamples 50,40
set samples 100,100
set view ,,0.7,1.7
set ticslevel 0.50
set xzeroaxis
set yzeroaxis
set xlabel "x"
set ylabel "y"
set zlabel "z"
set xtics( "-4" -4, "-2" -2, "0" 0, "2" 2, "4" 4)
set ytics( "-4" -4, "-2" -2, "0" 0, "2" 2, "4" 4)
set ztics( "0" 0, "1/(4pi)" 1./(4.*pi), "1/(2pi)" 1./(2.*pi) )
f(x,y) = exp(-(x*x+y*y)/2.)/(2.*pi)
splot f(x,y) with lines palette

set terminal postscript eps noenhanced color solid linewidth 2
set output "gaus.eps"
replot
set output
set terminal windows

