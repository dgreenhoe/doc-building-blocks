#==========================================================================
# gnuplot command file
# plot pyramid from data file
# daniel greenhoe
# in gnuplot, enter at the command line
#    > load "<this file>"
#==========================================================================
reset
set title
set key off
set hidden3d
#set contour base
set nokey
set size square
set noparametric
set hidden3d
set samples 100,100
#set view ,,0.7,1.3
set view 102,40
#set ticslevel 0.5
set ticslevel 0
set isosamples 20,20
set data style lines
a = 1
h = 1/sqrt(2)
set xrange [0:a]
set yrange [0:a]
#set xtics("-a" 0, "0" 0.5, "a" a);
#set ytics("-a" 0, "0" 0.5, "a" a);
#set ztics("0" 0, "h/2" h/2., "h" h );
#set xlabel "x"
#set ylabel "y"
#set zlabel "z"
 unset border
 unset xtics
 unset ytics
 unset ztics
 set palette color
 set palette defined (0 "#101010", 0.15 "#0000ff", 0.7 "#00ff00", 0.8 "#e0e000", 1.5 "#ff0000")
splot "pyr_sq.dat" with line palette

set terminal postscript eps noenhanced color solid linewidth 2
set output "pyrsqdjg_uncropped.eps"
replot
set output
set terminal windows

