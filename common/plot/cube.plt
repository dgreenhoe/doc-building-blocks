#==========================================================================
# gnuplot command file
# plot cube from cube.dat data file
# Daniel Greenhoe
# in gnuplot, type "load cube.plt"
#==========================================================================
#---------------------------------------
# Set plot environment
#---------------------------------------
 reset
 set title
 set hidden3d
#set contour base
 set key off
 set noparametric
 set hidden3d
 set samples 100,100
 set view 60,30
# set view 45,315
 set size square
 set isosamples 20,20
 set data style lines
 set xrange [0:1]
 set yrange [0:1]
#set ticslevel 0
#set xtics("-a" 0, "0" 0.5, "a" 1);
#set ytics("-a" 0, "0" 0.5, "a" 1);
#set ztics("-a" 0, "0" 0.5, "a" 1);
#set xlabel "x"
#set ylabel "y"
#set zlabel "z"
 unset border
 unset xtics
 unset ytics
 unset ztics

#---------------------------------------
# plot figure to monitor
#---------------------------------------
 set palette color
 set palette defined (0 "#101010", 0.15 "#0000ff", 0.7 "#00ff00", 0.8 "#e0e000", 1 "#ff0000")
 splot  "cube.dat" with line palette

#---------------------------------------
# generate eps file
#---------------------------------------
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "cubedjg_cab.eps"
 replot
 set output
 set terminal pop

