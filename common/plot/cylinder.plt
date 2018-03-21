#==========================================================================
# Cylinder 
# gnuplot command file
# daniel greenhoe
# in gnuplot:
#    > load "cylinder.plt"
#==========================================================================

#---------------------------------------
# parameters
#---------------------------------------
h = 1
r = 1

#---------------------------------------
# Set plot environment
#---------------------------------------
 reset
 set key off
 set title
 set hidden3d
 set view 60,30
 set size square
#set contour base
 set parametric
 set samples 10,10
 set isosamples 40,40
 set data style lines
 set dummy theta,v
 set ticslevel 0
#set xtics("-r" -r, "0" 0, "r" r);
#set ytics("-r" -r, "0" 0, "r" r);
#set ztics("0" 0, "h/2" h/2., "h" h);
 set urange [0:2*pi]
 set vrange [0:h]
# set xrange [-r:r]
# set yrange [-r:r]
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
 splot r*v*cos(theta),r*v*sin(theta),1 with line palette, \
       r*cos(theta),r*sin(theta),v  with line palette

#---------------------------------------
# generate eps file
#---------------------------------------
 set terminal push
set terminal postscript eps noenhanced color solid linewidth 2
set output "cylinderdjg_uncropped.eps"
replot
set output
set terminal pop

