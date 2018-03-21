#==========================================================================
# Sphere 
# gnuplot command file
# Daniel J. Greenhoe
# in gnuplot, enter
#    >load "<this file name>"
# reference:
# Alfred Gray, 
# Modern Differential Geometry of Curves and Surfaces... 2nd edition
# 1998, page 327
#==========================================================================

#---------------------------------------
# parameters
#---------------------------------------
 a = 4
 
#---------------------------------------
# Set plot environment
#---------------------------------------
 reset
 set parametric
 set title
 set data style lines
 set hidden3d
#set contour base
 set key off
 set isosamples 60,40
 set samples 100,100
 set view 30,30
 set size square
#set size noratio
 set ticslevel 0
#set xzeroaxis
#set yzeroaxis
 
#set xtics("-r" -r, "0" 0, "r" r);
#set ytics("-r" -r, "0" 0, "r" r);
#set ztics("-r" -r, "0" 0, "r" r);
 unset border
 unset xtics
 unset ytics
 unset ztics
  
 set urange [0:2*pi]
 set vrange [0:2*pi]
#set xrange [-r:r]
#set yrange [-r:r]
#set xlabel "x"
#set ylabel "y"
#set zlabel "z"

#---------------------------------------
# plot figure to monitor
#---------------------------------------
 set palette color
 set palette defined (-10 "#000000", -1 "#101010", -0.5 "#0000ff", 0 "#00ff00", 0.5 "#e0e000", 1 "#ff0000", 11 "#ff0000")
 splot (a+cos(u/2)*sin(v)-sin(u/2)*sin(2*v))*cos(u), \
       (a+cos(u/2)*sin(v)-sin(u/2)*sin(2*v))*sin(u), \
       sin(u/2)*sin(v)+cos(u/2)*sin(2*v) \
       with line palette
      

#---------------------------------------
# plot to eps
#---------------------------------------
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "klein_8i.eps"
 replot
 set output
 set terminal pop

