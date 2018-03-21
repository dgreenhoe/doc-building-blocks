#==========================================================================
# Sphere 
# gnuplot command file
# Daniel J. Greenhoe
# in gnuplot, enter
#    >load "torus.plt"
#==========================================================================

#---------------------------------------
# parameters
#---------------------------------------
 R = 4    # outer radius
 r = 1  # inner radius
 
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
 set hidden3d
 set isosamples 50,40
 set samples 100,100
 set view 30,30
#set size square
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
 set palette defined (-5 "#000000", -1 "#101010", -0.5 "#0000ff", 0 "#00ff00", 0.5 "#e0e000", 1 "#ff0000", 5 "#ff0000")
 splot (R+r*cos(v))*cos(u), (R+r*cos(v))*sin(u), r*sin(v) with line palette

#---------------------------------------
# generate eps file
#---------------------------------------
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "torus.eps"
 replot
 set output
 set terminal pop

