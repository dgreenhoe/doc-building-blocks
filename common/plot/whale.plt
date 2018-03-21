#==========================================================================
# Whale
# gnuplot command file
# Daniel J. Greenhoe
# in gnuplot, enter
#    >load "sphere.plt"
#==========================================================================
#---------------------------------------
# Set plot environment
#---------------------------------------
 reset
 set hidden3d
 set key off
#set ticslevel 0.25
#set contour base
 unset border
 unset xtics
 unset ytics
 unset ztics

#---------------------------------------
# plot figure to monitor
#---------------------------------------
 set palette color
 set palette defined (0 "#101010", 0.40 "#0000ff", 0.7 "#00ff00", 0.8 "#e0e000", 1 "#ff0000")
 splot "whale.dat" with line palette


#---------------------------------------
# generate eps file
#---------------------------------------
 set terminal push
set terminal postscript eps noenhanced color solid
set output "whale.eps"
replot
set output
set terminal pop

