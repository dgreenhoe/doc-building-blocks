#==========================================================================
# Cone 
# gnuplot command file
# Daniel J. Greenhoe
# in gnuplot, enter
#    >load "<this file name>"
#==========================================================================
 
#---------------------------------------
# parameters
#---------------------------------------
h = 1
R = 1

#---------------------------------------
# Set plot environment
#---------------------------------------
 reset
 set parametric
 set title
 set data style lines
 set urange [0:2*pi]
 set vrange [0:1]
#set contour base
 set key off
 set hidden3d
 set isosamples 30,20
 set samples 20,20
 set size square
 #set view 74,30
set view 110,0
#set ticslevel 0
#set xzeroaxis
#set yzeroaxis
#set xlabel "x"
#set ylabel "y"
#set zlabel "z"
#set xrange [-R*1.000001:R]
#set xtics("-R" -R, "-R/2" -R/2., "0" 0, "R/2" R/2., "R" R)
#set ytics("-R" -R, "-R/2" -R/2., "0" 0, "R/2" R/2., "R" R)
#set xtics("-R" -R,  "0" 0,  "R" R)
#set ytics("-R" -R,  "0" 0,  "R" R)
#set ztics("0" 0, "h/2" h/2., "h" h)
 unset border
 unset xtics
 unset ytics
 unset ztics

#---------------------------------------
# plot figure to monitor
#---------------------------------------
 set palette color
# set palette defined (0 "#101010", 0.40 "#0000ff", 0.7 "#00ff00", 0.8 "#e0e000", 1 "#ff0000")
 set palette defined (0 "#101010", 0.15 "#0000ff", 0.7 "#00ff00", 0.8 "#e0e000", 1 "#ff0000")
splot R*(1-v)*cos(u), R*(1-v)*sin(u), h*v with line palette
# splot R*(1-v)*cos(u), R*(1-v)*sin(u), h*v with line palette, \
#       R*v*cos(u), R*v*sin(u), h*(v-1) with line palette

#---------------------------------------
# generate eps file
#---------------------------------------
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "conedjg_uncropped.eps"
 replot
 set output
 set terminal pop

 