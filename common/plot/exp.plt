#==========================================================================
# Sphere 
# gnuplot command file
# Daniel J. Greenhoe
# in gnuplot, enter
#    >load "sphere.plt"
#==========================================================================

 reset
#---------------------------------------
# parameters
#---------------------------------------
 k = 1    # sphere morphing factor (sphere: k=1)
 r = 1    # sphere radius
 
#---------------------------------------
# free variables
#---------------------------------------
 set parametric                        # parametric mode
 set data style lines
 set urange [-pi/2:+pi/2]
 
#---------------------------------------
# rendering
#---------------------------------------
 set hidden3d
 set isosamples 50,40
 set samples 100,100
 set size square
 set palette color
 set palette defined (0 "#101010", 0.40 "#0000ff", 0.7 "#00ff00", 0.8 "#e0e000", 1 "#ff0000")
 #                        grey            blue           green          yellow       red           

#---------------------------------------
# support labeling
#---------------------------------------
#set title                             # 
#set hidden3d
#set contour base
 set key off
#set ticslevel 0
#set xzeroaxis
#set yzeroaxis
#set xtics("-r" -r, "0" 0, "r" r);
#set ytics("-r" -r, "0" 0, "r" r);
#set ztics("-r" -r, "0" 0, "r" r);
 unset border
 unset xtics
 unset ytics
 unset ztics
#set xlabel "x"
#set ylabel "y"
#set zlabel "z"
  
#---------------------------------------
# define functions
#---------------------------------------
 x(t) = r*cos(t)         # x component of ball
 y(t) = r*sin(t)         # y component of ball
 
#---------------------------------------
# plot figure
#---------------------------------------
 plot x(t), y(t) with lines
# splot x(u,v), y(u,v), z(u) with line palette
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "exp.eps"
 replot
 set output
 set terminal pop

