#==========================================================================
# Sphere 
# gnuplot command file
# Daniel J. Greenhoe
# in gnuplot, enter
#    >load "<this file name>"
#==========================================================================

 reset
#---------------------------------------
# parameters
#---------------------------------------
 a = 1    # x,y,z in [-a,+a]
 
#---------------------------------------
# free variables
#---------------------------------------
 set parametric                        # parametric mode
 set data style lines
 set urange [-pi/2:+pi/2]
 set vrange [-pi:+pi]
 
#---------------------------------------
# rendering
#---------------------------------------
 set hidden3d
 set isosamples 50,50
 set samples 400,400
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
#set xtics("-a" -a, "0" 0, "a" a);
#set ytics("-a" -a, "0" 0, "a" a);
#set ztics("-a" -a, "0" 0, "a" a);
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
 b = 1*sqrt(3)*a
 xr(u,v) = b*cos(u)*cos(v)               
 yr(u,v) = b*cos(u)*sin(v)
 zr(u)   = b*sin(u)

 x(u,v) = (xr(u,v)>=a)? a : (xr(u,v)<=-a)? -a : xr(u,v)
 y(u,v) = (yr(u,v)>=a)? a : (yr(u,v)<=-a)? -a : yr(u,v)
 z(u)   = (zr(u)  >=a)? a : (zr(u)  <=-a)? -a : zr(u)
 
#---------------------------------------
# plot figure
#---------------------------------------
 splot x(u,v), y(u,v), z(u) with line palette
#splot x(u,v), y(u,v), z(u) with points
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "cube2.eps"
 replot
 set output
 set terminal pop

