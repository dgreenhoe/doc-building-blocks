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
 a = .4    # x,y,z in [-a,+a]
 
#---------------------------------------
# free variables
#---------------------------------------
 set parametric                        # parametric mode
 set data style lines
 set urange [-pi/2:+pi/2]
 set vrange [0:2*pi]
 
#---------------------------------------
# rendering
#---------------------------------------
 set hidden3d
 set isosamples 50,40
 set samples 100,100
# set size square
 set size noratio
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
 xs(u,v) = cos(u)*cos(v)         # x component of sphere
 ys(u,v) = cos(u)*sin(v)         # y component of sphere
 zs(u)   = 10*sin(u)+16                # z component of sphere

 b = 1*sqrt(3)*a
 xr(u,v) = b*cos(u)*cos(v)               
 yr(u,v) = b*cos(u)*sin(v)
 zr(u)   = b*sin(u)

 xc(u,v) = (xr(u,v)>=a)? a : (xr(u,v)<=-a)? -a : xr(u,v)
 yc(u,v) = (yr(u,v)>=a)? a : (yr(u,v)<=-a)? -a : yr(u,v)
 zc(u)   = (zr(u)  >=a)? a : (zr(u)  <=-a)? -a : zr(u)
 

#---------------------------------------
# plot figure
#---------------------------------------
 splot xs(u,v), ys(u,v), zs(u) with line palette, \
       xc(u,v), yc(u,v), zc(u) with line palette
#      "cube.dat" with line palette



 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "bookcvr.eps"
 replot
 set output
 set terminal pop

