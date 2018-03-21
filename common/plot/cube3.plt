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
 set urange [0:4]
 set vrange [0:2]
 
#---------------------------------------
# rendering
#---------------------------------------
 set hidden3d
 set isosamples 20,40
 set samples 1000,1000
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
set xzeroaxis
set yzeroaxis
 set border
# unset xtics
# unset ytics
# unset ztics
#set xtics("-a" -a, "0" 0, "a" a);
#set ytics("-a" -a, "0" 0, "a" a);
#set ztics("-a" -a, "0" 0, "a" a);
set xlabel "x"
set ylabel "y"
set zlabel "z"
  
#---------------------------------------
# define functions
#---------------------------------------
 x(u,v) = (u>=0 && u<=1)? (v>=1 && v<=2)? a : 2*a*(v-0.5) :  \
          (u>=2 && u<=3)? (v>=0 && v<=1)? -a : 2*a*(v-1.5) :  \
          (u>=1 && u<2)? 2*a*(1.5-u) :                     \
          2*a*(u-3.5)

 y(u,v) = (u>=1 && u<=2)? (v>=1 && v<=2)? +a : 2*a*(v-0.5) :   \
          (u>=3 && u<=4)? (v>=0 && v<=1)? -a : 2*a*(v-1.5) :   \
          (u>=0 && u<=1)? 2*a*(u-0.5) : 2*a*(2.5-u) 

 z(u,v) = (u>=0 && u<=2)? (v>=0 && v<=1)? -a : 2*a*(v-1.5) :  \
          (u>=2 && u<=3)? (v>=1 && v<=2)? +a : 2*a*(v-0.5) :  \
          (v>=1 && v<=2)? +a : 2*a*(v-0.5) 
 
#---------------------------------------
# plot figure
#---------------------------------------
 splot x(u,v), y(u,v), z(u,v) with line palette
#splot x(u,v), y(u,v), z(u,v) with points
#splot u,v, x(u,v)
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "cube2.eps"
 replot
 set output
 set terminal pop

