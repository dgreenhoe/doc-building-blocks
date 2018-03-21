#==========================================================================
# Mandelbrot 
# gnuplot command file
# daniel greenhoe
# in gnuplot, enter
#    >load "mandel.plt"
#==========================================================================
#---------------------------------------
# Set plot environment
#---------------------------------------
 reset
 set title
 set hidden3d
 set nokey
 set samples 100,100
 set view ,,0.7,1.7
 set ticslevel 0.25

# mandelbrot demo
 set nopar
 set mapp cart
#set size square
 set view 60,30,1,1
 set auto
 set title "" 0,0
 set isosamples 60
 set hidden3d
# unset border
# unset xtics
# unset ytics
# unset ztics

#---------------------------------------
# define functions
#---------------------------------------
 compl(a,b)=a*{1,0}+b*{0,1}
 mand(z,a,n) = n<=0 || abs(z)>100 ? 1:mand(z*z+a,a,n-1)+1

#---------------------------------------
# plot figure to monitor
#---------------------------------------
 set palette color
 set palette defined (0 "#101010", 0.10 "#0000ff", 0.4 "#00ff00", 0.6 "#e0e000", 1 "#ff0000")
 splot [-2:1][-1.5:1.5] mand({0,0},compl(x,y),30) with line palette
# TANAKA Masaki  (Tokyo Institute of technology)
#                 masaki@isea.is.titech.ac.jp

#---------------------------------------
# generate eps file
#---------------------------------------
 set terminal push
 set terminal postscript eps noenhanced color solid
 set output "mandel.eps"
 replot
 set output
 set terminal pop

