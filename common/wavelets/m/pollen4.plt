#==========================================================================
# Sphere 
# gnuplot command file
# Daniel J. Greenhoe
# in gnuplot, enter
#    >load "<this filename>"
#==========================================================================

 reset

#---------------------------------------
# parameters
#---------------------------------------
 k = 1    # sphere morphing factor (sphere: k=1)
 r = 1    # sphere radius
 a = .4    # x,y,z in [-a,+a]
 
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
set hidden3d
#set contour base
 set key off
#set ticslevel 0
#set xzeroaxis
#set yzeroaxis
 unset border
 set xtics
 set ytics
 unset ztics
set xlabel "alpha"
set ylabel "t"
#set zlabel "psi(t)"
set view 45,109
  

#---------------------------------------
# plot figure
#---------------------------------------
splot 'p4psi_256x256.dat' every 5:8:0:0:255:255 with line
#splot 'p4psi_256x256.dat' every 5:8:0:0:255:255 with line palette
#splot 'p4psi_256x256.dat' every 5:8:38:0:170:256 with line palette

 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "p4psi_djg.eps"
 replot
 set output
 set terminal pop

