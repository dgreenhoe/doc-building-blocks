#==========================================================================
# gnuplot command file
# plot pyramid from data file
# Daniel J. Greenhoe
# in gnuplot, enter at the command line
#    > load "<this file>"
#==========================================================================

#---------------------------------------
# parameters
#---------------------------------------
 a = 1
 h = 1/sqrt(2.)

#---------------------------------------
# Set plot environment
#---------------------------------------
 reset
 set noparametric
 set title
 set data style lines
 set xrange [0:a]
 set yrange [0:a]
#set contour base
 set key off
 set isosamples 50,40
 set samples 20,20
 set size square
 set view 60, 30
#set ticslevel 0.5
#set xlabel "x"
#set ylabel "y"
#set zlabel "z"
#set xtics("-a" 0, "-2a/3" 1./6., "0" 0.5, "2a/3" 5./6., "a" a);
#set ytics("-a/sqrt(3)" 0, "0" 1./(2.*sqrt(3)), "2a/sqrt(3)" 1./sqrt(3));
#set ztics("0" 0, "h/2" h/2., "h" h );
 set hidden3d
 unset border
 unset xtics
 unset ytics
 unset ztics

#---------------------------------------
# plot figure to monitor
#---------------------------------------
 set palette color
 set palette defined (0 "#101010", 0.40 "#0000ff", 0.7 "#00ff00", 0.8 "#e0e000", 1 "#ff0000")
 splot "pyr_tri.dat" with line palette

#---------------------------------------
# generate eps file
#---------------------------------------
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "pyr_tri.eps"
 replot
 set output
 set terminal pop

