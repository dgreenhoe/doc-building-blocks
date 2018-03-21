#==========================================================================
# Chebychev functions 
# gnuplot command file
# Daniel J. Greenhoe
# in gnuplot, enter
#    >load "<this file name"
#==========================================================================

 reset
#---------------------------------------
# parameters
#---------------------------------------
 twscale=2./16.;
 textwidth=160.;  # text width in LaTeX document in mm
 lwidth = 3*(5*25.4)/(textwidth*twscale)

#---------------------------------------
# rendering
#---------------------------------------
 set samples 1024,1024
 set size square
 set data style lines
 set style line 1 linetype  1 linewidth lwidth pointtype 13 pointsize 1
 set style line 2 linetype 10 linewidth lwidth pointtype 13 pointsize 1
 set style line 3 linetype  3 linewidth lwidth pointtype 13 pointsize 1
 red=1
 green=2
 blue=3

#---------------------------------------
# support labeling
#---------------------------------------
 set key off
#set ticslevel 0
 set xzeroaxis
 set yzeroaxis
#set xtics("-r" -r, "0" 0, "r" r);
#set ytics("-r" -r, "0" 0, "r" r);
#set ztics("-r" -r, "0" 0, "r" r);
 unset border
 unset xtics
 unset ytics
#set xlabel "x"
#set ylabel "y"
 
#---------------------------------------
# Chebyshev polynomial Tn(x)
#   x = cos(theta)
#   Tn(x) = cos(nx)
#---------------------------------------
 set xrange [-1:+1]
 set yrange [-1.2:+1.2]
 theta(x) = acos(x);
 Tn(x,n) = cos(n*theta(x))

#---------------------------------------
# plot figure
#---------------------------------------
 n=10
 plot Tn(x,n) linestyle (n==0)? blue : red
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "Tn.eps"
 replot
 set output
 set terminal pop

