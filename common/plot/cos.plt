#==========================================================================
# cosine functions 
# gnuplot command file
# Daniel J. Greenhoe
# in gnuplot, enter
#    >load "<this file name"
#==========================================================================

 reset
#---------------------------------------
# parameters
#---------------------------------------
 set xrange [-1:+1]
 set yrange [-1.2:+1.2]
 twscale=2./16.;
 textwidth=160.;  # text width in LaTeX document in mm
 lwidth = 2*(5*25.4)/(textwidth*twscale)
 
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
# plot figure
#---------------------------------------
 n=7
#plot cos(n*2*pi*x) linestyle (n==0)? blue : red
 plot (cos(2*pi*x))**n linestyle (n==0)? blue : red
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
#set output "cos_nhz.eps"
 set output "cos_en.eps"
 replot
 set output
 set terminal pop

