#==========================================================================
# gnuplot command file
# Daniel J. Greenhoe
# in gnuplot, enter
#    >load "<this file name"
#==========================================================================

 reset
#---------------------------------------
# parameters
#---------------------------------------
 set xrange [-3:+3]
 set yrange [-3:+3]
 twscale=2./16.;
 textwidth=160.;  # text width in LaTeX document in mm
 lwidth = 2*(5*25.4)/(textwidth*twscale)
 Exx = 1
 Eyy = 1
 Exy = 0.80
 Exy2 = 0.95
 Ex  = 0
 Ey  = 0
  
#---------------------------------------
# rendering
#---------------------------------------
 set samples 32
 set isosamples 32
 set size square
 set view 73,119
 set data style lines
 set style line 1 linetype  1 linewidth lwidth pointtype 13 pointsize 1
 set style line 2 linetype 10 linewidth lwidth pointtype 13 pointsize 1
 set style line 3 linetype  3 linewidth lwidth pointtype 13 pointsize 1
 red=1
 green=2
 blue=3
 set contour surface
 set cntrparam bspline
 set cntrparam levels auto 8
#unset clabel
 set palette color
 set palette defined (0 "#c0c0ff", 0.20 "#0000ff", 0.5 "#00ff00", 0.8 "#e0e000", 1.4 "#ff0000")
 #                        grey            blue           green          yellow       red           

#---------------------------------------
# support labeling
#---------------------------------------
 set key off
#set ticslevel 0
#set xzeroaxis
#set yzeroaxis
#set xtics("-r" -r, "0" 0, "r" r);
#set ytics("-r" -r, "0" 0, "r" r);
#set ztics("-r" -r, "0" 0, "r" r);
#set xtics("Ex" 0);
#set ytics("Ey" 0);
 unset border
 unset xtics
 unset ytics
 unset ztics
#set xlabel "x"
#set ylabel "y"

#---------------------------------------
# functions
#---------------------------------------
 detM = Exx*Eyy-Exy*Exy;
 p(x,y) = (1/(2*pi)*sqrt(detM))*exp(-((x-Ex)**2*Eyy -2*(x-Ex)*(y-Ey)*Exy*Exy + (y-Ey)**2*Exx)/(2*detM))
 detM2 = Exx*Eyy-Exy2*Exy2;
 p2(x,y) = (1/(2*pi)*sqrt(detM2))*exp(-((x-Ex)**2*Eyy -2*(x-Ex)*(y-Ey)*Exy2*Exy2 + (y-Ey)**2*Exx)/(2*detM2))
 q(x,y) = p(x,y) - p2(x,y)

#---------------------------------------
# plot figure
#---------------------------------------
# splot p(x,y)  with line palette
# splot p(x,y)  with line green
 splot q(x,y)  with line green
 set terminal push
 set terminal postscript eps noenhanced color solid linewidth 2
 set output "normxy_80-95.eps"
 replot
 set output
 set terminal pop

