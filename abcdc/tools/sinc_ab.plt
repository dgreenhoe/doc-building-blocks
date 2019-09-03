# ============================================ 
#| Daniel Greenhoe                            |
# ============================================ 

#----------------------------------------------------------------------------
# GNU-Plot Setup
#----------------------------------------------------------------------------
clear 
reset 
unset key 
unset grid
set   terminal windows 
#set   style data lines 1 linewidth 3
set   style data lines 
set   xzeroaxis
set   yzeroaxis
set   ticslevel 0
unset label
set yzeroaxis
unset border

unset title
set   xlabel "f"
set   ylabel
set   xrange [-2:2]
set   yrange [-0.5:2.2]
set key top right box
unset key
set xtics("-2/b" -2, "-3/2b" -1.5, "-1/b" -1, "-1/2b" -0.5, "0" 0, "1/2b" 0.5, "1/b" 1, "3/2b" 1.5, "2/b" 2);
set ytics("0" 0, "2ab" 2)
a=1;
b=1;

F(x) = 2*a*b*sin(2*pi*x*b)/(2*pi*x*b)
plot  F(x)

set terminal postscript portrait enhanced color
set output "sinc_ab.eps"
replot

set terminal latex
set output "sinc_ab.lpl"
replot

set output
set terminal windows





