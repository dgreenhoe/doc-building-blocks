printf(" ============================================ \n");
printf("| Cubic Spline Wavelet; m=3                  |\n");
printf("|                                            |\n");
printf("| Daniel Greenhoe  9013605                   |\n");
printf("| NCTU -- National Chiao-Tung University     |\n");
printf("| using GNU Octave   http://www.octave.com   |\n");
printf("|                                            |\n");
printf(" ============================================ \n");

#----------------------------------------------------------
# Reference:
#     Stephane G. Mallat 
#     A Wavelet Tour of Signal Processing 2nd edition. 
#     Academic Press. 1999 
#     ISBN 0-12-466606-X 
#     http://www.apnet.com/fatbrain/mallat_ch1.pdf
#     page 235
#----------------------------------------------------------

#======================================
# Initialization
#======================================
clear;

#======================================
# Function Declarations
#======================================

#function: unwrap
#-----------------------------------
function y = unwrap(x)
   n = length(x);
   y = zeros(size(x));
   y(1:n/2) = x(n/2+1:n);
   y(n/2+1:n) = x(1:n/2);
#   y = [ x(n/2+1:n); x(1:n/2)];
endfunction

#function: Box window function
#-----------------------------------
function y = chi(t)
   n = length(t);
   for i=1:n
      if     (t(i)< 0)    y(i)=0;
      elseif (t(i)> 1)    y(i)=0;
      else                y(i)=1;
      endif
      endfor
endfunction

#function: Cubic Spline m 
#-----------------------------------
function y = spline(m,t)
   n  = length(t);
   y0 = ones(1,n);
   ym = y0;
   for i=1:m
      ym = conv(y0,ym);
   endfor
   
   odd = rem(m,2);
   if(odd)
      t_hi = 0   + (m+1)/2;
      t_lo = 0   - (m+1)/2;
   else
      t_hi = 1/2 + (m+1)/2;
      t_lo = 1/2 - (m+1)/2;
   endif
   
   ym = ym/max(ym);
   nm = length(ym);
   
   for i=1:n;
      im = nm*(t(i)-t_lo)/(t_hi-t_lo);
      if( im<=0 ) y(i)= 0;
      elseif ( im>nm ) y(i) = 0;
      else y(i) = ym( im );
      endif
   endfor
endfunction

#function: Plot eps file
#-----------------------------------
function eps_plot(filename)
   graw("set terminal postscript eps enhanced color blacktext solid linewidth 3 rounded \n");
   if( length(findstr(filename,"."))==0)
      filename = [filename,".eps"];
   endif
   output_cmd = sprintf("set output '%s' \n",filename);
   printf("   output file name = \"%s\" \n",filename);
   graw(output_cmd); 
   replot
   gset output
   gset terminal windows
endfunction


#----------------------------------------------------------------------------
# GNU-Plot Setup
#----------------------------------------------------------------------------
function gnuplot_setup()
   graw("\n\n");
   graw("#-----------------------------------\n");
   graw("# New plot                          \n");
   graw("#-----------------------------------\n");
   graw("clear\n");
   graw("reset\n");
   graw("unset key \n");
   graw("unset grid \n");
   gset style data lines
   gset style function lines
   gset style line 1 linetype 3 linewidth 1 pointtype 3 pointsize 1
   gset ticslevel 0
   graw("unset label \n")
   graw("unset xzeroaxis \n")
   graw("unset yzeroaxis \n")
   gset view 60,35
   gset xlabel "n"
   gset ylabel "alpha"
   graw("unset hidden3d \n")
   graw("unset contour \n")
   graw("unset border \n")
endfunction

#======================================
# Procedural Blocks
#======================================
gnuplot_setup();


# Parameters
#-----------------------------------
delay_sec=2;
M = 41;
N = 512;

h = zeros(M,1);
h(1:(M-1)/2+1) = [            \
            0.000146098;
           -0.000232304;
           -0.000285414;
            0.000462093;
            0.000559952;
           -0.000927187;
           -0.001103748;
            0.001882120;
            0.002186714;
           -0.003882426;
           -0.004353840;
            0.008201477;
            0.008685294;
           -0.017982291;
           -0.017176331;
            0.042068328;
            0.032080869;
           -0.110036987;
           -0.050201753;
            0.433923147;
            0.766130398;
          ];
          
h((M-1)/2+2:M) = h((M-1)/2:-1:1);
g = h;
for i=1:M  g(i) = (-1)**(i-1)*h(i);  endfor
          
#-----------------------------------
plot_title = sprintf("h(n) for cubic spline m=3")
#-----------------------------------
data = fopen("cf.dat","w");
for i=1:M  fprintf(data,"%10f  %10f \n",i-(M-1)/2-1, h(i));  endfor
fclose(data);

title(plot_title);
graw(sprintf("plot [][-1:1] 'cf.dat' with impulses linestyle 1 \n"));
eps_plot("sp_3_h");
pause(delay_sec);

#-----------------------------------
plot_title = sprintf("g(n) for cubic spline m=3")
#-----------------------------------
data = fopen("cf.dat","w");
for i=1:M  fprintf(data,"%10f  %10f \n",i-(M-1)/2-1, g(i));  endfor
fclose(data);

title(plot_title);
graw(sprintf("plot [][-1:1] 'cf.dat' with impulses linestyle 1 \n"));
eps_plot("sp_3_g");
pause(delay_sec);

#-----------------------------------
plot_title = sprintf("chi function")
#-----------------------------------
t = 2*(2*[0:(N-1)]/(N-1)-ones(1,N));
y = chi(t);

data = fopen("cf.dat","w");
for i=1:N  fprintf(data,"%10f  %10f \n",t(i), y(i));  endfor
fclose(data);

title(plot_title);
gset xlabel "t"
graw(sprintf("plot [][-1:1.5] 'cf.dat' with lines linestyle 1 \n"));
eps_plot("chi");
pause(delay_sec);


for m=0:5
#-----------------------------------
#m=2;
plot_title = sprintf("cubic spline m=%d scaling function",m)
#-----------------------------------
#maxmin = (m+1)/2 + 1;
maxmin = 4;
t = maxmin*(2*[0:(N-1)]/(N-1)-ones(1,N));
y = spline(m,t);

data = fopen("cf.dat","w");
for i=1:N  fprintf(data,"%10f  %10f \n",t(i), y(i));  endfor
fflush(data);
fclose(data);

title(plot_title);
gset xlabel "t"
graw(sprintf("plot [][-0.5:1.5] 'cf.dat' with lines linestyle 1 \n"));
file_name = sprintf("sp_%d_s",m);
eps_plot(file_name);
pause(delay_sec);

endfor
 
#======================================
# End Processing
#======================================



