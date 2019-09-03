printf(" ============================================ \n");
printf("| Daniel Greenhoe                            |\n");
printf("| National Chiao-Tung University             |\n");
printf("| http://www.nctu.edu.tw                     |\n");
printf(" ============================================ \n");

#======================================
# Initialization
#======================================
clear;

#======================================
# Function Declarations
#======================================

#----------------------------------------------------------------------------
#function: plot sinc^m function
#  N:   number of points
#  m:   power
#----------------------------------------------------------------------------
function plot_sinc(N,m)
   T = 1;
   x = [-3.5*T: 6/N : +3.5*T];
   y = (2/T)*sin(2*pi*x/T)./(2*pi*x/T);
   y = y.^m;
   graw( sprintf("set xtics('-3T' %d, '-2T' %d, '-T' %d, '0' 0, 'T' %d, '2T' %d, '3T' %d ) \n",-3*T,-2*T,-T,T,2*T,3*T) );

   graw( "set ytics(");
   for i=1:20
      graw( sprintf("'%dT' %d, ",i,i));
   endfor
   graw( sprintf("'%dT' %d) \n",i,i));

   plot_fx(x, y, sprintf("sinc^%d",m), sprintf("sinc_p%d",m));
   gset xtics autofreq;
   gset ytics autofreq;
endfunction

#----------------------------------------------------------------------------
#function: plot characteristic function
#  N:   number of points
#----------------------------------------------------------------------------
function plot_chi(N)
   T = 1;
   x = [-1.5*1/T: 3/N : +1.5*1/T];
   for i=1:N
      if( (x(i)>=-1/(2*T)) & (x(i)<1/(2*T)) ) y(i)=1;
      else                                    y(i)=0;
      endif
   endfor
   x=x(1:N);
   y=y(1:N);
   gset xlabel "f"
   graw( sprintf("set xtics('-3/2T' %f, '-1/T' %f, '-1/2T' %f, '0' 0, '1/2T' %f, '1/T' %f, '3/2T' %f ) \n",-3/(2*T),-1/T,-1/(2*T),1/(2*T),1/T,3/(2*T) ) );
   plot_fx(x, y, sprintf("pulse"), sprintf("pulse"));
   gset xtics autofreq;
   gset ytics autofreq;
endfunction

#----------------------------------------------------------------------------
#function: plot characteristic function shifted sum
#  N:   number of points
#----------------------------------------------------------------------------
function plot_chi_ss(N)
   T = 1;
   x = [-5/2: 5/N : +5/2];
   y = ones(1:N);
   y(round(0*N/5+1)) = 0;
   y(round(1*N/5  )) = 0;
   y(round(2*N/5  )) = 0;
   y(round(3*N/5  )) = 0;
   y(round(4*N/5  )) = 0;
   y(round(5*N/5  )) = 0;
   x=x(1:N);
   y=y(1:N);
   gset xlabel "f"
   graw( sprintf("set xtics('-3/2T' %f, '-1/T' %f, '-1/2T' %f, '0' 0, '1/2T' %f, '1/T' %f, '3/2T' %f ) \n",-3/(2*T),-1/T,-1/(2*T),1/(2*T),1/T,3/(2*T) ) );
   plot_fx(x, y, sprintf("pulse_ss"), sprintf("pulse_ss"));
   gset xtics autofreq;
   gset ytics autofreq;
endfunction


#----------------------------------------------------------------------------
#function: Plot eps file
#----------------------------------------------------------------------------
function eps_plot(filename)
   graw("set terminal postscript eps enhanced color blacktext solid linewidth 3 rounded \n");
   if( length(findstr(filename,"."))==0)
      filename = [filename,".eps"];
   endif
   printf("   output file name = \"%s\" \n",filename);
   fflush(stdout);
   graw(sprintf("set output '%s' \n",filename)); 
   graw("replot\n");
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
   graw("set xzeroaxis \n")
   graw("unset yzeroaxis \n")
   gset view 60,35
   gset xlabel "n"
   gset ylabel "alpha"
   graw("unset hidden3d \n")
   graw("unset contour \n")
   graw("unset border \n")
endfunction




#-----------------------------------
# Plot f(x)
#-----------------------------------
function plot_fx(x, fx, f_title, f_file)
   N = length(fx);
   if( length(x)!= N )
      x = [0:N-1];
   endif
   max_fx = max(fx);
   min_fx = min(fx);
   max_plot = max_fx*1.2;
   min_plot = min([-0.2 min_fx])*1.2;
   data  = fopen("cf.dat","w");
   for i=1:N
      fprintf(data ,"%10f  %10f  \n",x(i),fx(i) ); 
   endfor
   fflush(data);
   fclose(data);

   title(f_title);
   #gset xlabel "t"
   gset ylabel
   graw(sprintf("plot [][%f:%f] 'cf.dat' with lines linestyle 1 \n",min_plot,max_plot));
   eps_plot(f_file);
endfunction

#======================================
# Procedural Blocks
#======================================
gnuplot_setup();

#-----------------------------------
# Parameters
#-----------------------------------
delay_sec=1;
N = 200                                 #data size

#-----------------------------------
# Plot
#-----------------------------------
plot_sinc(N,4);
#plot_chi(N);
#plot_chi_ss(4*N);
#plot_Bspline(N,2);


#======================================
# End Processing
#======================================










