printf(" ============================================ \n");
printf("| Daniel Greenhoe                            |\n");
printf(" ============================================ \n");

#======================================
# Initialization
#======================================
clear;

#======================================
# Function Declarations
#======================================

#----------------------------------------------------------------------------
#function: unwrap
#----------------------------------------------------------------------------
function y = unwrap(x)
   n = length(x);
   y = zeros(size(x));
   y(1:n/2) = x(n/2+1:n);
   y(n/2+1:n) = x(1:n/2);
endfunction

#----------------------------------------------------------------------------
#function: Plot eps file
#----------------------------------------------------------------------------
function eps_plot(filename)
   gset terminal postscript portrait enhanced color
   if( length(findstr(filename,"."))==0)
      filename = [filename,".eps"];
   endif
   output_cmd = sprintf("set output '%s' \n",filename);
   printf("   output file name = \"%s\" \n",filename);
   fflush(stdout);
   graw(output_cmd); 
   graw("replot \n");
   graw("set output \n");
   graw("set terminal windows \n");
endfunction

#----------------------------------------------------------------------------
# GNU-Plot Setup
#----------------------------------------------------------------------------
function gnuplot_setup()
   graw("\n\n");
   graw("#-----------------------------------\n");
   graw("# New plot                          \n");
   graw("#-----------------------------------\n");
   graw("clear \n");
   graw("reset \n");
   graw("unset key \n");
   graw("unset grid \n");
   graw("set terminal windows \n");
   graw("set style data lines 1 linewidth 3\n");
   gset xzeroaxis
   gset yzeroaxis
   gset ticslevel 0
   graw("unset label \n");
   graw("unset yzeroaxis \n");
   gset view 43,72
   gset xlabel "f"
   graw("unset hidden3d \n");
   graw("unset contour \n");

   graw("unset title \n");
   gset xlabel "t (time)"
   gset ylabel "w (outcome)"
   gset zlabel "x(t,w)"

   gset key top right box
   graw("unset key \n");

endfunction

#-----------------------------------
# Plot random process x(t,w)
#-----------------------------------
function plot_rand_norm()
   N=50;
   w=[1 : 6];
   x=[0:N-1];
   M=length(w);
   randn("seed",1.234);

   data = fopen("dataf","w");
   for i=1:M
      y = randn(1,N);
      for j=1:N
         fprintf(data,"%10f  %10f  %10f \n", x(j), w(i), y(j) ); 
      endfor
      fprintf(data,"\n" ); 
      fprintf(data,"\n" ); 
   endfor

   fflush(data);
   fclose(data);

   graw("splot [][][] 'dataf' every 1 with lines \n"); 
 
   eps_plot("x(t,w)");
   gset xtics autofreq
   gset ytics autofreq
endfunction


#======================================
# Procedural Blocks
#======================================
gnuplot_setup();
plot_rand_norm();

#======================================
# End Processing
#======================================





