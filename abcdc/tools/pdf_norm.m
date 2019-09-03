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
#function: Normal distribution
#               ______
#              /      \
#          ___/        \___
#  u  = mean
#  v  = variance
#----------------------------------------------------------------------------
function [y,x] = f_pdf_norm(u,v,N)
   x=[-3 : 6/N : 3];
   x=x(1:N);
   y = (1/sqrt(2*pi*v))*exp(-(x-u).^2/(2*v));
#   y = (1/sqrt(2*pi*v))*ones(1,N);
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
   graw("set style data lines \n");
   gset xzeroaxis
   gset yzeroaxis
   gset ticslevel 0
   graw("unset label \n");
   graw("unset yzeroaxis \n");
   gset view 60,35
   gset xlabel "f"
   graw("unset hidden3d \n");
   graw("unset contour \n");
endfunction

#-----------------------------------
# Plot Normal distribution
#-----------------------------------
function plot_pdf_norm()
   N=200;
   v=[0.1 : 0.2 : 2];
   M=length(v)

   data = fopen("dataf","w");
   for i=1:M
      [y,x] = f_pdf_norm(0,v(i),N);
      for j=1:N
         fprintf(data,"%10f  %10f  %10f \n", x(j), v(i), y(j) ); 
      endfor
      fprintf(data,"\n" ); 
   endfor

   fflush(data);
   fclose(data);
 
#   title("Gaussian pdf");
   graw("unset title \n");
   gset xlabel "x"
   gset ylabel "variance"
   gset zlabel ""

   gset key top right box
   graw("unset key \n");
   graw("splot [][][] 'dataf' every 2 with lines \n"); 

   eps_plot("pdf_norm");
   gset xtics autofreq
   gset ytics autofreq
endfunction


#======================================
# Procedural Blocks
#======================================
gnuplot_setup();
plot_pdf_norm();

#======================================
# End Processing
#======================================





