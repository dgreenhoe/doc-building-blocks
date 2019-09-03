printf(" ============================================ \n");
printf("| Daniel Greenhoe                            |\n");
printf("| Daubechies-4 plots                         |\n");
printf(" ============================================ \n");

#======================================
# Initialization
#======================================
clear;

#======================================
# Function Declarations
#======================================

#----------------------------------------------------------
# References:
# 
#    C. Sidney Burrus, Ramesh A. Gopinath, and Haitao Guo
#    Introduction to Wavelets and Wavelet Transforms; A Primer
#    Prentice Hall, 1998
#    ISBN 0-13-489600-9
#
#    Stephane G. Mallat 
#    A Wavelet Tour of Signal Processing 2nd edition. 
#    Academic Press. 1999 
#    ISBN 0-12-466606-X 
#    http://www.apnet.com/fatbrain/mallat_ch1.pdf
#
#    Raghuveer M. Rao, Ajit S. Bopardikar 
#    Wavelet Transforms; Introduction to Theory and Applications
#    Addison-Wesley.  1998 
#    ISBN 0-201-63463-5 
#    http://www.awl.com/cseng 
#----------------------------------------------------------

#======================================
# Initialization
#======================================
clear;

#======================================
# Function Declarations
#======================================

#----------------------------------------------------------------------------
#function: downsample x(n) by L
#----------------------------------------------------------------------------
function y = downSample(x,L)
   n = length(x);
   m = ceil(n/L);
   y = 0*x(1:m);
   for i=0:(m-1)
      y(i+1) = x(i*L+1);
   endfor
endfunction

#----------------------------------------------------------------------------
#function: upsample x(n) by L
#----------------------------------------------------------------------------
function y = upSample(x,L)
   n = length(x);
   m = n*L;
   y = zeros(1,m);
   for i=0:(n -1)
      y(i*L+1) = x(i+1);
   endfor
endfunction

#----------------------------------------------------------------------------
#function: Generate scaling function from scaling filter coefficients h(n)
#   p = scaling function
#   h = scaling filter coefficients h(n)
#   k = number of iterations
#   d = density = number phi(t) samples per h(n) sample
#
# reference: Burrus page 67
#    MatLab code:    Appendix C page 258
#    Key equation:   page 67
# 
# reference: Rao
#
#----------------------------------------------------------------------------
function p = gen_phi(h,k,d)
   h  = 2*h/sum(h);                  #scale to sqrt(2)*sqrt(2) 
   m  =  length(h);                  #
   p  = ones(1,d*m)/m;               #<p,1>=1   ref: Rao p.53
   hu = upSample(h,d);               #upsample h(n) to match phi(t) density
   printf("\n");
   for i = 0:(k-1)                   #iterate
      printf("%d ",i);
      fflush(stdout);
      ph  = conv(hu,p);              #convolve
      p   = downSample(ph,2);        #downsample
   endfor  
   p = p(1:m*d); 
   printf("\n");
endfunction

#----------------------------------------------------------------------------
#function: Generate wavelet function psi(t) from scaling function phi(t) 
#          and wavelet filter coefficients g(n)
#   phi = scaling function phi(t)
#   g = wavelet filter coefficients g(n)
#   d = density = number phi(t) samples per g(n) sample
# reference: Burrus page 15
# 
#----------------------------------------------------------------------------
function psi = gen_psi(phi,g,d)
   g   = sqrt(2)*g;                  #scale to sqrt(2)*sqrt(2) 
   gu  = upSample(g,d);              #upsample h(n) to match phi(t) density
   pg  = conv(gu,phi);               #convolve
   psi = downSample(pg,2);           #downsample
endfunction

#----------------------------------------------------------------------------
#function: Generate wavelet coefficients g(n) from scaling coefficients h(n)
#  references: 
#     Mallat page 238
#     Burrus page  15
#     Burrus page  79
#----------------------------------------------------------------------------
function g = h2g_coefs(h)
   g = h;
   N = length(h);
   for n=0:(N-1)  
      g(n +1) = (-1)**(n)*h(N-1-n +1);  
   endfor
endfunction

#----------------------------------------------------------------------------
#function: Generate filter coefficients from wavelet/scalar coefficients
#   _
#   h(n) = h(-n) with shift by N where N = length of h(n)
#   _
#   g(n) = g(-n) with shift by N
#----------------------------------------------------------------------------
function hbar = h2hbar(h)
   hbar = zeros(size(h));              #allocate memory
   N    = length(h)     ;              #N = length of h(n)
   for n = 0:(N-1)                     #
      hbar(n +1) = h(N-1-n +1);        #hbar(n) = h(N-1-n)
   endfor
endfunction

#----------------------------------------------------------------------------
#function: Fast Wavelet Transform                             Daniel Greenhoe
#   x = input sequence of length n (power of 2)
#   w = output sequence
#   h = scaling coefficients
#   W = array of output wavelet functions w0, w1, w2, ..., wm-1; m = log2(n)
#       (upper triangular matrix)
#   v0= scalar function
#
# m rows               <-- n/2 columns -->
# ------ ________________________________________________________________
#   0   | w0 |<-- length=1       0 ... 0                                 |
#       |____|____                                                       |
#   1   | w1      |<-- length=2                                          |
#       |_________|                      0 ... 0                         |
#       |     o                                                          |
#       |     o                                                          |
#       |     o                                                          |
#       |___________________                 0 ... 0                     |
#   m-3 | w_m-3  length=n/8 |                                            |
#       |___________________|____________                                |
#   m-2 | w_m-2  length=n/4              |           0 ... 0             |
#       |________________________________|_______________________________|
#   m-1 | w_m-1  length=n/2                                              |
#       |________________________________________________________________|
#
#----------------------------------------------------------------------------
function [W,v0] = fwt(y,h)
   N = length(y);                      #n  = length of input and output vectors
   m = log(N)/log(2);                  #m  = log2(n)
   g = h2g_coefs(h);                   #gen wavelet coefs g(n) from scalar coefs h(n)
   hbar = h2hbar(h);                   #generate scalar  filter coefficients
   gbar = h2hbar(g);                   #generate wavelet filter coefficients
   
   W = zeros(m,N/2);                   #initialize W with 0s: m x n/2 matrix
   if((N-2**m)!=0)
      printf("ERROR: length must be power of 2 \n");
      return;
      endif
   v = y;
   row = m;                 
   n = N/2;
   L = 1;
   while(n>=1)
      w = conv(gbar,v);                #filter with wavelet coefficients g(n)
      w = w(1:2:length(w));            #downsample by 2 
      w = w(1:n);                      #truncate length
      W(row,1:n) = w;                  #update output vector
      v = conv(hbar,v);                #filter with scaling coefficients h(n)
      v = v(1:2:length(v));            #downsample by 2
      v = v(1:n);                      #truncate length
      n  = n/2;                        #number of data points gets smaller
      row = row-1;                     #decrement row
      L = L*2;                         #upsample factor
   endwhile
   v0 = v(1);
endfunction

#----------------------------------------------------------------------------
#function: Generate length-4 scaling coefficients h(n) from parameter alpha
#  note: if alpha=pi/3, Daubechies-4 scaling coefficients are generated
#  reference:  Burrus page  66
#----------------------------------------------------------------------------
function h = gen_h4(alpha)
h =  1/(2*sqrt(2)) * [\
     1 - cos(alpha) + sin(alpha) ;
     1 + cos(alpha) + sin(alpha) ;
     1 + cos(alpha) - sin(alpha) ;
     1 - cos(alpha) - sin(alpha) ;
    ];
endfunction

#----------------------------------------------------------------------------
#function: Generate scaling length-2 coefficients h(n)
#  reference:  Burrus 
#----------------------------------------------------------------------------
function h = gen_h2(alpha)
   h = (sqrt(2)/2)*[1 1];
endfunction

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
   graw("set terminal postscript eps enhanced color blacktext solid linewidth 3 rounded \n");
   if( length(findstr(filename,"."))==0)
      filename = [filename,".eps"];
   endif
   output_cmd = sprintf("set output '%s' \n",filename);
   printf("   output file name = \"%s\" \n",filename);
   fflush(stdout);
   graw(output_cmd); 
   replot
   gset output
   gset terminal windows
endfunction

#----------------------------------------------------------------------------
#function: f(x)
#              ___.__
#          ___|      |      ____
#                    |_____|
#----------------------------------------------------------------------------
function y = f1(x)
   y = x;
   N = length(x);
   for i=1:N
      z = x(i);
      if    ( z<0 ) y(i) =   0.0;
      elseif( z<3 ) y(i) =  19.8;
      elseif( z<8 ) y(i) =  20.0;
      elseif( z<12) y(i) = -25.0;
      else          y(i) =   0.0;
      endif
   endfor
endfunction

#----------------------------------------------------------------------------
#function: f(x)  square wave
#              ___________
#          ___|           |____
#----------------------------------------------------------------------------
function y = f_sq(x)
   y = x;
   N = length(x);
   for i=1:N
      z = x(i);
      if    ( z<0  ) y(i) =   0.0;
      elseif( z<10 ) y(i) =  10;
      else          y(i) =   0.0;
      endif
   endfor
endfunction

#----------------------------------------------------------------------------
#function: f(x)  triangle
#              /\
#          ___/  \___
#----------------------------------------------------------------------------
function y = f_tri(x)
   y = x;
   N = length(x);
   for i=1:N
      z = x(i);
      if    ( z<0 ) y(i) =  0.0;
      elseif( z<5 ) y(i) =    2*z;
      elseif( z<10) y(i) = 20-2*z;
      else          y(i) =  0.0;
      endif
   endfor
endfunction

#----------------------------------------------------------------------------
#function: f(x)  parabola
#          f(x) = (-2/5)x^2 + 4x + 0
#----------------------------------------------------------------------------
function y = f_par(x)
   y = x;
   N = length(x);
   for i=1:N
      z = x(i);
      if    ( z<0 ) y(i) =  0.0;
      elseif( z<10) y(i) =  -0.4*z*z + 4*z + 0;
      else          y(i) =  0.0;
      endif
   endfor
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
   gset xzeroaxis
   gset ticslevel 0
   graw("unset label \n")
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
   min_plot = min([-1/1.2 min_fx])*1.2;
   data  = fopen("cf.dat","w");
   for i=1:N
      fprintf(data ,"%10f  %10f  \n",x(i),fx(i) ); 
   endfor
   fflush(data);
   fclose(data);

   title(f_title);
   gset xlabel "x"
   gset ylabel
   graw(sprintf("plot [][%f:%f] 'cf.dat' with lines linestyle 1 \n",min_plot,max_plot));
   eps_plot(f_file);
endfunction

#-----------------------------------
# Plot fft
#-----------------------------------
function plot_fft(x, fx, f_title, f_file)
   N = length(fx);
   if( length(x)!= N )
      x = [0:N-1];
   endif
   max_fx = max(fx);
   min_fx = min(fx);
   max_plot = max_fx*1.2;
   min_plot = min([-1/1.2 min_fx])*1.2;
   data  = fopen("cf.dat","w");
   for i=1:N
      fprintf(data ,"%10f  %10f  \n",x(i),fx(i) ); 
   endfor
   fflush(data);
   fclose(data);

   title(f_title);
   gset xlabel "frequency"
   gset ylabel 
   graw(sprintf("plot [][%f:%f] 'cf.dat' with impulses \n",min_plot,max_plot));
   eps_plot(f_file);
endfunction

#-----------------------------------
# Plot coefficients h(n)
#-----------------------------------
function plot_hn(h, h_title, h_file )
   max_h = max(h);
   min_h = min(h);
   max_plot = ceil(max([ abs(max_h) abs(min_h)] )*1.1);

   data  = fopen("cf.dat","w");
   Nh = length(h);
   for i=1:Nh
      fprintf(data ,"%10f  %10f  \n",i-1,h(i) ); 
   endfor
   fflush(data);
   fclose(data);

   title(h_title);
   gset xlabel "n"
   gset ylabel
   gset noyzeroaxis
   gset xtics 1
   for n=1:Nh
      graw(sprintf("set label '%f' at %f,%f \n",h(n),(n-1)+0.1,h(n)*1));
   endfor
   graw(sprintf("plot [-1:%d][%f:%f] 'cf.dat' with impulses,    \
                                     'cf.dat' with points  \n", \
                                     Nh,-max_plot,max_plot));
   eps_plot(h_file);
   gset nolabel
   gset xtics autofreq
endfunction

#-----------------------------------
# Plot Wavelet transform
#-----------------------------------
function plot_W(x, fx, W, v0, f_title, f_file)
   [M,N] = size(W);
   n=1;
   L=N;

   x  = x( 1:2:length( x));
   fx = fx(1:2:length(fx));   
   data2 = fopen("f(x)","w");
   for j=-1:M-1
   for k=1:length(fx)/32:length(fx)
      fprintf(data2,"%10f  %d  %10f \n",x(k), j, fx(k) ); 
   endfor
   
   fprintf(data2,"\n" ); 
   endfor
   fflush(data2);
   fclose(data2);
   
   data = fopen("Wjk","w");
   for j=0:M-1
      for k=1:n
         fprintf(data,"%10f  %3d  %10f \n",x((k-1)*L+1), j, sqrt(2**j)*abs(W(j+1,k)) ); 
         for i=1:(L-1)
            fprintf(data,"%10f  %3d  0 \n",x((k-1)*L+1+i), j ); 
         endfor         
      endfor
      fprintf(data,"\n" ); 
      n=n*2;
      L=L/2;
   endfor

   fprintf(data,"%10f  -1  %10f \n",x(1), v0 ); 
   for k=2:N
      fprintf(data,"%10f  -1  0 \n",x(k) ); 
   endfor
   fflush(data);
   fclose(data);

   title(f_title);
   gset xlabel "Shift by x"
   gset ylabel "Subspace"
   gset zlabel "W(j,k)"
   
   gset xtics 1
   graw("set ytics (");
   for j=0:M-1
      graw( sprintf("'W%d' %d, ",j,j) );
   endfor
   graw( "'V0' -1 ) \n" );

   gset key top right box
   graw("splot [][][] 'Wjk' with impulses, 'f(x)' with lines \n"); 
   eps_plot(f_file);
   gset xtics autofreq
   gset ytics autofreq
   gset nokey
endfunction


#-----------------------------------
# Write wavelet transform data to file
#-----------------------------------
function data_W(x, fx, W, v0, f_title, f_file)
   [M,N] = size(W);

   x  = x( 1:2:length( x));
   fx = fx(1:2:length(fx));   
   
   data = fopen(f_file,"w");
   fprintf(data,"%s\n\n",f_title);
   fprintf(data,"v0 = %10f \n\n", v0 ); 

   fprintf(data,"    x          f(x)   ");
   for j=0:M-1
      fprintf(data,"   W%d     ",j);
   endfor
   printf("\n");
   
   
   for k=1:N
      fprintf(data,"%10f %10f", x(k), fx(k) ); 
      for j=0:M-1
         if(W(j+1,k)==0) fprintf(data,"    -     "); 
         else fprintf(data,"%10f", sqrt(2**j)*abs(W(j+1,k)) ); 
         endif
      endfor
      fprintf(data,"\n" ); 
   endfor

   fflush(data);
   fclose(data);

   printf("   output file name = \"%s\" \n",f_file);
   fflush(stdout);
endfunction



#======================================
# Procedural Blocks
#======================================
gnuplot_setup();

#-----------------------------------
# Parameters
#-----------------------------------
delay_sec=1;
N = 1024                                 #data size

#-----------------------------------
# Generate f(x)
#-----------------------------------
x  = 18*[0:N-1]/(N-1)-4;
#x  = 16*[0:N-1]/(N-1)-2;

#-----------------------------------
# Settings
#-----------------------------------
fx    = "psi";
basis = "D4";

#y  = f_par(x);

h   = gen_h4(pi/3);                   #Daubechies-4 scaling coefficients
#h   = gen_h2();                       #Haar scaling coefficients h(n)


#-----------------------------------
# Generate coefficients
#-----------------------------------
g    = h2g_coefs(h);                   #gen wavelet coefs g(n) from scalar coefs h(n)
hbar = h2hbar(h);                      #generate scalar  filter coefficients
gbar = h2hbar(g);                      #generate wavelet filter coefficients
d = N/length(h);                       #density=samples/unit(1)
phi = gen_phi(h,  16,d);               #generate psi(x) from h(n) -- scaling function
psi = gen_psi(phi,g,d);                #generate psi(x) from phi(x)
y = psi;
[W,v0] = fwt(y,h);                     #wavelet transform
M = length(phi);
xx = (M/d)*[0:M-1]/(M-1);
fft_psi = abs(fft(psi,8*N))(1:300);
fft_phi = abs(fft(phi,8*N))(1:300);

if(0)
#plot_fx(x, y, fx, fx);
 plot_hn(h,    sprintf("%s Scaling Coefficients h(n)",basis), sprintf("%s_hn",basis) );             pause(delay_sec);
 plot_hn(g,    sprintf("%s Wavelet Coefficients g(n)",basis), sprintf("%s_gn",basis) );             pause(delay_sec);
 plot_hn(hbar, sprintf("%s Scaling Filter Coefficients hbar(n)",basis), sprintf("%s_hbar",basis) ); pause(delay_sec);
 plot_hn(gbar, sprintf("%s Wavelet Filter Coefficients gbar(n)",basis), sprintf("%s_gbar",basis) ); pause(delay_sec);
 plot_fx(xx, phi, sprintf("%s phi(x) from h(n)",basis), sprintf("%s_phi",basis)); pause(delay_sec); pause(delay_sec); 
 plot_fx(xx, psi, sprintf("%s psi(x) from h(n)",basis), sprintf("%s_psi",basis)); pause(delay_sec); pause(delay_sec);
 plot_fft(xx, fft_psi, sprintf("Fourier Transform of %s psi(x)",basis), sprintf("%s_psi_fft",basis)); pause(delay_sec);
endif

 plot_fx(xx, phi, sprintf("%s scaling function phi(x)",basis), sprintf("%s_phi",basis)); pause(delay_sec);  
 plot_fx(xx, psi, sprintf("%s wavelet function psi(x)",basis), sprintf("%s_psi",basis)); pause(delay_sec);

# plot_W(x, y, W, v0, sprintf("Basis:%s,  f(x):%s,  coefs scaled by sqrt(2^j)",basis,fx), sprintf("%s_%s_fwt",basis,fx) );  pause(delay_sec);
#W(3,4)=0;
# plot_W(x, y, W, v0, sprintf("Basis:%s,  f(x):%s, W(3,4)=0, coefs scaled by sqrt(2^j)",basis,fx), sprintf("%s_%s_fwt0",basis,fx) ); pause(delay_sec);
#
#data_W(x, y, W, v0, sprintf("Basis:%s,  f(x):%s,  coefs scaled by sqrt(2^j)",basis,fx), sprintf("%s_%s.dat",basis,fx));




#======================================
# End Processing
#======================================





