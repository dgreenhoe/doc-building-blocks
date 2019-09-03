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
   h = 2*h/sum(h);                  #scale to sqrt(2)*sqrt(2) 
   m =  length(h);                  #
   p = ones(1,d*m)/m;               #<p,1>=1   ref: Rao p.53
   hu = upSample(h,d);              #upsample h(n) to match phi(t) density
   printf("\n");
   for i = 0:(k-1)                  #iterate
      printf("%d ",i);
      fflush(stdout);
      ph  = conv(hu,p);             #convolve
      p   = downSample(ph,2);       #downsample
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
#function: Generate Daubechies scaling coefficients
#  N:   number of non-zero scaling coefficients h(n)
#  h:   scaling coefficients
#  rR:  roots of R(z)R(z^-1) (the largest N/2 of these are discarded)
#  rH:  roots of H(z)        (N-1 roots)
#  reference:  Burrus page 260, Appendix C
#----------------------------------------------------------------------------
function [h,rR,rH] = gen_h(N)
   a=1;             #
   p2=1;            #
   q=1;             #
   H=sqrt(2)*[1 1]/2; # H(z) = sqrt(2)*(z+1)/2
   p=N/2;           # p = number of zeros at z=-1 (w=pi)

   for k=1:(p-1)
      H=conv(H,[1,1]/2);       # H(z) = H(z) x (z+1)/2
      a = -a*0.25*(k+p-1)/k;
      p2=conv(p2,[1,-2,1]);
      q=[0 q 0] + a*p2;
   endfor

   rR  = sort(roots(q));
   H  = conv(H,real(poly(rR(1:p-1))));
   rH = roots(H);
   h = H*sqrt(2)/sum(H);
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


#======================================
# Procedural Blocks
#======================================
gnuplot_setup();

#-----------------------------------
# Parameters
#-----------------------------------
delay_sec=1;
N = 1024                                 #data size
x  = 18*[0:N-1]/(N-1)-4;
basis = "D";
Nh=16;                    # N: number of non-zero coefficients in h(n)

#-----------------------------------
# Generate coefficients
#-----------------------------------
[h,rR,rH] = gen_h(Nh);    # generate scaling coefficients h(n)
g         = h2g_coefs(h); # generate wavelet coefficients g(n)
d = N/length(h);                       #density=samples/unit(1)
phi = gen_phi(h,  16,d);               #generate phi(x) from h(n) -- scaling function
psi = gen_psi(phi,g,d);                #generate psi(x) from g(n) -- wavelet function

#data_W(h,rR,rH, sprintf("%s%d",basis,Nh), sprintf("%s%d.dat",basis,Nh)  );
M = length(phi);
xx = (M/d)*[0:M-1]/(M-1);

plot_hn(h,  sprintf("%s%d Scaling Coefficients h(n)",basis,Nh), sprintf("%s%d_hn",basis,Nh) );   pause(delay_sec);
plot_hn(g,  sprintf("%s%d Wavelet Coefficients g(n)",basis,Nh), sprintf("%s%d_gn",basis,Nh) );   pause(delay_sec);
plot_fx(xx, phi, sprintf("%s%d phi(x) from h(n)",basis,Nh), sprintf("%s%d_phi",basis,Nh)); pause(delay_sec);
plot_fx(xx, psi, sprintf("%s%d psi(x) from g(n)",basis,Nh), sprintf("%s%d_psi",basis,Nh)); pause(delay_sec);


#======================================
# End Processing
#======================================










