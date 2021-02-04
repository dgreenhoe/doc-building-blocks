#============================================================================
# Daniel J. Greenhoe
# R script file
# setwd("c:/dan/personal/r/R");
# source("plantGrowth.R");
# Reference: https://math.stackexchange.com/questions/3990086/
#============================================================================
#---------------------------------------
# packages
#---------------------------------------
#install.packages("stats");
#install.packages("R.utils");
#install.packages("rootSolve");
 require(stats);
 require(R.utils);
 require(rootSolve);
 rm(list=objects());

#---------------------------------------
# Data
#---------------------------------------
 tdata = c(0:10)
 ydata = c(18,33,56,90,130,170,203,225,239,247,251)
 t     = seq( from=min(tdata), to=max(tdata), length=1000 )

#---------------------------------------
# Estimate Function N(t)
#---------------------------------------
 N0 = ydata[1]
 N = function(t,N0,Nh,a0)
 {
   result = Nh / ( 1 + (Nh/N0-1)*exp(-a0*t) )
 }

#---------------------------------------
# Cost Function
#---------------------------------------
 cost = function(N0,Nh,a0)
 {
   summ = 0;
   for (i in c(1:length(tdata)))
   {
     summ = summ + ( N(tdata[i],N0,Nh,a0) - ydata[i] )^2
   }
   result = summ
 }

#---------------------------------------
# Partial derivative with respect to a0 of Cost Function
#---------------------------------------
 Pcosta0 = function(N0, Nh, a0)
 {
   summ = 0;
   for (i in c(1:length(tdata)))
   {
     summ = summ + ( N(tdata[i],N0,Nh,a0) )^2 *
                   ( N(tdata[i],N0,Nh,a0) - ydata[i] ) *
                   ( tdata[i] * exp(-a0*tdata[i]) )
   }
   result = summ
 }

#---------------------------------------
# Partial derivative with respect to Nh of Cost Function
#---------------------------------------
 PcostNh = function(N0, Nh, a0)
 {
   summ = 0;
   for (i in c(1:length(tdata)))
   {
     summ = summ + ( 1 - exp(-a0*tdata[i]) ) *
                   ( N(tdata[i],N0, Nh, a0) )^2 *
                   ( N(tdata[i],N0, Nh, a0) - ydata[i] )
   }
   result = summ
 }

#---------------------------------------
# Partial derivative with respect to Nh of Cost Function
#---------------------------------------
 PcostN0 = function(N0, Nh, a0)
 {
   summ = 0;
   for (i in c(1:length(tdata)))
   {
     summ = summ + ( exp(-a0*tdata[i]) ) *
                   ( N(tdata[i],N0, Nh, a0) )^2 *
                   ( N(tdata[i],N0, Nh, a0) - ydata[i] )
   }
   result = summ
 }

#---------------------------------------
# Partial derivative vector of cost
#---------------------------------------
Pcost = function(x)
{
   N0 = x[1]
   Nh = x[2]
   a0 = x[3]
   F1 = Pcosta0( N0, Nh, a0 );
   F2 = PcostNh( N0, Nh, a0 );
   F3 = PcostN0( N0, Nh, a0 );
   result = c(F1, F2, F3);
}

#---------------------------------------
# Calculate roots
#---------------------------------------
 Roots = multiroot( f=Pcost, start=c(ydata[1], ydata[11], 0.6) );
 N0 = Roots$root[1]
 Nh = Roots$root[2]
 a0 = Roots$root[3]

#---------------------------------------
# Display
#---------------------------------------
 printf("(N0, Nh, a0) = (%.10f, %.10f, %.10f) with estim.precis=%.2e\n", N0, Nh, a0, Roots$estim.precis )
 colors = c( "red" , "blue" );
 traces = c( "N(t)", "data" );
 plot ( t , N(t, N0, Nh, a0), col=colors[1], lwd=2, type='l', xlab="t", ylab="y", ylim=c(0,max(ydata)+10) )
 lines( tdata, ydata        , col=colors[2], lwd=5, type='p' )
 legend("topleft", legend=traces, col=colors, lwd=3, lty=1:1)
 grid()

