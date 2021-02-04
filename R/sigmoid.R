#============================================================================
# Daniel J. Greenhoe
# R script file
# setwd("c:/dan/personal/r/R"); 
# source("sigmoid.R");
# Reference: https://math.stackexchange.com/questions/3939703/
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
 x = seq( from=-2.5, to=2.5, length=1000 )

#---------------------------------------
# Logistic sigmoid function
#---------------------------------------
linear = function(x,slope,yintercept)
{
   result = slope * x + yintercept
}

#---------------------------------------
# Logistic sigmoid function
#---------------------------------------
sigmoid_logistic = function(x,slope)
{
   b = 4 * slope
   result = 1 / (1 + exp(-b*x))
}

#---------------------------------------
# Logistic sigmoid function
#---------------------------------------
sigmoid_tanh = function(x,slope)
{
   b = 4 * slope
   result = (exp(x) - exp(-x)) / (exp(x) + exp(-x))
}

#---------------------------------------
# Logistic sigmoid function
#---------------------------------------
sigmoid_intlogp = function(x,slope,p)
{
   a = slope
  #b = 2 * slope / p
   b = slope;
   result = 1 - (1 / (1 + a*exp(-b*x)))^p 
}

#---------------------------------------
# Logistic sigmoid to power p
#---------------------------------------
sigmoid_logp = function(x,slope,p)
{
   a = 2^(1/p)-1
   b = slope * (1 + a)^(2 + p -1) / (a * p )
   result = (1 / (1 + a*exp(-b*x)))^p 
  #result = 1-(1 / (1 + a*exp(b*x)))^p 
}

#---------------------------------------
# Logistic sigmoid function
#---------------------------------------
sigmoid_tanhp = function(x,b,p)
{
   result = 0.5 * (((exp(b*x) - exp(-b*x)) / (exp(b*x) + exp(-b*x)))^p + 1)
}

#---------------------------------------
# Logistic sigmoid function
# https://books.google.com/books?id=OeUKAAAAYAAJ&pg=PA15
# https://archive.org/details/integralstable00peirrich/page/59/
# https://archive.org/details/integralstable00peirrich/page/81/
#---------------------------------------
sigmoid_inttanhp = function(x,b,p)
{
  #result = 1 - ((exp(b*x) - exp(-b*x)) / (exp(b*x) + exp(-b*x)))^p
   result = 1 - (tanh(b*x))^p
}

#---------------------------------------
# Logistic sigmoid function
# https://en.wikipedia.org/wiki/List_of_integrals_of_hyperbolic_functions
#---------------------------------------
sigmoid_inttanhpx = function(x,b,p)
{
   intb = 0
   for( n in seq(from=p, to=2, by=-2))
   {
     intb = intb + 1/(b*(n-1))*(-1)^(n-1)
   }
   total = 0
   for( n in seq(from=p, to=2, by=-2))
   {
     total = total + 1/(b*(n-1))*(tanh(b*x))^(n-1)
   }
   total = total-intb
   result = total
}

#---------------------------------------
# Logistic sigmoid function
# https://ia800806.us.archive.org/7/items/GradshteinI.S.RyzhikI.M.TablesOfIntegralsSeriesAndProducts/Gradshtein_I.S.%2C_Ryzhik_I.M.-Tables_of_integrals%2C_series_and_products.pdf
# 2.424 (3) page 119
#---------------------------------------
inttanhpxx = function(x,b,n)
{
   p = 2*n
   inta = 0
   for( k in c(1:n) )
   {
     inta = inta + (-1) / (2*n-2*k+1)
   }
   inta = inta / b
   intx = 0
   for( k in c(1:n) )
   {
     intx = intx + (tanh(b*x))^(2*n-2*k+1) / (2*n-2*k+1)
   }
   intx = intx / b
   result = intx - inta
}

#---------------------------------------
# Display
#---------------------------------------
 colors = c( "blue"   , "red", "orange", "black", "purple","green");
 traces = c( "sigmoid", "sigmoid^(1/2)", "sigmoid^(1/3)", "sigmoid^(1/4)", "sigmoid^(1/5)" );
# plot ( x, sigmoid_inttanhp(x,2,13) , col=colors[1], lwd=3, type='l', xlab="x", ylab="y" )
 plot ( x, sigmoid_inttanhp(x,1/(0.55-0.45),3090) , col=colors[1], lwd=3, type='l', xlab="x", ylab="y", ylim=c(-0.5,1.5) )
#lines ( x, sigmoid_inttanhpx(x,1/(0.55-0.45),6000) , col=colors[2], lwd=3, type='l' )
 lines ( x, inttanhpxx(x,1/(0.55-0.45),3090) , col=colors[3], lwd=3, type='l' )
# plot ( x, sigmoid_inttanhpx(x,10,10) , col=colors[2], lwd=3, type='l')
#lines ( x, sigmoid_tanhp(x,1,3) , col=colors[2], lwd=3, type='l', xlab="x", ylab="y" )
 #lines( x, linear(x,1,0)       , col=colors[3], lwd=2, type='l' )
# plot ( x, sigmoid_tanh(x,1) , col=colors[1], lwd=3, type='l', xlab="x", ylab="y", ylim=c(-1,1), xaxp  = c(-1, 1, 20), yaxp=c(0,1,10))
#plot ( x, sigmoid_logistic(x,1) , col=colors[1], lwd=3, type='l', xlab="x", ylab="y", ylim=c(0,1), xaxp  = c(-1, 1, 20), yaxp=c(0,1,10))
#lines( x, sigmoid_logp(x,1,1/2) , col=colors[2], lwd=2, type='l' )
#lines( x, sigmoid_logp(x,1,1/3) , col=colors[3], lwd=2, type='l' )
#lines( x, sigmoid_logp(x,1,1/4) , col=colors[4], lwd=2, type='l' )
#lines( x, sigmoid_logp(x,1,1/10), col=colors[5], lwd=2, type='l' )
 legend("topleft", legend=traces, col=colors, lwd=3, lty=1:1)
 grid(col="green", lwd=2, lty="dashed",nx=NULL)
# grid(nx=11, ny=11)

