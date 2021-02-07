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
# Logistic sigmoid to power p
#---------------------------------------
sigmoid_logp = function(x,slope,p)
{
   a = 2^(1/p)-1
   b = slope * (1 + a)^(2 + p -1) / (a * p )
   result = (1 / (1 + a*exp(-b*x)))^p 
}

#---------------------------------------
# tanh^p(bx)
#---------------------------------------
#tanhp = function(x,b,p)
#{
#   result = 0.5 * (((exp(b*x) - exp(-b*x)) / (exp(b*x) + exp(-b*x)))^p + 1)
#}

#---------------------------------------
# Logistic sigmoid function
# https://books.google.com/books?id=OeUKAAAAYAAJ&pg=PA15
# https://archive.org/details/integralstable00peirrich/page/59/
# https://archive.org/details/integralstable00peirrich/page/81/
#---------------------------------------
tanhp = function(x,b,p)
{
  #result = 1 - ((exp(b*x) - exp(-b*x)) / (exp(b*x) + exp(-b*x)))^p
   result = 1 - (tanh(b*x))^p
}

#---------------------------------------
# Logistic sigmoid function
# https://ia800806.us.archive.org/7/items/GradshteinI.S.RyzhikI.M.TablesOfIntegralsSeriesAndProducts/Gradshtein_I.S.%2C_Ryzhik_I.M.-Tables_of_integrals%2C_series_and_products.pdf
# 2.424 (3) page 119
#---------------------------------------
inttanhp = function(x,b,p)
{
   n = p/2
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
 colors = c( "blue", "red", "orange", "black", "purple","green");
#plot ( x, tanhp(x,2,1) , col=colors[1], lwd=3, type='l', xlab="x", ylab="y" )
#plot ( x, tanhp(x,1/(0.55-0.45),6180) , col=colors[1], lwd=3, type='l', xlab="x", ylab="y", ylim=c(-0.1,1.1), xlim=c(-1.5,1.5))
 plot ( x, inttanhp(x,1/(0.55-0.45),6180) , col=colors[2], lwd=3, type='l', xlab="x", ylab="y" )
#lines( x, linear(x,1,1/2)       , col=colors[6], lwd=2, type='l' )
#lines( x, linear(x,0,1/2)       , col=colors[6], lwd=2, type='l' )
 legend("topleft", legend="f(x); b=10; p=6180", col="red", lwd=3, lty=1:1)
 grid()

