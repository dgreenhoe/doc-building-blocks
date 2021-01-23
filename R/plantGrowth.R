#============================================================================
# Daniel J. Greenhoe
# R script file
# setwd("c:/dan/personal/r/R");
# dir()
# source("plantGrowth.R");
# Reference: https://math.stackexchange.com/questions/3990086/
#============================================================================
#---------------------------------------
# install packages (perform once)
#---------------------------------------
#install.packages("bspec");
#install.packages("ramify");
#---------------------------------------
# load add-on packages
#---------------------------------------
 rm(list=objects());
 require(stats);
# require(graphics);
# require(datasets);
#---------------------------------------
# load add-on packages
#---------------------------------------

#---------------------------------------
# Data
#---------------------------------------
t  = seq( from=0, to=10, length=100 )
tn = c(0:10)
yn = c(18,33,56,90,130,170,203,225,239,247,251)

#---------------------------------------
# Estimate N(t) of Data
#---------------------------------------
N0=18
Nh=252
N = function(t,a) Nh / (1 + (Nh/N0-1)*exp(-a*t))
#g = function(x,n) 1 + 13*exp(-x*n)
#N = function(x,xn,yn) n * exp(-x*xn) / (1 + 13*exp(-x*xn))(yn-1)

#---------------------------------------
# Cost Function
#---------------------------------------
#cost = function(a0) sum((N(tn,a0)-yn)^2)
cost = function(a0) (N(tn[1],a0)-yn[1])^2 +(N(tn[2],a0)-yn[2])^2 +(N(tn[3],a0)-yn[3])^2 +(N(tn[4],a0)-yn[4])^2 +(N(tn[5],a0)-yn[5])^2 +(N(tn[6],a0)-yn[6])^2 +(N(tn[7],a0)-yn[7])^2 +(N(tn[8],a0)-yn[8])^2 +(N(tn[9],a0)-yn[9])^2 +(N(tn[10],a0)-yn[10])^2 +(N(tn[11],a0)-yn[11])^2
a0=seq(from=0.1, to=3.0, length=1000)
#plot(a0,cost(a0), col="blue", lwd=2, type='l')

#---------------------------------------
# Optimal a0 that minimizes Cost Function
#---------------------------------------
x=seq(from=0.3, to=0.9, length=1000)
Dcost = function(x) (N(tn[ 1],x))^2*(N(tn[ 1],x)-yn[ 1])*tn[1]*exp(-x*tn[ 1])+ (N(tn[ 2],x))^2*(N(tn[ 2],x)-yn[ 2])*tn[ 2]*exp(-x*tn[ 2])+ (N(tn[ 3],x))^2*(N(tn[ 3],x)-yn[ 3])*tn[ 3]*exp(-x*tn[ 3])+ (N(tn[ 4],x))^2*(N(tn[ 4],x)-yn[ 4])*tn[ 4]*exp(-x*tn[ 4])+ (N(tn[ 5],x))^2*(N(tn[ 5],x)-yn[ 5])*tn[5]*exp(-x*tn[ 5])+ (N(tn[ 6],x))^2*(N(tn[ 6],x)-yn[ 6])*tn[6]*exp(-x*tn[ 6])+ (N(tn[ 7],x))^2*(N(tn[ 7],x)-yn[ 7])*tn[7]*exp(-x*tn[ 7])+ (N(tn[ 8],x))^2*(N(tn[ 8],x)-yn[ 8])*tn[8]*exp(-x*tn[ 8])+ (N(tn[ 9],x))^2*(N(tn[ 9],x)-yn[ 9])*tn[9]*exp(-x*tn[ 9])+ (N(tn[10],x))^2*(N(tn[10],x)-yn[10])*tn[10]*exp(-x*tn[10])+ (N(tn[11],x))^2*(N(tn[11],x)-yn[11])*tn[11]*exp(-x*tn[11])
plot(x, Dcost(x), col="blue", lwd=2, type='l')
aOpt = uniroot( Dcost, c(0.3, 0.9) )
a0fixed=0.6631183
#---------------------------------------
# Graphics
#---------------------------------------
plot( t,  N(t,a0fixed), col="red" , lwd=2, type='l', xlab="t", ylab="y" ) 
lines(  tn, yn     , col="blue", lwd=3, type='p' )
#legend("topleft", legend=c("N(t)", "data"), col=c(3,4))
legend("topleft", legend=c("N(t)", "data"), col=c("red", "blue"), lwd=3, lty=1:1)
#legend("topleft", legend=c("N(t)", "data"))
grid()