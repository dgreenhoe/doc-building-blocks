#============================================================================
# Daniel J. Greenhoe
# R script file
# setwd("c:/dan/r/R");
# setwd("c:/dan/doc-building-blocks/R");
# dir()
# source("sunspots.R");
#============================================================================
#---------------------------------------
# install packages (perform once)
#---------------------------------------
#install.packages("bspec");
#---------------------------------------
# load add-on packages
#---------------------------------------
 rm(list=objects());
 require(stats);
 require(graphics);
 require(datasets);
 require(bspec);    # https://www.rdocumentation.org/packages/bspec/
 data(sunspots, package="datasets");
#=======================================
# function: data
#=======================================
mydata = function(x)
  {
  N = length(x);
  avg = sum(x) / N;
  print(avg);
  plot(x, col="blue", type="l");
  lines(c(start(x)[1],end(x)[1]), c(avg,avg), col="red", type="l");
  x_smooth = smooth(x, kind="3R");
  lines(x_smooth, col="green", type="l");
  }
#=======================================
# function: PSD
#=======================================
mypsd = function(x) 
  {
 class(x); # make sure x is of type time-series ("ts")
 #x_psd = bspec::welchPSD(x, length(x)/4);
#plot(x_psd$power, type="l", xlim=c(1850,2019), ylim=c(0,100), col="red")
  }
#=======================================
# function: pdf
#=======================================
mypdf = function(x){
  plot(density(x, bw="SJ"))
  }
#=======================================
# function: pdf recursive
#=======================================
pdfr = function(x)
{
  N = length(x);
  muhat = seq(from=0, to=0, length=N);
  muhat[1] = x[1];
  w = seq(from=2, to=N, length=N-1);
  muhat[2:N] = muhat[1:N-1] + (x[2:N] - muhat[1:N-1]) / w[1:N-1];
  print(muhat[1:10]);
}
#---------------------------------------
# load data
#---------------------------------------
 x = sunspot.month;
#---------------------------------------
# plot data
#---------------------------------------
 mydata(x);
#mypsd(x);
#mypdf(x);
#pdfr(x);
#---------------------------------------
# process data
#---------------------------------------
