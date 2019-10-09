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
 require(stats);
 require(graphics);
 require(bspec);    # https://www.rdocumentation.org/packages/bspec/
 data(sunspots);
#---------------------------------------
# load data
#---------------------------------------
 x = sunspot.month;
#---------------------------------------
# calculate average
#---------------------------------------
 N = length(x);
 avg = sum(x) / N;
 print(avg);
#---------------------------------------
# calculate average recursively
#---------------------------------------
 muhat = seq(from=0, to=0, length=N);
 muhat[1] = x[1];
 w = seq(from=2, to=N, length=N-1);
N=100;
 muhat[2:N] = muhat[1:N-1] + (x[2:N] - muhat[1:N-1]) / w[1:N-1];
print(muhat[1:N]);
#---------------------------------------
# plot data
#---------------------------------------
 plot(x, col="blue", type="l");
 lines(c(start(x)[1],end(x)[1]), c(avg,avg), col="red", type="l");
 lines(muhat, col="green", type="l");
#---------------------------------------
# load data
#---------------------------------------
 class(x); # make sure x is of type time-series ("ts")
# x_psd = bspec::welchPSD(x, N/4);
#plot(x_psd$power, type="l", xlim=c(1850,2019), ylim=c(0,100), col="red")
#---------------------------------------
# process data
#---------------------------------------
#x_smooth = smooth(x, kind="3R");
#---------------------------------------
# plot
#---------------------------------------
#plot(mydata, col="blue");
#lines(x, col="blue")
#lines(x_smooth, col="red")

