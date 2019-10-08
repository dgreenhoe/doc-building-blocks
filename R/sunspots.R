#============================================================================
# Daniel J. Greenhoe
# R script file
# setwd("c:/dan/r/R");
# setwd("c:/dan/doc-building-blocks/R");
# dir()
# source("sunspots.R")
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
 N = length(x);
 avg = sum(x) / N;
 print(avg);
 mu = mean(x);
 plot(x, col="blue", type="l");
 lines(c(start(x)[1],end(x)[1]), c(avg,avg), col="red", type="l");
#---------------------------------------
# load data
#---------------------------------------
 class(x); # make sure x is of type time-series ("ts")
 x_psd = bspec::welchPSD(x, N/4);
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

