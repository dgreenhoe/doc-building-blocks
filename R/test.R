#============================================================================
# Daniel J. Greenhoe
# R script file
# setwd("c:/dan/r/R");
# dir()
# source("test.R")
#============================================================================
require(stats);
require(graphics);
data(BJsales);
install.packages("car")
x = sunspot.month;
plot(x, col="white");
x_smooth = smooth(x, kind="3R");
lines(x_smooth, col="red")

