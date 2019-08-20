#============================================================================
# Daniel J. Greenhoe
# R script file
# setwd("c:/dan/r/R");
# setwd("c:/dan/doc-building-blocks/R");
# dir()
# source("test.R")
#============================================================================
#---------------------------------------
# load add-on packages
#---------------------------------------
require(stats);
require(graphics);
#install.packages("car");
require("car");
#install.packages("bspec");
require("bspec")
#---------------------------------------
# load data
#---------------------------------------
data(BJsales);
x = sunspot.month;
#---------------------------------------
# load data
#---------------------------------------
# https://www.kaggle.com/brijbhushannanda1979/bigmart-sales-data
#mydata = read.csv(file="bigmart-sales-data-train.csv");
#x=mydata$Item_Outlet_Sales;

#https://www.kaggle.com/shenba/time-series-datasets/version/1
#mydata = read.csv(file="sales-of-shampoo-over-a-three-ye.csv");
#x = mydata$Sales.of.shampoo.over.a.three.year.period

#https://www.kaggle.com/shenba/time-series-datasets/version/1
# mydata = read.csv(file="Electric_Production.csv");
# x = mydata$IPG2211A2N
#x_ts = as.ts(x)
#x_psd = welchPSD(x_ts, length(x_ts)/4)
#plot(x_psd$power, type="l", ylim=c(0,100), col="blue")

#https://www.kaggle.com/nodarokroshiashvili/time-series
 mydata = read.csv(file="TS.csv");
 x = mydata$price
 x = x-mean(x);
 x_ts = as.ts(x)
 x_psd = welchPSD(x_ts, length(x_ts)/4)
 plot(x_psd$power, type="l", ylim=c(0,2.5e9), col="green")

 print(names(mydata))
#---------------------------------------
# process data
#---------------------------------------
x_smooth = smooth(x, kind="3R");
#---------------------------------------
# plot
#---------------------------------------
#plot(mydata, col="blue");
#lines(x, col="blue")
#lines(x_smooth, col="red")

