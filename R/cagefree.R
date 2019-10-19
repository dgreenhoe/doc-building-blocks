#============================================================================
# Daniel J. Greenhoe
# R script file
# setwd("c:/dan/r/R");
# setwd("c:/dan/doc-building-blocks/R");
# dir()
# source("eggs.R");
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
 require(graphics);
 require(datasets);
 require(ramify);
 require(bspec);     # https://www.rdocumentation.org/packages/bspec/
 require(lubridate); # https://www.rdocumentation.org/packages/lubridate/
#---------------------------------------
# load data
#---------------------------------------
 x    = read.csv("../common/datasets/cagefree-ratios_osf-6hty8.csv", comment.char = "#");
 ratio = x$ratio_hens;
 year  = lubridate::decimal_date(lubridate::ymd(x$observed_month));
 plot( year, ratio, col="blue", type='o' );

# for(n in c(1:length(y$ratio_hens))){
#   print(sprintf("  (%8f,   %8f)", lubridate::decimal_date(lubridate::ymd(y$observed_month[n])), y$ratio_hens[n]));
#   }
#---------------------------------------
# plot data
#---------------------------------------
# mydata(x);
# estT = mypsd(x, numSegments=3, dataDump=FALSE);
# print(sprintf("estimated period = %.16f years", estT));


#mypdf(x);
#pdfr(x);
#---------------------------------------
# process data
#---------------------------------------
