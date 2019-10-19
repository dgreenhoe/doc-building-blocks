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
mypsd = function(x, numSegments=4, dataDump=FALSE, dataPlot=TRUE, dataFile="mypsd.dat") 
  {
  xts       = as.ts(as.vector(x)); # year indices seems to confuse welchPSD
  N         = length(xts);         # length of time series
  Fs        = 12;                  # sample rate = 12 samples per year
  estMean   = mean(xts);           # estimated mean
  segLength = N / numSegments;     # segment length
  xpsd = bspec::welchPSD(xts - estMean, seglength = segLength);
  psdMax = max(xpsd$power);
  binMax = ramify::argmax(as.matrix(xpsd$power), rows = FALSE);
  freqMax = xpsd$frequency[binMax] * Fs; # dominate non-DC frequency
  periodT = 1 / freqMax;           # estimated period
  if(dataPlot){ 
    plot(xpsd$power, type="b", xlim=c(1,50), ylim=c(0,max(xpsd$power)), col="blue");
    }
  if(dataDump){
    sink(dataFile);
    print(sprintf("%%============================================================================="));
    print(sprintf("%% Daniel J. Greenhoe "                                                         ));
    print(sprintf("%% PSD data file suitable for use by LaTeX PStricks"                            ));
    print(sprintf("%% number of segments = %f", numSegments                                        ));
    print(sprintf("%% estimated period   = %12.6f", periodT                                        ));
    print(sprintf("%%============================================================================="));
    print(sprintf("["                                                                              ));
    for(i in 1:length(xpsd$power)) 
      print(sprintf("(%10.6f, %10.6f)", xpsd$frequency[i], xpsd$power[i]/psdMax               ));
    print(sprintf("]"                                                                              ));
    sink();
    }
  periodT;                         # return estimated period
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
 range_all = c(39:74);
 range_cagefree_org  = c(75:110);
 range_cagefree_non  = c(111:146);
 range_cagefree      = c(111:146);
 x    = read.csv("../common/datasets/egg-production_osf-hmfn3.csv", comment.char = "#");
 y    = read.csv("../common/datasets/cagefree-ratios_osf-6hty8.csv", comment.char = "#");

 year = lubridate::decimal_date(lubridate::ymd(x$observed_month))

 year_all = year[range_all];
 year_cagefree_org = year[range_cagefree_org];
 year_cagefree_non = year[range_cagefree_non];
 print(year_all - year_cagefree_org)
 print(year_cagefree_org - year_cagefree_non)

 ratio_all           = x$n_eggs[range_all] / x$n_hens[range_all];
 ratio_cagefree_org  = x$n_eggs[range_cagefree_org] / x$n_hens[range_cagefree_org];
 ratio_cagefree_non  = x$n_eggs[range_cagefree_non] / x$n_hens[range_cagefree_non];
 ratio_cagefree      = (ratio_cagefree_org + ratio_cagefree_non) / 2;


 plot( year[range_all], ratio_all, type="o", ylim=c(20.5,25));
#lines(year[range_cagefree_org], ratio_cagefree_org, type="o", col="red");
 lines(year[range_cagefree], ratio_cagefree, type="o", col="red");
 #lines(year[range_cagefree_non], (x$n_eggs[range_cagefree_non] / x$n_hens[range_cagefree_non]),type="o", col="blue");


 #for(n in c(1:length(year_all))){
 #  print(sprintf("  (%8f,   %8f)", year_all[n], ratio_all[n]));
 #  }

 for(n in c(1:length(y$ratio_hens))){
   print(sprintf("  (%8f,   %8f)", lubridate::decimal_date(lubridate::ymd(y$observed_month[n])), y$ratio_hens[n]));
   }
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
