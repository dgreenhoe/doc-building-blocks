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
#install.packages("ramify");
#---------------------------------------
# load add-on packages
#---------------------------------------
 rm(list=objects());
 require(stats);
 require(graphics);
 require(datasets);
 require(ramify);
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
 x = sunspot.month;
#---------------------------------------
# plot data
#---------------------------------------
 mydata(x);
 estT = mypsd(x, numSegments=3, dataDump=FALSE);
 print(sprintf("estimated period = %.16f years", estT));


#mypdf(x);
#pdfr(x);
#---------------------------------------
# process data
#---------------------------------------
