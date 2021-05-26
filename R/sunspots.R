#============================================================================
# Daniel J. Greenhoe
# R script file
# setwd("c:/dan/r/R");
# setwd("c:/dan/doc-building-blocks/R");
# setwd("c:/dan/personal/r/R");
# source("sunspots.R");
#============================================================================
#---------------------------------------
# install packages (perform once)
#---------------------------------------
#install.packages("bspec");
#install.packages("ramify");
#install.packages("freqdom");
#install.packages("matlib");
#install.packages("R.utils");
#---------------------------------------
# load add-on packages
#---------------------------------------
 rm(list=objects());
# require(stats);
# require(graphics);
# require(ramify);
# require(freqdom);
  require(bspec);    # welchPSD
# require(matlib);
  require(R.utils);  # printf
#data(sunspots, package="datasets"); # deprecated in favor of Silso data
#---------------------------------------
# Global parameters
#---------------------------------------
 author   = "Daniel J. Greenhoe"
 thisfile = "sunspots.R"
 baseName = "sunspots"
#------------------------------------------------------------------------------
# \brief   Estimate Auto-Correlation Function (ACF) of sunspot data minus estimated mean
# \returns data table
#------------------------------------------------------------------------------
sunspots_getData = function( dataDump    = FALSE,
                             dataPlot    = TRUE, 
                             dataFileIn  = "../data/silso_SN_m_tot_V2.0_20210524.csv", 
                             dataFileOut = "tex/sunspots.dat"
                           )
{
 #x = sunspot.month;
  x       = read.csv(file=dataFileIn, header=TRUE,sep=";", comment.char="#", strip.white=TRUE);
  xts     = as.ts(as.vector(x$count));
  dvect   = as.vector(x$date)
  cvect   = as.vector(x$count)
  if(dataPlot)
  {
    plot(x=dvect, y=cvect, lwd=2, type='l', main="Sunspot Estimation", col="blue", xlab="year", ylab="count");
    abline(h=seq(from=0,   to= 400, by=100), lty='dashed', col = "green")
    abline(v=seq(from=1750,to=2020, by= 10), lty='dashed', col = "green")
  }
  if(dataDump)
  {
    sink(dataFileOut);
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% Sunspot monthly mean data file suitable for use with LaTeX PStricks\n"           );
    printf("%% For an example, see \"%s.tex\"\n", baseName                                      );
    printf("%% This file auto-generated with \"%s\"\n", thisfile                                );
    printf("%% using data from \"%s\"\n", dataFileIn                                            );
    printf("%% Hand-editing not recommended\n"                                                  );
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(i in 1:length(dvect))
      printf("  (%12.8f, %12.8f)\n", dvect[i], cvect[i]                                         );
    printf("]\n"                                                                                );
    sink();
  }
  return(x);
}

#------------------------------------------------------------------------------
# \brief   Estimate Auto-Correlation Function (ACF) of sunspot data minus estimated mean
# \returns ACF data table
#------------------------------------------------------------------------------
sunspots_ACF = function( dataDump    = FALSE,
                         dataPlot    = TRUE, 
                         nLag        = 2000, 
                         dataIn, 
                         dataFileOut = "tex/sunspots_acf.dat"
                       )
{
 #x       = read.csv(file=dataFileIn, header=TRUE,sep=";", comment.char="#", strip.white=TRUE);
  x       = dataIn;
  xts     = as.ts(as.vector(x$count));
  estMean = mean(xts);
  a       = stats::acf(xts - estMean, type="correlation", lag.max=nLag, plot=FALSE)
  avect   = as.vector(a$acf)
  lvect   = as.vector(a$lag)/12

  if(dataPlot)
  {
    plot(x=lvect, y=avect, lwd=2, type='l', xaxp=c(0,160,16), yaxp=c(-0.4,1.0,14), col="blue",
         main="Sunspot Auto-Correlation Estimation", xlab="lag in years", ylab="magnitude");
    abline(h=seq(from=0, to=1.0, by= 0.1), lty='dashed', col = "green")
    abline(v=seq(from=0, to=160, by=10.0), lty='dashed', col = "green")
  }
  if(dataDump)
  {
    sink(dataFileOut);
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% Sunspot auto-correlation function (ACF) data file suitable for use with LaTeX PStricks\n" );
    printf("%% For an example, see \"%_acf.tex\"\n", baseName                                   );
    printf("%% This file auto-generated using \"%s\" --- hand-editing not recommended\n", thisfile);
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(i in 1:length(avect))
      printf("  (%12.8f, %12.8f)\n", lvect[i], avect[i]                                         );
    printf("]\n"                                                                                );
    sink();
  }
  return(a);
}

#------------------------------------------------------------------------------
# \brief   Estimate sunspot period using Welch Estimate of PSD
# \returns Estimated period
#------------------------------------------------------------------------------
sunspots_PSD = function( dataDump    = FALSE,
                         dataPlot    = TRUE, 
                         numSegments = 4, 
                         dataIn, 
                         dataFileOut = "tex/sunspots_psd.dat"
                       )
{
  xts       = as.ts(as.vector(dataIn$count));
  N         = length(xts);         # length of time series
  Fs        = 12;                  # sample rate = 12 samples per year
  estMean   = mean(xts);           # estimated mean
  segLength = N / numSegments;     # segment length
  xpsd = bspec::welchPSD(xts - estMean, seglength = segLength);
  psdMax = max(xpsd$power);
  binMax = ramify::argmax(as.matrix(xpsd$power), rows = FALSE);
  freqMax = xpsd$frequency[binMax] * Fs; # dominate non-DC frequency
  periodT = 1 / freqMax;           # estimated period
  printf("sunspots_PSD(x):\n");
  resultStr = sprintf("numSegments = %d\nf = %12.8f samples/year\nperiod = %12.8f years\n", numSegments, freqMax, periodT);
  printf("%s\n", resultStr);

  if(dataPlot)
  {
    titleStr = sprintf("0-average Sunspot Power Spectral Density (PSD) Estimate using welchPSD with numSegments=%d",numSegments);
    plot(x=xpsd$frequency*Fs, y=sqrt(xpsd$power/N), type="h", lwd=3, col="blue", xaxp=c(0,Fs/2,60), yaxp=c(0,20,5), xlab="", ylab="", main="", sub="");
    title(main=titleStr, xlab="samples/year (max=Fs/2=6 samples/year)", ylab="Gain (sqrt of power)");
    title(main=resultStr, line=-10);
    abline(v=seq(from=0, to=Fs/2, by=0.1), lty="dotted", col="green")
    abline(v=seq(from=0, to=Fs/2, by=1.0), lty="solid",  col="green")
    abline(h=seq(from=0, to=20.0, by=4.0), lty="dashed", col="green")
  }
  if(dataDump)
  {
    sink(dataFileOut);
    printf("%%=============================================================================\n");
    printf("%% %s \n", author                                                                 );
    printf("%% PSD data file suitable for use by LaTeX PStricks\n"                            );
    printf("%% number of segments          = %d\n",                  numSegments              );
    printf("%% estimated maximum frequency = %12.6f samples/year\n", freqMax                  );
    printf("%% estimated period            = %12.6f years\n",        periodT                  );
    printf("%% This file auto-generated using \"%s\" --- hand-editing not recommended\n", thisfile);
    printf("%%=============================================================================\n");
    printf("[\n"                                                                              );
    for(i in 1:length(xpsd$power))
      printf("  (%12.8f, %12.8f)\n", xpsd$frequency[i]*Fs, xpsd$power[i]/psdMax               );
    printf("]\n"                                                                              );
    sink();
  }
  return(periodT);
}

#------------------------------------------------------------------------------
# \brief Estimate sunspot period using PCA and Welch Estimate of PSD of PCA
# PCA
# https://cran.r-project.org/web/packages/matlib/vignettes/inv-ex1.html
# https://stat.ethz.ch/R-manual/R-patched/library/base/html/eigen.html
# https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix
#------------------------------------------------------------------------------
sunspots_PCA_PSD = function(x, numSegments=4)
{
  #x = read.csv(file='data/sunspots_silso_20210522.csv', header=TRUE,sep=",", comment.char="#", strip.white=TRUE)
  x = read.csv(file='../data/silso_SN_m_tot_V2.0_20210524.csv', header=TRUE,sep=";", comment.char="#", strip.white=TRUE)
  numSegments=4
  nLag = 2000
  dataDump = TRUE
  numSegments=4
  xts       = as.ts(as.vector(x$count)); # year indices seems to confuse welchPSD
  N         = length(xts);         # length of time series
  Fs        = 12;                  # sample rate = 12 samples per year
  estMean   = mean(xts);           # estimated mean
  segLength = N / numSegments;     # segment length
  a         = acf(xts - estMean, type="correlation", lag=nLag)
  avect     = as.vector(a$acf)
  lvect     = as.vector(a$lag)/12
  al        = cbind(lvect,avect)
  A         = stats::toeplitz(avect)
  Q         = eigen(A, symmetric=TRUE, only.values=FALSE)
  V         = Q$vectors
  L         = Q$values
  D         = diag(L)

  colors = c("blue", "red", "orange", "green", "purple", "brown", "black");
  traces = colors[1:5]
  for( n in 1:length(traces))
  {
    if(n==1) plot( L[n] * Q$vectors[,n], type='l', lwd=3, col=colors[n])
    else     lines(L[n] * Q$vectors[,n], type='l', lwd=3, col=colors[n])
    traces[n] = sprintf("Eigen Vector %2d", n)
  }
  legend("topleft", legend=traces, col=colors, lwd=3, lty=1:1)

  printf("sunspots_PCA_PSD(x) using Welch PSD:\n");
  for( n in 1:10 )
  {
    xpsd    = bspec::welchPSD(x=as.ts(L[n] * Q$vectors[,n]), seglength = segLength);
    psdMax  = max(xpsd$power);
    binMax  = ramify::argmax(as.matrix(xpsd$power), rows = FALSE);
    freqMax = xpsd$frequency[binMax] * Fs; # dominate non-DC frequency
    periodT = 1 / freqMax;                 # estimated period
    printf("Vector %2d (lambda=%10.6f) f=%8.6f samples/year period=%9.6f years\n", n, L[n], freqMax, periodT);
  }

  printf("sunspots_PCA_PSD(x) using DFT:\n");
  for( n in 1:20 )
  {
    x=as.ts(L[n] * Q$vectors[,n])
    xfft    = fft(x, inverse=FALSE);
    fftMax  = abs(xfft);
    binMax  = ramify::argmax(as.matrix(fftMax), rows = FALSE);
    freqMax = (binMax-1) * Fs / length(xfft);
    periodT = 1 / freqMax;                 # estimated period
    phase   = Arg(xfft[binMax]);
    degrees = phase / pi * 180
    printf("Vector %2d (lambda=%10.6f) f=%8.6f samples/year period=%9.6f years phase=%9.6f(%9.6f)\n", n, L[n], freqMax, periodT, phase, degrees);
  }

  dataFile = "tex/sunspots_acf.dat"
  if(dataDump)
  {
    sink(dataFile);
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% Sunspot auto-correlation function (ACF) data file suitable for use with LaTeX PStricks\n" );
    printf("%% For an example, see \"%s_acf.tex\"\n", baseName                                  );
    printf("%% This file auto-generated using \"%s\" --- hand-editing not recommended\n", thisfile);
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(i in 1:length(avect))
      printf("  (%12.8f, %12.8f)\n", lvect[i], avect[i]                                         );
    printf("]\n"                                                                                );
    sink();
  }
}

#------------------------------------------------------------------------------
# Main Processing
#------------------------------------------------------------------------------

#  A = matrix( c(1, 2, 3,
#                2, 5, 6,
#                3, 6, 9 ), nrow=3 );
#  Q     = eigen(A, symmetric=TRUE, only.values=FALSE)
#  V     = Q$vectors
#  L     = Q$values
#  D     = diag(L)
#  B = V %*% D %*% inv(V)
#  C = B - A

 spotData = sunspots_getData( dataDump=FALSE, dataPlot=TRUE                  );
 acfData  = sunspots_ACF(     dataDump=FALSE, dataPlot=TRUE,  dataIn=spotData );
#psdData  = sunspots_PSD(     dataDump=FALSE, dataPlot=TRUE,  dataIn=spotData, numSegments=1 );
 psdData  = sunspots_PSD(     dataDump=FALSE, dataPlot=FALSE, dataIn=spotData, numSegments=2 );
 psdData  = sunspots_PSD(     dataDump=FALSE, dataPlot=FALSE, dataIn=spotData, numSegments=3 );
 psdData  = sunspots_PSD(     dataDump=TRUE,  dataPlot=TRUE,  dataIn=spotData, numSegments=4 );
 psdData  = sunspots_PSD(     dataDump=FALSE, dataPlot=FALSE, dataIn=spotData, numSegments=5 );
 psdData  = sunspots_PSD(     dataDump=FALSE, dataPlot=FALSE, dataIn=spotData, numSegments=6 );
 psdData  = sunspots_PSD(     dataDump=FALSE, dataPlot=FALSE, dataIn=spotData, numSegments=7 );
 psdData  = sunspots_PSD(     dataDump=FALSE, dataPlot=FALSE, dataIn=spotData, numSegments=8 );
 psdData  = sunspots_PSD(     dataDump=FALSE, dataPlot=FALSE, dataIn=spotData, numSegments=9 );
 psdData  = sunspots_PSD(     dataDump=FALSE, dataPlot=FALSE, dataIn=spotData, numSegments=10);
# sunspots_PCA_PSD(x)


