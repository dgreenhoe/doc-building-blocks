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
  require(matlib);   # inv
  require(R.utils);  # printf
#data(sunspots, package="datasets"); # deprecated in favor of Silso data
#---------------------------------------
# Global parameters
#---------------------------------------
 author   = "Daniel J. Greenhoe"
 thisfile = "sunspots.R"
 baseName = "sunspots"
 colors = c("blue", "red", "orange", "green", "purple", "brown", "black");
 LaTeXstr = "Sunspot data file suitable for use with LaTeX PStricks / pst-plot"
 AutoGenStr = sprintf("This file auto-generated using \"%s\" --- hand-editing not recommended", thisfile);
#------------------------------------------------------------------------------
# \brief   Estimate Auto-Correlation Function (ACF) of sunspot data minus estimated mean
# \returns data table
#------------------------------------------------------------------------------
sunspots_tseries_data = function( 
  dataDump    = FALSE,
  dataPlot    = TRUE,
  dataFileIn  = "../data/silso_SN_m_tot_V2.0_20210524.csv",
  dataFileBase = "sunspots"
){
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
    sink(sprintf("tex/%s.dat", dataFileBase));
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% %s\n", LaTeXstr );
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
sunspots_tseries_acf = function( 
  dataDump     = FALSE,
  dataPlot     = TRUE,
  nLag         = 2000,
  dataSpots,
  dataFileBase = "sunspots_tseries_acf.dat"
){
  x       = dataSpots;
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
    sink(sprintf("tex/%s.dat", dataFileBase));
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% %s\n", LaTeXstr );
    printf("%% For an example, see \"%.tex\"\n", dataFileBase                                   );
    printf("%% %s\n", AutoGenStr);
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
# \returns Estimated PSD
#------------------------------------------------------------------------------
sunspots_psd_coefs = function( 
  dataDump     = FALSE,
  dataPlot     = TRUE,
  numSegments  = 4,
  dataSpots,
  dataFileBase = "sunspots_psd"
){
  xts       = as.ts(as.vector(dataSpots$count));
  N         = length(xts);         # length of time series
  Fs        = 12;                  # sample rate = 12 samples per year
  estMean   = mean(xts);           # estimated mean
  segLength = N / numSegments;     # segment length
  xpsd      = bspec::welchPSD(xts - estMean, seglength = segLength);
  psdMax    = max(xpsd$power);
  binMax    = ramify::argmax(as.matrix(xpsd$power), rows = FALSE);
  freqMax   = xpsd$frequency[binMax] * Fs; # dominate non-DC frequency
  periodT   = 1 / freqMax;           # estimated period
  resultStr = sprintf("numSegments = %d\nf = %12.8f samples/year\nperiod = %12.8f years\n", numSegments, freqMax, periodT);
  printf("sunspots_PSD(x): %s\n", resultStr);

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
    sink(sprintf("tex/%s.dat",dataFileBase));
    printf("%%=============================================================================\n");
    printf("%% %s \n", author                                                                 );
    printf("%% %s\n", LaTeXstr                                                                );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                );
    printf("%% number of segments          = %d\n",                  numSegments              );
    printf("%% estimated maximum frequency = %12.6f samples/year\n", freqMax                  );
    printf("%% estimated period            = %12.6f years\n",        periodT                  );
    printf("%% %s\n", AutoGenStr                                                              );
    printf("%%=============================================================================\n");
    printf("[\n");
    scaledFreq = xpsd$frequency*Fs;
    scaledGain = sqrt(xpsd$power)/sqrt(psdMax);
    for(i in 1:length(scaledFreq))
      printf("  (%12.8f, %12.8f)\n", scaledFreq[i], scaledGain[i]);
    printf("]\n");
    sink();
  }
  return(xpsd);
}

#------------------------------------------------------------------------------
# \brief Estimate sunspot period using Primary Component Analysis (PCA)
# https://cran.r-project.org/web/packages/matlib/vignettes/inv-ex1.html
# https://stat.ethz.ch/R-manual/R-patched/library/base/html/eigen.html
# https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix
#------------------------------------------------------------------------------
sunspots_eigen_basis = function( 
  dataDump     = FALSE,
  dataPlot     = TRUE,
  Fs           = 12,
  numSegments  = 4,
  nLag         = 2000,
  dataSpots,
  dataFileBase = "sunspots_eigen"
){
  x           = dataSpots;
  xts         = as.ts(as.vector(x$count));
  N           = length(xts);
  estMean     = mean(xts);
  segLength   = N / numSegments;
  a           = acf(xts - estMean, type="correlation", lag=nLag)
  avect       = as.vector(a$acf)
  lvect       = as.vector(a$lag)/12
  A           = stats::toeplitz(avect)
  Q           = eigen(A, symmetric=TRUE, only.values=FALSE)
  V           = Q$vectors
  L           = Q$values
  D           = diag(L)

  if(dataPlot)
  {
    colors = c("blue", "red", "orange", "green", "purple", "black","brown");
    traces = colors[1:6]
    for( n in 1:length(traces))
    {
      if(n==1) plot( L[n] * V[,n], type='l', lwd=3, col=colors[n])
      else     lines(L[n] * V[,n], type='l', lwd=3, col=colors[n])
      traces[n] = sprintf("Eigen Vector %2d", n)
    }
    legend("topleft", legend=traces, col=colors, lwd=3, lty=1:1)
  }

  printf("sunspots_eigen_basis(x) using Welch PSD:\n");
  segLength   = nLag / numSegments
  for( n in 1:10 )
  {
    xpsd    = bspec::welchPSD(x=as.ts(L[n] * Q$vectors[,n]), seglength = segLength);
    psdMax  = max(xpsd$power);
    binMax  = ramify::argmax(as.matrix(xpsd$power), rows = FALSE);
    freqMax = xpsd$frequency[binMax] * Fs; # dominate non-DC frequency
    periodT = 1 / freqMax;                 # estimated period
    printf("Vector %2d (lambda=%10.6f) f=%8.6f samples/year period=%9.6f years\n", n, L[n], freqMax, periodT);
  }

  numYears = 12
  printf("sunspots_eigen_basis(x) using DFT:\n");
  for( n in 1:numYears )
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

  for( numYears in 1:12 )
  {
    sum1 = 0;
    sumLambda = sum(L[1:numYears])
    for( n in 1:numYears )
    {
      x=as.ts(L[n] * Q$vectors[,n])
      xfft    = fft(x, inverse=FALSE);
      fftMax  = abs(xfft);
      binMax  = ramify::argmax(as.matrix(fftMax), rows = FALSE);
      freqMax = (binMax-1) * Fs / length(xfft);
      periodT = 1 / freqMax;                 # estimated period
      phase   = Arg(xfft[binMax]);
      degrees = phase / pi * 180
      sum1 = sum1 + periodT*L[n]
    }
    printf("weighted periodT over %02d years: %12.8f\n",numYears, sum1/sumLambda);
  }

  if(dataDump)
  {
    for(n in 1:8)
    {
      sink(sprintf("tex/%s_%d.dat",dataFileBase,n));
      printf("%%=============================================================================\n"  );
      printf("%% %s \n", author                                                                   );
      printf("%% %s\n", LaTeXstr );
      printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                  );
      printf("%% %s\n", AutoGenStr);
      printf("%%=============================================================================\n"  );
      printf("[\n"                                                                                );
      vect = L[n] * V[,n]
      for(i in 1:length(vect))
        printf("  (%12.8f, %12.8f)\n", i/12, vect[i]                                              );
      printf("]\n"                                                                                );
      sink();
    }
  }
  return(Q);
}

#------------------------------------------------------------------------------
# \brief Data synthesis using eigen vector quasi basis
#------------------------------------------------------------------------------
sunspots_eigen_coefs = function(
  dataDump     = FALSE,                    # dump data to file
  dataPlot     = TRUE,                     # plot data
  Length       = 105,                      # number of coefficients to dump/plot
  dataSpots    = spotData,                 # sunspot data
  dataEigen,                               # eigen-pair data
  dataFileBase = "sunspots_eigen_coefs"    # file base name
){
  V        = dataEigen$vectors             # V = eigen-vectors
  L        = dataEigen$values              # L = eigen-values
  D        = diag(L)                       # D = diagonal matrix of eigen-values
  N        = length(L);                    # N = number of eigen-pairs
  M        = length(dataSpots$count)       # M = number of sunspot values
  spots    = dataSpots$count[(M-N+1):M]    # N most recent sunspot values
  zeroMean = spots - mean(spots)           # zero-mean sunspot data
  stime    = spotData$date[(M-N+1):M]      # N most recent time values
  coefs    = as.numeric( zeroMean %*% V ); # coef[n] = projection of data onto eigen-vector n
  if( dataPlot )
  {
    plot(  coefs[1:Length], col=colors[1], type='h', lwd=5)
    lines( coefs[1:Length], col=colors[1], type='p', lwd=3)
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s.dat",dataFileBase));
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% Sunspot eigen coefficient data (%d coefficients)\n", Length                      );
    printf("%% %s\n", LaTeXstr );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                  );
    printf("%% %s\n", AutoGenStr);
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(n in 1:Length) printf("  (%3d, %14.8f)\n", n-1, coefs[n]                                );
    printf("]\n"                                                                                );
    sink();
  }
  return(coefs)
}

#------------------------------------------------------------------------------
# \brief Data synthesis using eigen vector quasi basis
 #A  = matrix( c(1, 2, 3, 2, 5, 6, 3, 6, 1 ), nrow=3 );
 #B  = V %*% D
 #B2 = V %*% D %*% inv(V)
 #C  = B - A
#------------------------------------------------------------------------------
sunspots_eigen_synth = function(
  dataDump     = FALSE,                     # dump data to file
  dataPlot     = TRUE,                      # plot data
  numCoefs     = 6,                         # number of coefficients to synthesis with
  dataSpots    = spotData,                  # sunspot data
  dataEigen,                                # eigen-pair data
  dataFileBase = "sunspots_eigen_syn"       # file base name
){
  V        = dataEigen$vectors              # eigen-vectors
  L        = dataEigen$values               # eigen-values
  D        = diag(L)                        # diagonal matrix of eigen-values
  N        = length(L);                     # number of eigen-pairs
  M        = length(dataSpots$count)        # number of sunspot data values
  spots    = dataSpots$count[(M-N+1):M]     # last N sunspot data values
  zeroMean = spots - mean(spots)            # zero-mean data
  stime    = spotData$date[(M-N+1):M]       # last N sunspot time values
  coefs    = as.numeric( zeroMean %*% V );  # projection coefficient of spots onto eigen-vector n
  G = sqrt(as.numeric(zeroMean%*%zeroMean)) # estimated energy of sunspot waveform
 #fsyn = coefs %*% V[,1:numCoefs]           # partially-synthesized vector
  fsyn     = 0 * c(1:N)
  g = 0;
  for( n in 1:numCoefs )
  {
    fsyn = fsyn + coefs[n] * V[,n]          # partially-synthesized vector
    g = g + (coefs[n])^2                    # energy of scaled eigen-vectors
  }
  fsyn      = ((G/sqrt(g)) * fsyn)          # scale vector to match energy of original data
  fsyn      = fsyn + mean(spots)            # restore mean
  errorVect = fsyn - spots                  # calculate error vector
  rmsError  = sqrt( (errorVect %*% errorVect))/N # RMS error
  printf("Total RMS synthesis error using %d eigen coefficients = %12.8f\n", numCoefs, rmsError );
  if( dataPlot )
  {
    plot( stime, spots, col=colors[1], type='l')
    lines(stime, fsyn,  col=colors[2], type='l', lwd=3)
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s_sunSpots.dat",dataFileBase));
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% %s\n", LaTeXstr );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                               );
    printf("%% %s\n", AutoGenStr);
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(i in 1:length(fsyn))
      printf("  (%12.8f, %12.8f)\n", stime[i], spots[i]                                          );
    printf("]\n"                                                                                );
    sink();
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s_numCoefs%d.dat",dataFileBase,numCoefs));
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% Sunspot eigen synthesis vector data using %d coefficients\n", numCoefs           );
    printf("%% %s\n", LaTeXstr                                                                  );
    printf("%% For an example, see \"%s_eigen_syn.tex\"\n", baseName                            );
    printf("%% %s\n", AutoGenStr                                                                );
    printf("%% Total RMS synthesis error using %d coefficients = %12.8f\n", numCoefs, rmsError  );
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(n in 1:length(fsyn)) printf("  (%12.8f, %12.8f)\n", stime[n], fsyn[n]                   );
    printf("]\n"                                                                                );
    sink();
  }
  return(fsyn)
}

#------------------------------------------------------------------------------
# \brief ACF of coefficient
#------------------------------------------------------------------------------
sunspots_eigen_acf = function(
  dataDump     = FALSE,                    # dump data to file
  dataPlot     = TRUE,                     # plot data
  Length       = 100,                      # number of ACF elements to dump/plot
  dataCoefs    = coefs,                    # coefficient data
  dataFileBase = "sunspots_coefs_acf"      # file base name
){
  N       = length(dataCoefs);             # number of coefficients provided
  a       = stats::acf(dataCoefs, type="correlation", lag.max=(N-1), plot=FALSE)
  avect   = as.vector(a$acf)               # acf value vector
  lvect   = as.vector(a$lag)               # acf lag vector
  if( dataPlot )
  {
    plot( lvect[1:Length], avect[1:Length], col=colors[1], type='h', lwd=3 );
    lines(lvect[1:Length], avect[1:Length], col=colors[1], type='p', lwd=3 );
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s.dat",dataFileBase));
    printf("%%=============================================================================\n");
    printf("%% %s \n", author                                                                 );
    printf("%% %s\n", LaTeXstr                                                                );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                );
    printf("%% %s\n", AutoGenStr                                                              );
    printf("%%=============================================================================\n");
    printf("[\n"                                                                              );
    for(n in 1:Length) printf("  (%3d, %12.8f)\n", lvect[n], avect[n]                         );
    printf("]\n"                                                                              );
    sink();
  }
  return(avect)
}

#------------------------------------------------------------------------------
# \brief Project data sequence onto sinusoidal basis vectors yielding 
#        a sequence of DFT coefficients
#------------------------------------------------------------------------------
sunspots_dft_coefs = function(
  dataDump     = FALSE,                    # dump data to file
  dataPlot     = TRUE,                     # plot data
  plotLength   = 1000,                     # number of coefficients to dump/plot
  evalLength   = 2000,                     # number of elements to evaluate
  dataSpots    = spotData,                 # sunspot data
  dataFileBase = "sunspots_dft_coefs",     # file base name
  Fs           = 12                        # sample rate in samples per year
){
  M        = length(dataSpots$count)        # number of data values
  N       =  evalLength                     # number data values to evaluate
  spots    = dataSpots$count[(M-N+1):M]     # last N sunspot data values
  zeroMean = (spots - mean(spots))/N        # scaled zero-mean data
  stime    = spotData$date[(M-N+1):M]       # last N sunspot time values
  Xfft     = fft(zeroMean)
  #nLag = evalLength
  #a        = stats::acf(abs(Xfft), type="correlation", lag.max=nLag, plot=FALSE)
  #avect    = as.vector(a$acf)
  #lvect    = as.vector(a$lag)
  f = Fs * c(0:(plotLength-1)) / N
  if( dataPlot )
  {
    plot(  f, Re(Xfft[1:plotLength]), col=colors[1], type='h', lwd=3, xlab="frequency (samples/year)" );
    lines( f, Re(Xfft[1:plotLength]), col=colors[1], type='p', lwd=3 );
    lines( f, Im(Xfft[1:plotLength]), col=colors[2], type='h', lwd=3 );
    lines( f, Im(Xfft[1:plotLength]), col=colors[2], type='p', lwd=3 );
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s_Real.dat",dataFileBase));
    printf("%%=============================================================================\n");
    printf("%% %s \n", author                                                                 );
    printf("%% %s\n", LaTeXstr                                                                );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                );
    printf("%% %s\n", AutoGenStr                                                              );
    printf("%%=============================================================================\n");
    printf("[\n"                                                                              );
    for(n in 1:plotLength) printf("  (%12.8f, %15.8f)\n", f[n], Re(Xfft[n])                   );
    printf("]\n"                                                                              );
    sink();
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s_Imag.dat",dataFileBase));
    printf("%%=============================================================================\n");
    printf("%% %s \n", author                                                                 );
    printf("%% %s\n", LaTeXstr                                                                );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                );
    printf("%% %s\n", AutoGenStr                                                              );
    printf("%%=============================================================================\n");
    printf("[\n"                                                                              );
    for(n in 1:plotLength) printf("  (%12.8f, %15.8f)\n", f[n], Im(Xfft[n])                   );
    printf("]\n"                                                                              );
    sink();
  }
  #if( dataDump )
  #{
  #  sink(sprintf("tex/%s_acf.dat",dataFileBase));
  #  printf("%%=============================================================================\n");
  #  printf("%% %s \n", author                                                                 );
  #  printf("%% %s\n", LaTeXstr                                                                );
  #  printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                );
  #  printf("%% %s\n", AutoGenStr                                                              );
  #  printf("%%=============================================================================\n");
  #  printf("[\n"                                                                              );
  #  for(n in 1:Length) printf("  (%12.8f, %12.8f)\n", lvect[n], avect[n]                      );
  #  printf("]\n"                                                                              );
  #  sink();
  #}
  return(Xfft)
}

#------------------------------------------------------------------------------
# \brief Data synthesis using DFT basis
#------------------------------------------------------------------------------
sunspots_dft_syn = function(
  dataDump     = FALSE,                     # dump data to file
  dataPlot     = TRUE,                      # plot data
  numCoefs     = 6,                         # number of coefficients to synthesis with
  dataSpots    = spotData,                  # sunspot data
  dataDFT,                                  # DFT coefficients
  plotLength   = 1000,                      # number of coefficients to dump/plot
  evalLength   = 2000,                      # number of elements to evaluate
  Fs           = 12,                        # sample rate in samples per year
  dataFileBase = "sunspots_dft_syn"         # file base name
){
  V = cos(2*pi*c(0:(numCoefs-1)));
  V        = dataEigen$vectors              # eigen-vectors
  L        = dataEigen$values               # eigen-values
  D        = diag(L)                        # diagonal matrix of eigen-values
  N        = length(L);                     # number of eigen-pairs
  M        = length(dataSpots$count)        # number of sunspot data values
  spots    = dataSpots$count[(M-N+1):M]     # last N sunspot data values
  zeroMean = spots - mean(spots)            # zero-mean data
  stime    = spotData$date[(M-N+1):M]       # last N sunspot time values
  coefs    = as.numeric( zeroMean %*% V );  # projection coefficient of spots onto eigen-vector n
  G = sqrt(as.numeric(zeroMean%*%zeroMean)) # estimated energy of sunspot waveform
  fsyn     = 0 * c(1:N)
  g = 0;
  for( n in 1:numCoefs )
  {
    fsyn = fsyn + coefs[n] * V[,n]          # partially-synthesized vector
    g = g + (coefs[n])^2                    # energy of scaled eigen-vectors
  }
  fsyn      = ((G/sqrt(g)) * fsyn)          # scale vector to match energy of original data
  fsyn      = fsyn + mean(spots)            # restore mean
  errorVect = fsyn - spots                  # calculate error vector
  rmsError  = sqrt( (errorVect %*% errorVect))/N # RMS error
  printf("Total RMS synthesis error using %d eigen coefficients = %12.8f\n", numCoefs, rmsError );
  if( dataPlot )
  {
    plot( stime, spots, col=colors[1], type='l')
    lines(stime, fsyn,  col=colors[2], type='l', lwd=3)
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s_sunSpots.dat",dataFileBase));
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% %s\n", LaTeXstr );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                               );
    printf("%% %s\n", AutoGenStr);
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(i in 1:length(fsyn))
      printf("  (%12.8f, %12.8f)\n", stime[i], spots[i]                                          );
    printf("]\n"                                                                                );
    sink();
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s_numCoefs%d.dat",dataFileBase,numCoefs));
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% Sunspot eigen synthesis vector data using %d coefficients\n", numCoefs           );
    printf("%% %s\n", LaTeXstr                                                                  );
    printf("%% For an example, see \"%s_eigen_syn.tex\"\n", baseName                            );
    printf("%% %s\n", AutoGenStr                                                                );
    printf("%% Total RMS synthesis error using %d coefficients = %12.8f\n", numCoefs, rmsError  );
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(n in 1:length(fsyn)) printf("  (%12.8f, %12.8f)\n", stime[n], fsyn[n]                   );
    printf("]\n"                                                                                );
    sink();
  }
  return(fsyn)
}

#------------------------------------------------------------------------------
# Main Processing
#------------------------------------------------------------------------------
 T = TRUE
 F = FALSE
 spotData  = sunspots_tseries_data( dataDump=F, dataPlot=T                                    );
 acfData   = sunspots_tseries_acf(  dataDump=F, dataPlot=T, dataSpots=spotData                );
 psdCoefs  = sunspots_psd_coefs(    dataDump=F, dataPlot=T, dataSpots=spotData, numSegments=4 );
 eigenPairs= sunspots_eigen_basis(  dataDump=F, dataPlot=T, dataSpots=spotData, nLag=2000     );
 coefs     = sunspots_eigen_coefs(  dataDump=F, dataPlot=F, dataSpots=spotData, dataEigen=eigenPairs, Length=105 );
 fsyn      = sunspots_eigen_synth(  dataDump=F, dataPlot=F, dataSpots=spotData, dataEigen=eigenPairs, numCoefs=105   );
 eigenACF  = sunspots_eigen_acf(    dataDump=F, dataPlot=F, dataCoefs=coefs,    Length=100 );
 dftACF    = sunspots_dft_coefs(    dataDump=F, dataPlot=F, dataSpots=spotData, evalLength=2000, plotLength=1001 );
