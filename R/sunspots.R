#============================================================================
# Daniel J. Greenhoe
# R script file
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
#install.packages("Hadamard.R")
#install.packages("pracma")
#install.packages("timma")
#install.packages('bit')
#install.packages('bitops')

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
  require(HadamardR)
  require(pracma)
  require(timma)   # grays(n): Generate gray code vector of length 2^n
  require(bitops)  # bitAnd, bitOr, bitXor, bitFlip, bitShiftL

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
# \brief   Get sunspot data
# \returns sunspot data
# data(sunspots, package="datasets"); # deprecated in favor of Silso data
#------------------------------------------------------------------------------
sunspots_tseries_data = function(
  verbose     = TRUE,
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
  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
    printf("  Data from \"%s\"\n", dataFileIn );
  }
  if( dataPlot )
  {
    plot(x=dvect, y=cvect, lwd=2, type='l', main="Sunspot Estimation", col="blue", xlab="year", ylab="count");
    abline(h=seq(from=0,   to= 400, by=100), lty='dashed', col = "green")
    abline(v=seq(from=1750,to=2020, by= 10), lty='dashed', col = "green")
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s.dat", dataFileBase ));
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
# \brief Calculate sunspot sample rate in samples/year
#------------------------------------------------------------------------------
calc_Fs = function( 
  verbose = FALSE, 
  dataIn = spotData 
){
  numSamples = length(spotData$year);
  numYears   = max(spotData$year)-min(spotData$year)+1;
  sampleRate = numSamples / numYears;
  Fs         = round( sampleRate );
  error      = (Fs - sampleRate) / sampleRate;
  if( verbose )
  {
    printf("calc_Fs(...)\n");
    printf("  numPoints  = %d\n"  , numSamples );
    printf("  numYears   = %d\n"  , numYears   );
    printf("  sampleRate = %f\n"  , sampleRate );
    printf("  Fs         = %d\n"  , Fs         );
    printf("  error      = %f%%\n", error*100  );
  }
  return( Fs )
}

#------------------------------------------------------------------------------
# \brief   Grid
#------------------------------------------------------------------------------
abgrid = function(xmin, xmax, ymin, ymax, xstep, ystep)
{
  y = seq( from=ymin, to=ymax, by=ystep );
  x = seq( from=xmin, to=xmax, by=xstep );
  abline( h=y, lty='dashed', col = "green");
  abline( v=x, lty='dashed', col = "green");
  L = list(x=x, y=y);
  return( L );
}

#------------------------------------------------------------------------------
# \brief   Estimate Auto-Correlation Function (ACF) of sunspot data minus estimated mean
# \returns ACF data table
#------------------------------------------------------------------------------
sunspots_tseries_acf = function(
  verbose      = TRUE,
  dataDump     = FALSE,
  dataPlot     = TRUE,
  nLag         = 2000,
  dataIn       = spotData,                 # sunspot data
  dataFileBase = "sunspots_tseries_acf"
){
  x       = dataIn;
  xts     = as.ts(as.vector(x$count));
  estMean = mean(xts);
  a       = stats::acf(xts - estMean, type="correlation", lag.max=nLag, plot=FALSE)
  avect   = as.vector(a$acf)
  lvect   = as.vector(a$lag)/12

  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
  }
  if( dataPlot )
  {
    plot(x=lvect, y=avect, lwd=2, type='l', xaxp=c(0,160,16), yaxp=c(-0.4,1.0,14), col="blue",
         main="Sunspot Auto-Correlation Estimation", xlab="lag in years", ylab="magnitude");
    abgrid( xmin=0, xmax=170, ymin=-0.4, ymax=1.0, xstep=10, ystep=0.1 );
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s.dat", dataFileBase ));
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
  verbose      = TRUE,
  dataDump     = FALSE,
  dataPlot     = TRUE,
  numSegments  = 4,
  dataIn       = spotData,
  dataFileBase = "sunspots_psd"
){
  Fs        = calc_Fs( dataIn = spotData);    # sample rate in samples/year
  xts       = as.ts(as.vector(dataIn$count)); #
  N         = length(xts);                    # length of time series
  estMean   = mean(xts);                      # estimated mean
  segLength = N / numSegments;                # segment length
  xpsd      = bspec::welchPSD(xts - estMean, seglength = segLength);
  psdMax    = max(xpsd$power);
  binMax    = ramify::argmax(as.matrix(xpsd$power), rows = FALSE);
  freqMax   = xpsd$frequency[binMax] * Fs; # dominate non-DC frequency
  periodT   = 1 / freqMax;           # estimated period
  resultStr = sprintf("numSegments = %d\nf = %12.8f samples/year\nperiod = %12.8f years\n", numSegments, freqMax, periodT);
  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
  }
  if( dataPlot )
  {
    titleStr = sprintf("0-average Sunspot Power Spectral Density (PSD) Estimate using welchPSD with numSegments=%d",numSegments);
    plot(x=xpsd$frequency*Fs, y=sqrt(xpsd$power/N), type="h", lwd=3, col="blue", xaxp=c(0,Fs/2,60), yaxp=c(0,20,5), xlab="", ylab="", main="", sub="");
    title(main=titleStr, xlab="samples/year (max=Fs/2=6 samples/year)", ylab="Gain (sqrt of power)");
    title(main=resultStr, line=-10);
    abgrid( xmin=0, xmax=Fs/2, xstep=0.1, ymin=0, ymax=20, ystep=2 );
  }
  if( dataDump )
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
  verbose      = TRUE,
  dataDump     = FALSE,
  dataPlot     = TRUE,
  numSegments  = 4,
  windowLength = 2048,
  dataIn       = spotData,
  dataFileBase = "sunspots_eigen"
){
  Fs          = calc_Fs( dataIn = spotData);
  nLag        = windowLength - 1;
  x           = dataIn;
  xts         = as.ts(as.vector(x$count));
  N           = length(xts);
  estMean     = mean(xts);
  segLength   = N / numSegments;
  a           = stats::acf(xts - estMean, type="correlation", lag=nLag, plot=FALSE)
  avect       = as.vector(a$acf)
  lvect       = as.vector(a$lag)/12
  A           = stats::toeplitz(avect)
  Q           = eigen(A, symmetric=TRUE, only.values=FALSE)
  V           = Q$vectors
  L           = Q$values
  D           = diag(L)

  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
  }
  if( dataPlot )
  {
    traces = colors[1:6]
    for( n in 1:length(traces))
    {
      if(n==1) plot( L[n] * V[,n], type='l', lwd=3, col=colors[n])
      else     lines(L[n] * V[,n], type='l', lwd=3, col=colors[n])
      traces[n] = sprintf("Eigen Vector %2d", n)
    }
    legend("topleft", legend=traces, col=colors, lwd=3, lty=1:1)
    abgrid( xmin=0, xmax=nLag, xstep=100, ymin=-8, ymax=8, ystep=2 );
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

  if( dataDump )
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
  verbose      = TRUE,
  dataDump     = FALSE,                    # dump data to file
  dataPlot     = TRUE,                     # plot data
  Length       = 105,                      # number of coefficients to dump/plot
  dataIn       = spotData,                 # sunspot data
  basis,                                   # eigen-pair data
  dataFileBase = "sunspots_eigen_coefs"    # file base name
){
  V        = basis$vectors                 # V = eigen-vectors
  L        = basis$values                  # L = eigen-values
  D        = diag(L)                       # D = diagonal matrix of eigen-values
  N        = length(L);                    # N = number of eigen-pairs
  M        = length(dataIn$count)          # M = number of sunspot values
  spots    = dataIn$count[(M-N+1):M]       # N most recent sunspot values
  zeroMean = spots - mean(spots)           # zero-mean sunspot data
  stime    = spotData$date[(M-N+1):M]      # N most recent time values
  coefs    = as.numeric( zeroMean %*% V ); # coef[n] = projection of data onto eigen-vector n
  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
  }
  if( dataPlot )
  {
    plot(  coefs[1:Length], col=colors[1], type='h', lwd=5)
    lines( coefs[1:Length], col=colors[1], type='p', lwd=3)
    abgrid( xmin=0, xmax=120, xstep=5, ymin=-700, ymax=2000, ystep=100 );
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
  verbose      = TRUE,
  dataDump     = FALSE,                     # dump data to file
  dataPlot     = TRUE,                      # plot data
  numCoefs     = 6,                         # number of coefficients to synthesis with
  dataIn       = spotData,                  # sunspot data
  basis,                                    # eigen-pair data
  dataFileBase = "sunspots_eigen_synth"     # file base name
){
  V        = basis$vectors                  # eigen-vectors
  L        = basis$values                   # eigen-values
  D        = diag(L)                        # diagonal matrix of eigen-values
  N        = length(L);                     # number of eigen-pairs
  M        = length(dataIn$count)           # number of sunspot data values
  spots    = dataIn$count[(M-N+1):M]        # last N sunspot data values
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
  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
    printf("  Total RMS synthesis error using %d eigen coefficients = %12.8f\n", numCoefs, rmsError );
    printf("  Energy of zero mean vector         G = %12.8f\n", G)
    printf("  Energy of synthesized vector sqrt(g) = %12.8f\n", sqrt(g))
    printf("  Energy ratio G/sqrt(g)               = %12.8f\n", G/sqrt(g))
  }
  if( dataPlot )
  {
    plot( stime, spots, col=colors[1], type='l')
    lines(stime, fsyn,  col=colors[2], type='l', lwd=3)
    abgrid( xmin=1850, xmax=2020, xstep=10, ymin=0, ymax=350, ystep=50 );
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
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                  );
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
  verbose      = TRUE,
  dataDump     = FALSE,                    # dump data to file
  dataPlot     = TRUE,                     # plot data
  Length       = 100,                      # number of ACF elements to dump/plot
  coefs        = coefs,                    # coefficient data
  dataFileBase = "sunspots_coefs_acf"      # file base name
){
  N       = length(coefs);                 # number of coefficients provided
  a       = stats::acf(coefs, type="correlation", lag.max=(N-1), plot=FALSE)
  avect   = as.vector(a$acf)               # acf value vector
  lvect   = as.vector(a$lag)               # acf lag vector
  if( dataPlot )
  {
    plot( lvect[1:Length], avect[1:Length], col=colors[1], type='h', lwd=3 );
    lines(lvect[1:Length], avect[1:Length], col=colors[1], type='p', lwd=3 );
    abgrid( xmin=0, xmax=Length, xstep=Length/10, ymin=-0.2, ymax=1, ystep=0.1 );
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
# \brief Calculate Discrete Fourier basis
#------------------------------------------------------------------------------
sunspots_dft_basis = function(
  verbose      = TRUE,                    # verbose
  dataDump     = FALSE,                   # dump data to file
  dataPlot     = TRUE,                    # plot 
  Fs           = 12,                      # samples / year
  windowLength = 2048,                    # evaluation window length
  numVectors   = 5,                       # plot vectors 0 -- (numVectors-1)
  dataIn       = spotData,                # sunspot data
  dataFileBase = "sunspots_dft_basis"     # base file name
){                                        #
  N = windowLength                        # N = length of data/basis vectors
  M = ceiling(N/2)  #(N-1)/2              # M = number of basis vectors
  n = matrix(seq(from=0, by=1, length=N));# 
  k = matrix(seq(from=0, by=1, length=M));# 
  f = k/N*Fs                              # 
  V = cos(2*pi*(n %*% t(k))/N);           # 
  W = sin(2*pi*(n %*% t(k))/N);           # 

  if( dataPlot )
  {
    traces = colors[1:numVectors]
    for( n in 1:length(traces))
    {
      v = V[c(1:M), n]
      if(n==1) plot( f, v, type='p', lwd=2, col=colors[n], ylim=c(-1.2,1.2))
      else     lines(f, v, type='p', lwd=2, col=colors[n])
      traces[n] = sprintf("DFT basis vector k=%2d", n-1)
    }
    abgrid( xmin=0, xmax=Fs/2, xstep=0.2, ymin=-1, ymax=1, ystep=0.2 );
    legend("topleft", legend=traces, col=colors, lwd=3, lty=1:1)
  }
  if( dataDump )
  {
    for(n in 1:numVectors)
    {
      vect = V[,n]
      sink(sprintf("tex/%s_cos_%d.dat",dataFileBase,n-1));
      printf("%%=============================================================================\n");
      printf("%% %s \n", author                                                                 );
      printf("%% %s\n", LaTeXstr                                                                );
      printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                );
      printf("%% Note: The Euclidean norm of this vector is sqrt(N/2)=sqrt(%d/2)=%.8f\n", N, sqrt(N/2) );
      printf("%%       To normalize, scale by 1/sqrt(N/2)=%.8f\n", 1/sqrt(N/2)                  );
      printf("%% %s\n", AutoGenStr                                                              );
      printf("%%=============================================================================\n");
      printf("[\n"                                                                              );
      for(i in 1:length(vect))
        printf("  (%12.8f, %12.8f)\n", f[i], vect[i]                                            );
      printf("]\n"                                                                              );
      sink();
    }
    for(n in 1:numVectors)
    {
      vect = W[,n]
      sink(sprintf("tex/%s_sin_%d.dat",dataFileBase,n-1));
      printf("%%=============================================================================\n");
      printf("%% %s \n", author                                                                 );
      printf("%% %s\n", LaTeXstr                                                                );
      printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                );
      printf("%% Note: The Euclidean norm of this vector is sqrt(N/2)=sqrt(%d/2)=%.8f\n", N, sqrt(N/2) );
      printf("%%       To normalize, scale by 1/sqrt(N/2)=%.8f\n", 1/sqrt(N/2)                  );
      printf("%% %s\n", AutoGenStr                                                              );
      printf("%%=============================================================================\n");
      printf("[\n"                                                                              );
      for(i in 1:length(vect))
        printf("  (%12.8f, %12.8f)\n", f[i], vect[i]                                            );
      printf("]\n"                                                                              );
      sink();
    }
  }
  L = list( f=f, V=V, W=W, N=N, M=M );
  return( L );
}

#------------------------------------------------------------------------------
# \brief Project data sequence onto sinusoidal basis vectors yielding
#        a sequence of DFT coefficients
#------------------------------------------------------------------------------
sunspots_dft_coefs = function(
  verbose      = TRUE,
  dataDump     = FALSE,                    # dump data to file
  dataPlot     = TRUE,                     # plot data
  basis        = dftBasis,                 # basis vectors
  plotLength   = 1024,                     # number of coefficients to dump/plot
  dataIn       = spotData,                 # sunspot data
  dataFileBase = "sunspots_dft_coefs",     # file base name
  Fs           = 12                        # sample rate in samples per year
){
  Fs       = calc_Fs( dataIn = spotData);  # sample rate in samples/year
  N        = length(basis$V[,1])
  M        = length(dataIn$count)          # number of data values
  spots    = dataIn$count[(M-N+1):M]       # last N sunspot data values
  estMean  = mean(spots)                   # estimated mean
  zeroMean = (spots - estMean)/N           # scaled zero-mean data
  stime    = spotData$date[(M-N+1):M]      # last N sunspot time values
  V        = basis$V
  W        = basis$W
  coefsV   = as.numeric( zeroMean %*% V );    # coef[n] = projection of data onto basis-vector n
  coefsW   = as.numeric( zeroMean %*% (-W) ); # coef[n] = projection of data onto basis-vector n
  f = Fs * c(0:(plotLength-1)) / N
  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
    printf("  Window Length = %d\n", N);
    printf("  Number of Coefficients = %d\n", plotLength);
  }
  if( dataPlot )
  {
    plot(  f, coefsV, col=colors[1], type='h', lwd=3, xlab="frequency (samples/year)" );
    lines( f, coefsV, col=colors[1], type='p', lwd=3 );
    lines( f, coefsW, col=colors[2], type='h', lwd=3 );
    lines( f, coefsW, col=colors[2], type='p', lwd=3 );
    abgrid( xmin=0, xmax=Fs/2, xstep=0.1, ymin=-25, ymax=25, ystep=5 );
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
    for(n in 1:plotLength) printf("  (%12.8f, %15.8f)\n", f[n], coefsV[n]                     );
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
    for(n in 1:plotLength) printf("  (%12.8f, %15.8f)\n", f[n], coefsW[n]                     );
    printf("]\n"                                                                              );
    sink();
  }
  L = list( f=f, V=coefsV, W=coefsW );
  return( L );
}

#------------------------------------------------------------------------------
# \brief Data synthesis using DFT basis
#------------------------------------------------------------------------------
sunspots_dft_synth = function(
  verbose      = TRUE,
  dataDump     = FALSE,                     # dump data to file
  dataPlot     = TRUE,                      # plot data
  numCoefs     = 17,                        # number of coefficient pairs to synthesize with
  dataIn       = spotData,                  # sunspot data
  basis,                                    # DFT coefficients
  dataFileBase = "sunspots_dft_synth"       # file base name
){
  Fs       = calc_Fs( dataIn = spotData);   # sample rate in samples/year
  N        = length(basis$V[,1])            # N = length of single basis vector
  M        = length(dataIn$count)           # number of sunspot data values
  V        = basis$V / sqrt(N/2)            # normalize V basis vectors such that ||V[,n]||=1
  W        = basis$W / sqrt(N/2)            # normalize W basis vectors such that ||W[,n]||=1
  V[,1]    = V[,1]/sqrt(2)                  # special case V[,1]: make ||V[,1]||=1
  spots    = dataIn$count[(M-N+1):M]        # last N sunspot data values
  estMean  = mean(spots)                    # estimated mean
  zeroMean = spots - estMean                # zero-mean data
  stime    = spotData$date[(M-N+1):M]       # last N sunspot time values
  coefsV   = as.numeric( zeroMean %*% V );  # projection coefficient of spots onto real basis
  coefsW   = as.numeric( zeroMean %*% W );  # projection coefficient of spots onto imag basis
  G = sqrt(as.numeric(zeroMean%*%zeroMean)) # estimated energy of sunspot waveform
  f1       = 0 * c(1:N)
  f2       = 0 * c(1:N)
  fsyn     = 0 * c(1:N)
  g1       = 0;
  g2       = 0;
  for( n in 1:numCoefs )
  {
    f1 = f1 + (coefsV[n] * V[,n])
    f2 = f2 + (coefsW[n] * W[,n])
    g1 = g1 + (coefsV[n]^2)                   # energy of scaled V vectors
    g2 = g2 + (coefsW[n]^2)                   # energy of scaled W vectors
  }
  g=g1+g2
  fsyn      = f1 + f2                         # add synthesis due to cos and sin bases
  if(g>1e-9) fsyn      = ((G/sqrt(g)) * fsyn) # scale vector to match energy of original data
  fsyn      = fsyn + estMean                  # restore mean
  errorVect = fsyn - spots                    # calculate error vector
  rmsError  = sqrt( (errorVect %*% errorVect))/N # RMS error
  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
    printf("  Total RMS synthesis error using 2 x %d = %d dft coefficients = %12.8f\n", numCoefs, 2*numCoefs, rmsError );
    printf("  Energy of zero mean vector         G = %12.8f\n", G)
    printf("  Energy of synthesized vector sqrt(g) = %12.8f\n", sqrt(g))
    printf("  Energy ratio G/sqrt(g)               = %12.8f\n", G/sqrt(g))
  }
  if( dataPlot )
  {
    plot( stime, spots, col=colors[1], type='l')
    lines(stime, fsyn,  col=colors[2], type='l', lwd=3)
    abgrid( xmin=1850, xmax=2020, xstep=10, ymin=0, ymax=350, ystep=50 );
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s_sunSpots.dat",dataFileBase));
    printf("%%=============================================================================\n");
    printf("%% %s \n", author                                                                 );
    printf("%% %s\n", LaTeXstr                                                                );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                );
    printf("%% %s\n", AutoGenStr                                                              );
    printf("%%=============================================================================\n");
    printf("[\n"                                                                              );
    for(i in 1:length(fsyn))
      printf("  (%12.8f, %12.8f)\n", stime[i], spots[i]                                       );
    printf("]\n"                                                                              );
    sink();
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s_numCoefs%d.dat",dataFileBase,2*numCoefs));
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% Sunspot DFT synthesis vector data using 2 x %d = %d coefficients\n", numCoefs, 2*numCoefs);
    printf("%% %s\n", LaTeXstr                                                                  );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                  );
    printf("%% %s\n", AutoGenStr                                                                );
    printf("%% Total RMS synthesis error using %d coefficient pairs = %12.8f\n", numCoefs, rmsError);
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(n in 1:length(fsyn)) printf("  (%12.8f, %12.8f)\n", stime[n], fsyn[n]                   );
    printf("]\n"                                                                                );
    sink();
  }
  L = list( N=N, M=M, coefsV=coefsV, coefsW=coefsW, V=V, W=W, fsyn=fsyn, z=zeroMean, estMean=estMean, G=G, g=g, g1=g1, g2=g2 );
  return( L );
}

#------------------------------------------------------------------------------
# \brief ACF of coefficient
#------------------------------------------------------------------------------
acfComplex = function(z, dataPlot=FALSE)
{
  N = length(z);
  x = Re(z);
  y = Im(z);
  A = stats::ccf(x, x, type="correlation", lag.max=(N-1), plot=FALSE)
  B = stats::ccf(y, y, type="correlation", lag.max=(N-1), plot=FALSE)
  C = stats::ccf(x, y, type="correlation", lag.max=(N-1), plot=FALSE)
  D = stats::ccf(y, x, type="correlation", lag.max=(N-1), plot=FALSE)
  a = A$acf[N:(2*N-1)];               # acf value vector
  b = B$acf[N:(2*N-1)];               # acf value vector
  c = C$acf[N:(2*N-1)];               # acf value vector
  d = D$acf[N:(2*N-1)];               # acf value vector
  acfReal = a + b;
  acfImag = d - c;
  acfMag = 0.5*sqrt(acfReal^2 + acfImag^2)
  avect   = A$acf[N:(2*N-1)];               # acf value vector
  lvect   = A$lag[N:(2*N-1)];               # acf lag vector
  if( dataPlot )
  {
    plot( lvect[1:Length], acfMag[1:Length], col=colors[1], type='h', lwd=3 );
    lines(lvect[1:Length], acfMag[1:Length], col=colors[1], type='p', lwd=3 );
   #abgrid( xmin=0, xmax=Fs/2, xstep=0.1, ymin=-25, ymax=25, ystep=5 );
  }
  L = list(N=N, a=a, b=b, c=c, d=d, mag=acfMag, lag=lvect);
  return( L );
}

#------------------------------------------------------------------------------
# \brief ACF of coefficient
#------------------------------------------------------------------------------
sunspots_dft_acf = function(
  verbose      = TRUE,
  dataDump     = FALSE,                    # dump data to file
  dataPlot     = TRUE,                     # plot data
  Length       = 100,                      # number of ACF elements to dump/plot
  coefs        = dftCoefs,                 # coefficient data
  dataFileBase = "sunspots_dft_acf"        # file base name
){
  V = coefs$V
  W = coefs$W
  dataComplex = complex( real=V, imaginary=W )
  N       = length(dataComplex);             # number of coefficients provided
  acfData = acfComplex(dataComplex);
  acfMag  = acfData$mag;
  lvect   = acfData$lag;
  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
    printf("  Number of coefficient pairs = %d\n", N);
  }
  if( dataPlot )
  {
    plot( lvect[1:Length], acfMag[1:Length], col=colors[1], type='h', lwd=3 );
    lines(lvect[1:Length], acfMag[1:Length], col=colors[1], type='p', lwd=3 );
    abgrid( xmin=0, xmax=Length, xstep=Length/10, ymin=0, ymax=1, ystep=0.1 );
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
    for(n in 1:Length) printf("  (%3d, %12.8f)\n", lvect[n], acfMag[n]                        );
    printf("]\n"                                                                              );
    sink();
  }
  L = list( lvect=lvect, acfMag=acfMag );
  return( acfMag );
}

#------------------------------------------------------------------------------
# \brief Bit reverse n-bit binary number |b_{n-1}|b_{n-2}|...|b_1|b_0|
# \returns binary number |b_0|b_1|...|b_{n-2}|b_{n-1}|
# \example: ncols=3
#   bitrev( c(0:7), 3 ) = [ 0 4 2 6 1 5 3 7 ]
#------------------------------------------------------------------------------
bitrev = function( avect, nBits )
{
  N = length( avect );
  bvect = 0 * c(1:N)
  for( n in 1:N )
  {
    a = avect[n]
    b = 0
    maska = 1
    maskb = bitwShiftL( 1, nBits-1 )
    for( i in 1:nBits )
    {
      if( bitwAnd( a, maska ) ) { b = bitwOr( b, maskb ) }
      maska = bitwShiftL( maska, 1 );
      maskb = bitwShiftR( maskb, 1 );
    }
    bvect[n] = b;
  }
  return( bvect )
}

#------------------------------------------------------------------------------
# \brief Generate ncols x ncols Bit-Reversal Matrix in which each column n
#        has row=bitrev(n-1,ncols)+1 set equal to 1 
#        and all other row positions set equal to 0.
# \example ncols=3
#   B = bitReverseMatrix(8)
#        _ col1 col2 col3 col4 col5 col6 col7 col8 _
#       |    1    0    0    0    0    0    0    0   | row1
#       |    0    0    0    0    1    0    0    0   | row2
#       |    0    0    1    0    0    0    0    0   | row3
#   B = |    0    0    0    0    0    0    1    0   | row4
#       |    0    1    0    0    0    0    0    0   | row5
#       |    0    0    0    0    0    1    0    0   | row6
#       |    0    0    0    1    0    0    0    0   | row7
#       |_   0    0    0    0    0    0    0    1  _| row8
#
#------------------------------------------------------------------------------
bitReverseMatrix = function( ncols )
{
  B = matrix( seq(from=0, to=0, length=ncols^2), ncol=ncols );
  numbits = log(ncols)/log(2)
  for( n in 1:ncols )
  {
    m = bitrev(n-1, numbits) + 1
    B[m,n] = 1
  }
  return( B )
}

#------------------------------------------------------------------------------
# \brief Generate ncols x ncols Gray-Code Matrix in which each column n
#        has row=grays(ncols)[n]+1 set equal to 1 
#        and all other row positions set equal to 0.
# \example ncols=3
#   grays( ncols )     = [ 0 1 3 2 6 7 5 4 ]
#   grays( ncols ) + 1 = [ 1 2 4 3 7 8 6 5 ]
#   G = grayCodeMatrix(8)
#        _ col1 col2 col3 col4 col5 col6 col7 col8 _
#       |    1    0    0    0    0    0    0    0   | row1
#       |    0    1    0    0    0    0    0    0   | row2
#       |    0    0    0    1    0    0    0    0   | row3
#   G = |    0    0    1    0    0    0    0    0   | row4
#       |    0    0    0    0    0    0    0    1   | row5
#       |    0    0    0    0    0    0    1    0   | row6
#       |    0    0    0    0    1    0    0    0   | row7
#       |_   0    0    0    0    0    1    0    0  _| row8
#------------------------------------------------------------------------------
grayCodeMatrix = function( ncols )
{
  G = matrix( seq(from=0, to=0, length=ncols^2), ncol=ncols );
  numbits = log(ncols)/log(2)
  grayCodes = grays( numbits )
  for( n in 1:ncols )
  {
    m = grayCodes[n] + 1
    G[m,n] = 1
  }
  return( G )
}

#------------------------------------------------------------------------------
# \brief Generate ncols x ncols Walsh Matrix in sequency ordering
# https://en.wikipedia.org/wiki/Walsh_matrix#Sequency_ordering
# Tam and Goulet (1992) "On Arithmetical Shift for Walsh Functions"
#   https://ieeexplore.ieee.org/document/1672117
# \example
#   W =  Walsh_seq_Matrix( 8 )
#       _ col1 col2 col3 col4 col5 col6 col7 col8 _
#      |    1    1    1    1    1    1    1    1   | row1
#      |    1    1    1    1   -1   -1   -1   -1   | row2
#      |    1    1   -1   -1   -1   -1    1    1   | row3
#   W= |    1    1   -1   -1    1    1   -1   -1   | row4
#      |    1   -1   -1    1    1   -1   -1    1   | row5
#      |    1   -1   -1    1   -1    1    1   -1   | row6
#      |    1   -1    1   -1   -1    1   -1    1   | row7
#      |_   1   -1    1   -1    1   -1    1   -1  _| row8
#       --> increasing number of sign changes -->
#------------------------------------------------------------------------------
Walsh_seq_Matrix = function( ncols )
{
  H = hadamard(ncols)            # generate Hadamard matrix
  B = bitReverseMatrix( ncols )  # generate bit-reversal matrix
  G = grayCodeMatrix( ncols )    # generate Gray-Code matrix
  W = (H %*% B) %*% G            # generate Walsh matrix
  return(W)
}

#------------------------------------------------------------------------------
# \brief Calculate Walsh basis
#------------------------------------------------------------------------------
sunspots_walsh_basis = function(
  verbose      = TRUE,
  dataDump     = FALSE,
  dataPlot     = TRUE,
  numVectors   = 5,
  windowLength = 8,
  dataIn       = spotData,                 # sunspot data
  dataFileBase = "sunspots_walsh_basis"
){
  Fs = calc_Fs( dataIn = spotData);   # sample rate in samples/year
  N  = windowLength
  f  = seq(from=0, by=(1/Fs), length=N)
  W  = Walsh_seq_Matrix( N )

  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
  }
  if( dataPlot )
  {
    traces = colors[1:numVectors]
    for( n in 1:length(traces))
    {
      if(n==1) plot(  f, W[,n], type='o', lwd=2, col=colors[n], ylim=c(-1.2,1.2))
      else     lines( f, W[,n], type='o', lwd=2, col=colors[n])
      traces[n] = sprintf("Walsh basis vector k=%2d", n-1)
    }
    abgrid( xmin=0, xmax=Fs/2, xstep=0.2, ymin=-1, ymax=1, ystep=0.2 );
    legend("topleft", legend=traces, col=colors, lwd=3, lty=1:1)
  }
  if( dataDump )
  {
    for(n in 1:numVectors)
    {
      vect = W[,n]
      sink(sprintf("tex/%s_%d.dat",dataFileBase,n-1));
      printf("%%=============================================================================\n");
      printf("%% %s \n", author                                                                 );
      printf("%% %s\n", LaTeXstr                                                                );
      printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                );
      printf("%% %s\n", AutoGenStr                                                              );
      printf("%%=============================================================================\n");
      printf("[\n"                                                                              );
      for(i in 1:length(vect))
        printf("  (%12.8f, %12.8f)\n", f[i], vect[i]                                            );
      printf("]\n"                                                                              );
      sink();
    }
  }
  L = list(f=f, W=W)
  return( L );
}

#------------------------------------------------------------------------------
# \brief Data synthesis using eigen vector quasi basis
#------------------------------------------------------------------------------
sunspots_walsh_coefs = function(
  verbose      = TRUE,
  dataDump     = FALSE,                    # dump data to file
  dataPlot     = TRUE,                     # plot data
  plotLength   = 105,                      # number of coefficients to plot
  dataIn       = spotData,                 # sunspot data
  basis,                                   # basis data
  dataFileBase = "sunspots_walsh_coefs"    # file base name
){
  V        = basis$W                       # V = eigen-vectors
  N        = length(V[,1]);                # N = length of each basis vector
  M        = length(dataIn$count)          # M = number of sunspot values
  V        = V / sqrt(N)                   # normalize basis vectors
  spots    = dataIn$count[(M-N+1):M]       # N most recent sunspot values
  zeroMean = spots - mean(spots)           # zero-mean sunspot data
  stime    = spotData$date[(M-N+1):M]      # N most recent time values
  coefs    = as.numeric( zeroMean %*% V ); # coef[n] = projection of data onto eigen-vector n
  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
  }
  if( dataPlot )
  {
    plot(  c(0:(plotLength-1)), coefs[1:plotLength], col=colors[1], type='h', lwd=3)
    lines( c(0:(plotLength-1)), coefs[1:plotLength], col=colors[1], type='p', lwd=2)
   #abgrid( xmin=0, xmax=120, xstep=5, ymin=-700, ymax=2000, ystep=100 );
  }
  if( dataDump )
  {
    sink(sprintf("tex/%s.dat",dataFileBase));
    printf("%%=============================================================================\n"  );
    printf("%% %s \n", author                                                                   );
    printf("%% Sunspot Walsh coefficient data (%d coefficients)\n", plotLength                  );
    printf("%% Basis vector length = %d\n", N                                                   );
    printf("%% %s\n", LaTeXstr                                                                  );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                  );
    printf("%% %s\n", AutoGenStr                                                                );
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(n in 1:plotLength) printf("  (%3d, %14.8f)\n", n-1, coefs[n]                            );
    printf("]\n"                                                                                );
    sink();
  }
  L = list( coefs=coefs, V=V, z=zeroMean )
  return(L)
}

#------------------------------------------------------------------------------
# \brief Data synthesis using Walsh vector basis
#------------------------------------------------------------------------------
sunspots_walsh_synth = function(
  verbose      = TRUE,
  dataDump     = FALSE,                     # dump data to file
  dataPlot     = TRUE,                      # plot data
  numCoefs     = 6,                         # number of coefficients to synthesis with
  dataIn       = spotData,                  # sunspot data
  basis,                                    # basis vectors
  dataFileBase = "sunspots_walsh_synth"     # file base name
){
  V        = basis$W                        # basis vectors
  N        = length(V[,1]);                 # number of eigen-pairs
  M        = length(dataIn$count)           # number of sunspot data values
  V        = V / sqrt(N)                    # normalize
  spots    = dataIn$count[(M-N+1):M]        # last N sunspot data values
  estMean  = mean(spots)                    # estimated mean
  zeroMean = spots - estMean                # zero-mean data
  stime    = spotData$date[(M-N+1):M]       # last N sunspot time values
  coefs    = as.numeric( zeroMean %*% V );  # projection coefficient of spots onto eigen-vector n
  G = sqrt(as.numeric(zeroMean%*%zeroMean)) # estimated energy of sunspot waveform
  fsyn     = 0 * c(1:N)
  g        = 0;
  for( n in 1:numCoefs )
  {
    fsyn = fsyn + coefs[n] * V[,n]          # partially-synthesized vector
    g = g + (coefs[n])^2                    # energy of scaled eigen-vectors
  }
  if(g>1e-12){fsyn = ((G/sqrt(g)) * fsyn)}  # scale vector to match energy of original data
  fsyn      = fsyn + estMean                # restore mean
  errorVect = fsyn - spots                  # calculate error vector
  rmsError  = sqrt( (errorVect %*% errorVect))/N # RMS error
  if( verbose )
  {
    printf("%s(...)\n", dataFileBase );
    printf("  Total RMS synthesis error using %d coefficients = %12.8f\n", numCoefs, rmsError );
    printf("  Energy of zero mean vector         G = %12.8f\n", G)
    printf("  Energy of synthesized vector sqrt(g) = %12.8f\n", sqrt(g))
    printf("  Energy ratio G/sqrt(g)               = %12.8f\n", G/sqrt(g))
  }
  if( dataPlot )
  {
    plot( stime, spots, col=colors[1], type='l', ylim=c(-25,375))
    lines(stime, fsyn,  col=colors[2], type='l', lwd=3)
    abgrid( xmin=1850, xmax=2020, xstep=10, ymin=0, ymax=350, ystep=50 );
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
    printf("%% %s\n", LaTeXstr                                                                  );
    printf("%% For an example, see \"%s.tex\"\n", dataFileBase                                  );
    printf("%% %s\n", AutoGenStr                                                                );
    printf("%% Total RMS synthesis error using %d coefficients = %12.8f\n", numCoefs, rmsError  );
    printf("%% Basis vector length = %d\n", N                                                   );
    printf("%%=============================================================================\n"  );
    printf("[\n"                                                                                );
    for(n in 1:length(fsyn)) printf("  (%12.8f, %12.8f)\n", stime[n], fsyn[n]                   );
    printf("]\n"                                                                                );
    sink();
  }
  L = list( z=zeroMean, V=V, N=N, M=M )
  return( L )
}

#------------------------------------------------------------------------------
# Main Processing
#------------------------------------------------------------------------------
 T = TRUE
 F = FALSE

 spotData  = sunspots_tseries_data( verbose=F, dataDump=F, dataPlot=T                                    );
 acfData   = sunspots_tseries_acf(  verbose=F, dataDump=F, dataPlot=T,            );
 psdCoefs  = sunspots_psd_coefs(    verbose=F, dataDump=F, dataPlot=T, numSegments=4 );
#eigenBasis= sunspots_eigen_basis(  verbose=F, dataDump=F, dataPlot=T, windowLength=2001     );
#eigenCoefs= sunspots_eigen_coefs(  verbose=F, dataDump=F, dataPlot=T, basis=eigenBasis, dataIn=spotData, Length=105 );
#eigenSynth= sunspots_eigen_synth(  verbose=F, dataDump=F, dataPlot=T, basis=eigenBasis, numCoefs=6   );
#eigenACF  = sunspots_eigen_acf(    verbose=F, dataDump=F, dataPlot=T, coefs=eigenCoefs,Length=100 );
 dftBasis  = sunspots_dft_basis(    verbose=F, dataDump=F, dataPlot=T, windowLength=2001,  numVectors=5 );
 dftCoefs  = sunspots_dft_coefs(    verbose=F, dataDump=F, dataPlot=T, basis=dftBasis, dataIn=spotData, plotLength=1001 );
 dftSynth  = sunspots_dft_synth(    verbose=F, dataDump=F, dataPlot=T, basis=dftBasis, numCoefs=17  );
 dftACF    = sunspots_dft_acf(      verbose=F, dataDump=F, dataPlot=T, coefs=dftCoefs, Length=100 );
 walshBasis= sunspots_walsh_basis(  verbose=F, dataDump=F, dataPlot=T, windowLength=2048,  numVectors=5 );
 walshCoefs= sunspots_walsh_coefs(  verbose=F, dataDump=F, dataPlot=T, basis=walshBasis, dataIn=spotData, plotLength=2048 );
 walshSynth= sunspots_walsh_synth(  verbose=T, dataDump=T, dataPlot=T, basis=walshBasis, dataIn=spotData, numCoefs=35  );
#f = dftBasis$f
#V = dftBasis$V
#W = dftBasis$W
#N = dftBasis$N
#M = dftBasis$M
#z = complex( real=V, imaginary=W )
 W = walshBasis$W
 f = walshBasis$f
N = walshSynth$N
M = walshSynth$M
z = walshSynth$z
V = walshSynth$V

