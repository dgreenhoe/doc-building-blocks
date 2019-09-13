//============================================================================
//! \author  Daniel J. Greenhoe
//! \class   Induction-Estimation
//! \brief   Perform ML-estimation of induction L (Henrys) in an AWGN system
//! \details An estimate of L can be used in a metal detector and/or in 
//!          an induction loop for vehicle detection at a traffic intersection.
//!          Example: A decrease in L indicates an increased probability
//!          of the event E of a vehicle being present 
//!          (due to increased eddy currents)
//! \code{.markdown} Vocabulary
//!      ML:   maximum likelihood
//!      AWGN: additive white Gaussian noise
//! \endcode
//!
//! \code{.markdown} Hardware block diagram
//!                                             20 kHz sinusoid
//!      _______     _______     _______     ________  |
//!     |       |   |       |   |analog |   |voltage |   resistor  inductor
//!     |Cortex-|   | DAC   |   |low    |   |follower|      R          L
//!     |M4F    |=>=|periph.|->-|pass   |->-|/ ampli-|---^v^v^v-----OOOOOOO--
//!     |       |   |(12bit)|   |filter |   |  fier  |                       |
//!     |_______|   |_______|   |_______|   |________|                   x(t)|
//!        ||                                            Gaussian            |
//!        ||        _______     _______     ________     noise    ________  |
//!       /||\      |       |   |analog |   |voltage |     v(t)   |        | |
//!        ||       | ADC   |   |low    |   |follower|      |     |squaring| |
//!        ||===<===|periph.|---|pass   |-<-|(high   |--<--(+)--<-|circuit |-
//!                 |(12bit)|   |filter |   |imped.) | y(t)       |        |
//!                 |_______|   |_______|   |________|            |________|
//! \endcode
//!
//! \see
//! \cite Leonidas Deligiannidis and Hamid R Arabnia (editors). (2014)
//!       Emerging Trends in Image Processing, Computer Vision and Pattern Recognition
//!       Chapter 4: Traffic Detectors
//!       Publisher: Morgan Kaufmann 
//!       ISBN10:012802092X, ISBN13:9780128020920. 640 pages
//!       https://books.google.com/books?id=2O59BAAAQBAJ&pg=PA233
//!   
//! \cite Swadesh Chaulya and G. M. Prasad. (2016)
//!       Sensing and Monitoring Technologies for Mines and Hazardous Areas: Monitoring and Prediction Technologies
//!       Section 2.2.1 Inductive-Loop Detector
//!       Publisher: Elsevier
//!       ISBN10: 0128031956, ISBN13: 9780128031957. 432 pages
//!       https://books.google.com/books?id=lvReBwAAQBAJ&pg=PA87
//!
//! \cite Mandyam D. Srinath and P.K. Rajasekaran and R. Viswanathan (1996)
//!       ISBN10: 013125295X
//!       Introduction to Statistical Signal Processing with Applications
//!       https://books.google.com/books?id=rIa9QgAACAAJ&pg=PA158
//!
//! \cite Daniel J. Greenhoe (2019)
//!       Theorem 11.7: ML amplitude estimation
//!       A Book Concerning Statistical Signal Processing
//!       https://github.com/dgreenhoe/pdfs/blob/master/abcstat.pdf
//=============================================================================
#include "estimate_amplitude.h"
//-----------------------------------------------------------------------------
//! \brief   Estimate amplitude A of sinusoid x(t)=A sin(2pi ft) in AWGN
//! \details By the "Sufficient Statistic Theorem" (Fisher 1922),
//!          a sufficient statistic for the ML-estimate Aest of A is the projection 
//!          of the received signal y(t) onto a basis for x(t), where
//!          x(t) is the 20 kHz sinusoid squared.
//!          By the Shannon Sampling Theorem, a basis for x(t) is samples
//!          (using the ADC) of x(t) at sample rate Fs >= 2x20kHz.
//!          In this application, let Fs = 10x2x20kHz = 200kHz.
//! \code{.markdown}
//!                   ____n=N                 
//!                   \       [               ]2
//!   Aest = argmin_A  \      [ y(n) - x(n;A) ]
//!                   /___n=1 [               ]
//!
//!                  (     ____n=N                       )
//!                  (  d  \       [               ]2    )
//!        = argmin_A( ---  \      [ y(n) - x(n;A) ]  =0 )
//!                  ( dA  /___n=1 [               ]     )
//!
//!              (    ____n=N      )
//!              ( 1  \            )
//!        = SQRT(---  \      y(n) )
//!              ( N  /___n=1      )
//!
//! \endcode
//-----------------------------------------------------------------------------
double estimate_amplitude( //! \return estimated amplitude
  const std::vector<double> &y   //! \param[in] y: input data
  )
  {
    const double sum  = std::inner_product( y.begin(), y.end(), y.begin(), 0.0 );
    const long   N    = y.size();
    const double Aest = sqrt(sum/(double)N);
    return Aest;
  }
