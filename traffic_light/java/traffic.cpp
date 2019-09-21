//============================================================================
//! \author  Daniel J. Greenhoe
//! \class   Traffic
//=============================================================================
#include <random>
#include "traffic.h"

//-----------------------------------------------------------------------------
//! \brief   Traffic arrival statistical model
//! \details Based on the Poisson distribution
//! \see
//! \cite<Guberinic2007>
//!      Slobodan Guberinic, Gordana Senborn, and Bratislav Lazic (2007)
//!      Optimal Traffic Control: Urban Intersections
//!      See page 19, Section "2.2.1. Modeling arrival flow with the Poisson Process"
//!      CRC Press. ISBN10: 1420062700 ISBN13: 9781420062700. 364 pages
//!      https://books.google.com/books?id=l89KVDzAzPYC&pg=PA19
//-----------------------------------------------------------------------------
void estimate_traffic_arrival_times(
  double mean,                    //! \param[in]  mean: mean of random process
  std::vector<double> & arrivals  //! \param[out] arrivals: PRNG random values
  )
{
  long i;
  const long N = arrivals.size();
  std::mt19937 prng; // Use Mersenne Twister type PRNG
  std::poisson_distribution<int> X(1000. * mean);
  for(i=0; i<N; i++) arrivals[i] = (double)X(prng) / 1000.;
}

//-----------------------------------------------------------------------------
//! \brief   Inductive Loop samples
//! \details Based on the Gaussian (Normal) distribution
//-----------------------------------------------------------------------------
void induction_loop_samples(
  double amplitude,               //! \param[in]  amplitude: amplitude of sinusoid
  double frequency,               //! \param[in]  frequency: frequency of sinusoid
  double mean,                    //! \param[in]  mean: mean of random process
  double stdDev,                  //! \param[in]  stdDev: standard deviation
  std::vector<double> & samples   //! \param[out] samples: induction loop samples
  )
{
  long i;
  const long N = samples.size();
  double sinusoid;
  std::random_device random_device{};
  std::mt19937 prng{random_device()}; // Use Mersenne Twister type PRNG
  std::normal_distribution<double> X{mean, stdDev};
  for(i=0; i<N; i++)
  { 
    sinusoid = amplitude * cos(2*M_PI*frequency*(double)i/1000.);
    samples[i] = X(random_device) + sinusoid*sinusoid*0.01;
  }
}

