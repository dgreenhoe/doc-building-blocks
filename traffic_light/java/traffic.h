//============================================================================
//! \author  Daniel J. Greenhoe
//=============================================================================

void estimate_traffic_arrival_times(
  double mean,                    //! \param[in]  mean: mean of random process
  std::vector<double> & arrivals  //! \param[out] arrivals: PRNG random values
  );
void induction_loop_samples(
  double amplitude,               //! \param[in]  amplitude: amplitude of sinusoid
  double frequency,               //! \param[in]  frequency: frequency of sinusoid
  double mean,                    //! \param[in]  mean: mean of random process
  double stdDev,                  //! \param[in]  stdDev: standard deviation
  std::vector<double> & samples   //! \param[out] samples: induction loop samples
  );
