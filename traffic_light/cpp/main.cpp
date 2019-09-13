//============================================================================
//! \author  Daniel J. Greenhoe
//! \class   Main
//! \brief   main execution handler
//=============================================================================
#include <cstdio>
#include "estimate_amplitude.h"
#include "estimate_induction.h"
#include "sequencer.h"
#include "traffic.h"

int main(void)
{
  const long N = 1000;
  const long Nsamples = N*100;
  const double mean=20;
  int i, j;
  char buf[256];
  trafficLight tl(red);
  const double amplitude  = 10.;
  const double Resistance = 10000;
  const double frequency  = 20e3;
  const double stdDev     = 10.;
  std::vector<double> arrival_times(N);
  std::vector<double> samples(Nsamples);
  std::vector<double> sample_set(100);
  double Aest = 0;
  double induction_L = 0;
  estimate_traffic_arrival_times( mean, arrival_times );
  induction_loop_samples(amplitude, frequency, 0.0, stdDev, samples);

  printf("Aest      = %lf\n",Aest);
  printf("induction = %lf Henrys\n", induction_L);
  printf("%-12s %-12s %-12s %-12s\n","master", "back", "left", "right");
  printf("----------------------------------------------------\n");
  for(i=0; i<N; i++){
    for(j=0; j<100; j++) sample_set[j] = samples[100*i+j];
    Aest = estimate_amplitude( sample_set );
    induction_L = estimate_induction(amplitude, Aest, Resistance, frequency);
    if( induction_L < 0.5 && arrival_times[i] < 20.0 ) 
    {
      tl++;  // sequence light
      printf("%-8s(%d)  ", tl.getStrM(buf), tl.getM());
      printf("%-8s(%d)  ", tl.getStrB(buf), tl.getB());
      printf("%-8s(%d)  ", tl.getStrL(buf), tl.getL());
      printf("%-8s(%d)  ", tl.getStrR(buf), tl.getR());
    }
    else printf("%52s","");
    printf("L=%6.3lf  Aest=%6.3lf arrival=%6.3lf", induction_L, Aest, arrival_times[i]);
    printf("\n");
  }
  return 0;
}

