//============================================================================
//! \author  Daniel J. Greenhoe
//============================================================================
#include <cmath>
double estimate_induction(//! \return estimated induction L
  const double A,    //! \param[in] A:    input amplitude as in A_in sin(2pi ft)
  const double Aest, //! \param[in] Aest: estimated/measured amplitude
  const double R,    //! \param[in] R:    resistance in test circuit
  const double f     //! \param[in] f:    frequency of test signal
  );
