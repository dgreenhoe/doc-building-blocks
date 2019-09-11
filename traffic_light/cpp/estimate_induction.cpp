//============================================================================
//! \author  Daniel J. Greenhoe
//! \class   Induction
//! \brief   Perform ML-estimation of induction L (Henrys) in an AWGN system
//=============================================================================
#include "estimate_induction.h"
//-----------------------------------------------------------------------------
//! \brief   Estimate induction L from estimated amplitude Aest
//!
//! 
//!     v(t) = i(t) x impedance = 
//!             d                                       V(s)
//!   v(t) = L ---- i(t)   V(s) = sL I(s)   impedance = ---- = sL
//!             dt                                      I(s)
//!                
//!                
//!                resistor  
//!               impedance=R 
//!           o----^v^v^v-----|------------------------o
//!             --> I(s)  --> |
//!                           |
//!          X(s)             O inductor             Y(s)
//!                           O impedance=sL    
//!           o               O                        o  
//!           |               |                        |  
//!           |               |                        |  
//!         -----           -----                    -----
//!          ---             ---                      --- 
//!           -               -                        -  
//!                 X(s)        X(s)    
//! Current I(s) = --------- = --------
//!                impedance    R + sL
//!
//!
//! output voltage Y(s) = I(s) x inductor impedance
//!                     = I(s) sL
//!                        X(s) sL
//!                     = --------
//!                        R + sL
//!
//!              Y(s)       sL
//! gain G(s) = ------ = --------
//!              X(s)     R + sL 
//!
//!     |    |2   |    sL    |2      |   iwL    |2       2|     1    |2   
//!     |G(w)|  = | -------- |     = | -------- |  = (wL) | -------- |    
//!     |    |    |  R + sL  |s=iw   |  R + iwL |         |  R + iwL | 
//!                                                        
//!                   2|     R - iwL   |2       (wL)^2
//!             = (wL) | ------------- |  = --------------
//!                    |  R^2 + (wL)^2 |     R^2 + (wL)^2
//! 
//! Solve for L:
//!          2      g^2 R^2               |    |2
//!         L   = ------------  where g = |G(w)|  
//!                (1-g^2) w^2            |    | 
//-----------------------------------------------------------------------------
double estimate_induction(//! \return estimated induction L
  const double A,         //! \param[in] A:    input amplitude as in A_in sin(2pi ft)
  const double Aest,      //! \param[in] Aest: estimated/measured amplitude
  const double R,         //! \param[in] R:    resistance in test circuit
  const double f          //! \param[in] f:    frequency of test signal
  )
  {
    const double gain        = Aest / A;        // gain
    const double w           = 2 * M_PI * f;    // frequency in radians per second
    const double Lsq         = (gain * gain * R * R ) / (w*w * (1 - gain*gain));
    const double induction_L = sqrt(Lsq);
    return induction_L;
  }

