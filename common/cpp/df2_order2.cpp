//===========================================================================
//! Daniel J. Greenhoe
//! \brief DFII order 1 filter
//! \code{.markdown}
//!                           p     |\b0
//!   x[n] -->--(+)---->------o-----| |---(+)-->-- y[n]
//!              |            |     |/     |
//!              |           \|/           |
//!              |          __|__          |
//!              |         | -1  |         |
//!             /|\        |z    |        /|\
//!              |         |_____|         |
//!              |            |s1          |
//!              |  -a1/|     |     |\b1   |
//!             (+)---| |--<--o-->--| |---(+)
//!              |     \|     |     |/     |
//!              |           \|/           |
//!              |          __|__          |
//!              |         | -1  |         |
//!             /|\        |z    |        /|\
//!              |         |_____|         |
//!              |            |s2          |
//!              |  -a2/|     |     |\b2   |
//!              o----| |--<--o-->--| |----o
//!                    \|           |/
//! \endcode
//===========================================================================
void df2_order1_filter(//! \return     Return state value
  const double a1,       //! \param[in]  a1     filter coefficient a1
  const double a2,       //! \param[in]  a2     filter coefficient a2
  const double b0,       //! \param[in]  b0     filter coefficient b0
  const double b1,       //! \param[in]  b1     filter coefficient b1
  const double b2,       //! \param[in]  b2     filter coefficient b2
        double *s1,      //! \param[in]  state  state of state-machine filter
        double *s2,      //! \param[in]  state  state of state-machine filter
  const float *x,        //! \param[in]  x      pointer to input  data
        float *y,        //! \param[out] y      pointer to output data
  const long   N         //! \param[in]  N      length of x, y;
  )
{
  long n;
  double xn,yn;
  for(n=0; n<N; n++)
  {
    xn    = (double)x[n];      // convert float to double
    p     = xn  -a1*s1 -a2*s2; 
    yn    = b0*p + b1*s1 + b2*s2;
    y[n]  = (float)y n;       // convert double to float
    *s2   = *s1;              // update state
    *s1   = p;                   // update state
  }
}
