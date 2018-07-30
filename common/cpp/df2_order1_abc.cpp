//===========================================================================
//! Daniel J. Greenhoe
//! \brief DFII order 1 filter
//! \code{.markdown}
//!                           p     |\a
//!   x[n] -->--(+)---->------o-----| |---(+)-->-- y[n]
//!              |            |     |/     |
//!              |           \|/           |
//!              |          __|__          |
//!              |         | -1  |         |
//!             /|\        |z    |        /|\
//!              |         |_____|         |
//!              |            |            |
//!              |   -c/|     |     |\b    |
//!              o----| |--<--o-->--| |----o
//!                    \|   state   |/
//! \endcode
//===========================================================================
double df2_order1_filter(//! \return     Return state value
  const double a,        //! \param[in]  a      filter coefficient a
  const double b,        //! \param[in]  b      filter coefficient b
  const double c,        //! \param[in]  c      filter coefficient c
  const double state,    //! \param[in]  state  state of state-machine filter
  const float *x,        //! \param[in]  x      pointer to input  data
        float *y,        //! \param[out] y      pointer to output data
  const long   N         //! \param[in]  N      length of x, y;
  )
{
  long n;
  double xn,yn;
  for(n=0; n<N; n++)
  {
    xn    = (double)x[n];    // convert float to double
    p     = xn  - c*state;
    yn    = a*p + b*state;
    y[n]  = (float)yn;       // convert double to float
    state = p;               // update state
  }
  return state;
}
