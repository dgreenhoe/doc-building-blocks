//=============================================================================
//! \author  Daniel J. Greenhoe
//! \class   Traffic-Light-Control
//! \brief   Sequence 4 face 5 light/face traffic light
//=============================================================================
#include "sequencer.h"
//-----------------------------------------------------------------------------
//! \brief Step master light state variable
//-----------------------------------------------------------------------------
faceState step_masterLight(
  faceState state  //! \param[in] state: current state of master light
  )
  {
  faceState next;
  switch(state)
  {
    case green  : next = yellow; break;
    case yellow : next = red   ; break;
    case red    : next = arrow ; break;
    case arrow  : next = green ; break;
    default     : next = red   ; break;
  }
  return next;
  }
