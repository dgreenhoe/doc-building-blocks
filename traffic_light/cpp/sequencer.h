//============================================================================
//! \author  Daniel J. Greenhoe
//============================================================================
#include <cmath>

enum faceState { red, yellow, green, arrow };

faceState step_masterLight(
  faceState state  //! \param[in] state: current state of master light
  );

