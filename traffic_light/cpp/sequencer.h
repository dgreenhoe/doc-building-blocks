//============================================================================
//! \author  Daniel J. Greenhoe
//============================================================================
#include <cmath>

enum faceState { red, yellow, green, arrow };

//-----------------------------------------------------------------------------
//! \brief A single traffic light face consisting of four lights of 
//!        enum type <faceState>
//-----------------------------------------------------------------------------
class trafficFace {
  private:
    faceState state;  //! \param[inout] state: current state of light
  public:
    trafficFace(void);      //! initialize state to default
    trafficFace(faceState); //! initialize state to <faceState>
    ~trafficFace(void);
    void step(void);        //! step current state of traffic light face
  };