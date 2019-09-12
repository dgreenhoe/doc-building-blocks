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
    ~trafficFace(void);     //! destructor
    int  getInt(void);      //! get current state in integer representation
    char* getStr(char *buf);//! Get current state in string representation
    void step(void);        //! step current state of traffic light face
    void operator++(void);  //! same as step()
  };