//============================================================================
//! \author  Daniel J. Greenhoe
//============================================================================
#include <cmath>

enum faceState { green, yellow1, red1, arrow, yellow2, red2 };
enum lightFace { master, back, left, right };

//-----------------------------------------------------------------------------
//! \brief A single traffic light face consisting of four lights of 
//!        enum type <faceState>
//-----------------------------------------------------------------------------
class trafficFace {
  private:
    faceState state;          //! \param[inout] state: current state of light
  public:
    trafficFace (void);      // constructor
    trafficFace (faceState); // constructor
   ~trafficFace (void);      // destructor
    faceState get(void);     // get current state
    char* getStr(char *buf); // Get current state in string representation
    void step   (void);      // step current state of traffic light face
    void operator++(int);    // same as step()
};

//-----------------------------------------------------------------------------
//! \brief A traffic light consisting of faces, each with four lights of 
//!        enum type <faceState>
//-----------------------------------------------------------------------------
class trafficLight : public trafficFace {
  public:
    trafficLight(void)             :  trafficFace()     {}; // constructor
    trafficLight(faceState state)  :  trafficFace(state){}; // constructor
    faceState getM(void);       // get current state of master face
};

