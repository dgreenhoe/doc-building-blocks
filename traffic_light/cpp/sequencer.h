//============================================================================
//! \author  Daniel J. Greenhoe
//============================================================================
#include <cmath>

enum faceState { green, yellow, red, arrow, yellow2, red2 };
enum lightFace { master, back, left, right };

//-----------------------------------------------------------------------------
//! \brief A single traffic light face consisting of four lights of 
//!        enum type <faceState>
//-----------------------------------------------------------------------------
class trafficFace {
  private:
    faceState state;          //! \param[inout] state: current state of light
  public:                     
    trafficFace  (void);      // constructor
    trafficFace  (faceState); // constructor
   ~trafficFace  (void);      // destructor
    faceState get(void);      // get current state
    char* getStr (faceState fs, char *buf); // get string for a given <fs>
    void step    (void);      // step current state of traffic light face
    void operator++(int);     // same as step()
};

//-----------------------------------------------------------------------------
//! \brief A traffic light consisting of faces, each with four lights of 
//!        enum type <faceState>
//-----------------------------------------------------------------------------
class trafficLight : public trafficFace {
  public:
    trafficLight     (void)             :  trafficFace()     {}; // constructor
    trafficLight     (faceState state)  :  trafficFace(state){}; // constructor
    faceState getM   (void);      // get state of MASTER face
    faceState getB   (void);      // get state of BACK   face
    faceState getL   (void);      // get state of LEFT   face
    faceState getR   (void);      // get state of RIGHT  face
    char*     getStrM(char *buf); // get state of MASTER as string
    char*     getStrB(char *buf); // get state of BACK   as string
    char*     getStrL(char *buf); // get state of LEFT   as string
    char*     getStrR(char *buf); // get state of RIGHT  as string
};

