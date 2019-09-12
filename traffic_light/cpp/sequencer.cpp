//=============================================================================
//! \author  Daniel J. Greenhoe
//! \class   Traffic-Light-Control
//! \brief   Sequence 4 face 5 light/face traffic light
//! \details The traffic light is assumed to have 4 faces:
//!          master, back, left, right.
//! \code{.markdown}
//!             _________________
//!            /.               /|
//!           / .              / |
//!          /  .             /  |
//!         /   .            /   |
//!        /________________/   D|
//!        |    .  ______   |  E |
//!        |    . (      )  | R  |
//!        | R  . (  RED )  |    |
//!        |    . (______)  |    |
//!        |    .           |   .|
//!        |    .  ______   |  L |
//!        |    . (      )  | E  |
//!        | Y  . (YELLOW)  |Y   |
//!        |    . (______)  |    |
//!        |    .           |   .|
//!        |    .  ______   |  E |
//!        |    . (      )  | R  |
//!        | G  . (GREEN )  |G   |
//!        |    . (______)  |    |
//!        |    .           |    |
//!        |    .  ______   |back|
//!        |    .-(ARROW )--|----
//!       l|eft.  ( <=== )  |   /
//!        |  .   (______)  |  /
//!        | .              | /right
//!        |._______________|/
//!              master
//!
//! \endcode
//=============================================================================
#include <string.h>
#include "sequencer.h"
//-----------------------------------------------------------------------------
//! \brief trafficFace initializers and destructor
//! \details By default, a traffic light face is initialized to <red>.
//-----------------------------------------------------------------------------
trafficFace:: trafficFace(faceState initialState){ state = initialState; }
trafficFace:: trafficFace(void)                  { trafficFace(red);     }
trafficFace::~trafficFace(void)                  {                       }

//-----------------------------------------------------------------------------
//! \brief Get current state in integer representation
//-----------------------------------------------------------------------------
faceState trafficFace::get(void){ return state; }

//-----------------------------------------------------------------------------
//! \brief Get current state in null-terminated string representation
//-----------------------------------------------------------------------------
char* trafficFace::getStr(faceState fs, char *buf)
{
  switch(fs)
  {
    case green   : strcpy( buf, "green"  ); break;
    case yellow  : strcpy( buf, "yellow" ); break;
    case yellow2 : strcpy( buf, "yellow" ); break;
    case red     : strcpy( buf, "red"    ); break;
    case red2    : strcpy( buf, "red"    ); break;
    case arrow   : strcpy( buf, "arrow"  ); break;
    default      : strcpy( buf, "ERROR"  ); break;
  }
  return buf;
}

//-----------------------------------------------------------------------------
//! \brief   Step traffic light master face state.
//! \details Note that an alternative is the overloaded ++ operator.
//-----------------------------------------------------------------------------
void trafficFace::step(void)
{
  faceState nextState;
  switch(state)
  {
    case green   : nextState = yellow ; break;
    case yellow  : nextState = red    ; break;
    case red     : nextState = arrow  ; break;
    case arrow   : nextState = yellow2; break;
    case yellow2 : nextState = red2   ; break;
    case red2    : nextState = green  ; break;
    default      : nextState = red    ; break; // Naive error handling 
  }                                           // (should throw exception)
  state = nextState;
}

//-----------------------------------------------------------------------------
//! \brief   Step operator ++ for traffic light master face.
//! \details Note that an alternative is the function step();
//! \cite    https://docs.microsoft.com/en-us/cpp/cpp/increment-and-decrement-operator-overloading-cpp
//-----------------------------------------------------------------------------
void trafficFace::operator++(int){ step(); }

//-----------------------------------------------------------------------------
//! \brief Get current state of MASTER face
//-----------------------------------------------------------------------------
faceState trafficLight::getM   (void)     { return get();                }

//-----------------------------------------------------------------------------
//! \brief Get current state of BACK face
//-----------------------------------------------------------------------------
faceState trafficLight::getB(void){ return getM(); } // same as master face

//-----------------------------------------------------------------------------
//! \brief Get current state of LEFT face
//-----------------------------------------------------------------------------
faceState trafficLight::getL(void)
{
  const faceState masterState = get(); // get master state
  faceState leftState;
  switch(masterState)
  {
    case green   : leftState = yellow ; break;
    case yellow  : leftState = red    ; break;
    case red     : leftState = arrow  ; break;
    case arrow   : leftState = yellow2; break;
    case yellow2 : leftState = red2   ; break;
    case red2    : leftState = green  ; break;
    default      : leftState = red    ; break; // Naive error handling 
  }
  return leftState;
}

//-----------------------------------------------------------------------------
//! \brief Get current state of RIGHT face
//-----------------------------------------------------------------------------
faceState trafficLight::getR(void){ return getL(); } // same as left face

//-----------------------------------------------------------------------------
//! \brief Get current state of faces as strings
//-----------------------------------------------------------------------------
char* trafficLight::getStrM(char *buf){ return getStr( getM(), buf ); }
char* trafficLight::getStrB(char *buf){ return getStr( getB(), buf ); }
char* trafficLight::getStrL(char *buf){ return getStr( getL(), buf ); }
char* trafficLight::getStrR(char *buf){ return getStr( getR(), buf ); }
