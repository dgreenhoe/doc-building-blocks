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
int trafficFace::getInt(void){ return state; }

//-----------------------------------------------------------------------------
//! \brief Get current state in null-terminated string representation
//-----------------------------------------------------------------------------
char* trafficFace::getStr(char *buf)
{
  switch(state)
  {
    case green  : strcpy( buf, "green"  ); break;
    case yellow : strcpy( buf, "yellow" ); break;
    case red    : strcpy( buf, "red"    ); break;
    case arrow  : strcpy( buf, "arrow"  ); break;
    default     : strcpy( buf, "ERROR"  ); break;
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
    case green  : nextState = yellow; break;
    case yellow : nextState = red   ; break;
    case red    : nextState = arrow ; break;
    case arrow  : nextState = green ; break;
    default     : nextState = red   ; break;
  }
  state = nextState;
}

//-----------------------------------------------------------------------------
//! \brief   Step operator ++ for traffic light master face.
//! \details Note that an alternative is the function step();
//-----------------------------------------------------------------------------
void trafficFace::operator++(void){ step(); }

