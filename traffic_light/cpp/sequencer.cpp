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
#include "sequencer.h"
//-----------------------------------------------------------------------------
//! \brief trafficFace initializers and destructor
//! \details By default, a traffic light face is initialized to <red>.
//-----------------------------------------------------------------------------
trafficFace:: trafficFace(faceState initialState){ state = initialState; }
trafficFace:: trafficFace(void)                  { trafficFace(red);     }
trafficFace::~trafficFace(void)                  {                       }

//-----------------------------------------------------------------------------
//! \brief Step master light state variable
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
