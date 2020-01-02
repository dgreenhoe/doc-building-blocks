//============================================================================
//! \author  Daniel J. Greenhoe
//! \class   Main
//! \brief   main execution handler
//=============================================================================
#include <cstdio>
#include <cstring>
#include "main.h"
#include "help.h"
#include "command.h"

//-----------------------------------------------------------------------------
//! \brief   main execution handler
//! \details entry point for executable code
//! \see
//!   https://www.tutorialspoint.com/cprogramming/c_command_line_arguments.htm
//-----------------------------------------------------------------------------
int main(const int argc, const char *argv[])
{
  command(argc, argv);
  if(argc<2) help(NULL);
  if(argc>10) defaultx(argc, argv);
  return 0;
}

