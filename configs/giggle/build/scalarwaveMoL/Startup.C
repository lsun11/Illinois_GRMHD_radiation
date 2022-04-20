/*
  Register banner that will be printed out at startup
*/
#include "cctk.h"

extern "C" int scalarwaveMoL_Startup(void)
{
   CCTK_RegisterBanner("ABE's Scalar Field Evolver, using MoL for timestepping");
   return 0;
}
