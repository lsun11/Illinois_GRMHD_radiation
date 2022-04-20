/*
  Register banner that will be printed out when the code starts.
  Another lame function. -Zach
*/
#include "cctk.h"

static char *rcsid = "$Header: Exp $";

CCTK_FILEVERSION(ABE_linearizedwave_Startup)

extern "C" int ABE_linearizedwave_Startup(void)
{
   CCTK_RegisterBanner("ABE's Linearized Wave Evolver");
   return 0;
}
