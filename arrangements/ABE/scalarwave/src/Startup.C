/*
  Register banner that will be printed out at startup
*/
#include "cctk.h"

static char *rcsid = "$Header: /home/astro/CVS/Cactus/arrangements/ABE/scalarwave/src/Startup.C,v 1.1.1.1 2006/02/17 21:31:26 zetienne Exp $";

CCTK_FILEVERSION(ABE_scalarwave_Startup)

extern "C" int ABE_scalarwave_Startup(void)
{
   CCTK_RegisterBanner("ABE's Scalar Field Solver");
   return 0;
}
