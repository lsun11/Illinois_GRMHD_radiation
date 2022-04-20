/*
  Register banner that will be printed out at startup
*/
#include "cctk.h"

static char *rcsid = "$Header: yeah, right. $";

CCTK_FILEVERSION(Lapse_Startup)

extern "C" int Lapse_Startup(void)
{
   CCTK_RegisterBanner("ABE's Lapse thorn is starting up.");
   return 0;
}
