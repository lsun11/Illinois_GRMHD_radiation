/*
  Register banner that will be printed out at startup
*/
#include "cctk.h"

static char *rcsid = "$Header: yeah, right. $";

CCTK_FILEVERSION(Shift_Startup)

extern "C" int Shift_Startup(void)
{
   CCTK_RegisterBanner("Shift thorn is starting up.");
   return 0;
}
