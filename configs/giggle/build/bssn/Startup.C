/*
  Register banner that will be printed out at startup
*/
#include "cctk.h"

static char *rcsid = "$Header: yeah, right. $";

CCTK_FILEVERSION(BSSN_Startup)

extern "C" int BSSN_Startup(void)
{
   CCTK_RegisterBanner("ABE's BSSN thorn is starting up.");
   return 0;
}
