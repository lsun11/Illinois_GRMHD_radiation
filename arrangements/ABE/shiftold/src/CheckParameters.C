/*
 *  Subroutines for checking that initial parameters are not crazy.
 */

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
//#include "cctk_WarningLevel.h"

static char *rcsid = "$Header: ...pretend like there's some CVS crap here...$";

CCTK_FILEVERSION(Shift_CheckParameters)

extern "C" void Shift_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(Spatial_Gauge==0)
  {
    CCTK_INFO("Evolving with zero shift...");
  }

  return;
}

