/*
 *  Subroutines for checking that initial parameters are not crazy.
 */

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
//#include "cctk_WarningLevel.h"

static char *rcsid = "$Header: ...pretend like there's some CVS crap here...$";

CCTK_FILEVERSION(Lapse_CheckParameters)

extern "C" void Lapse_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(CCTK_Equals(slicing_type,"geodesic"))
    {
      CCTK_INFO("Evolving with geodesic slicing...");
    }
  if(CCTK_Equals(slicing_type,"harmonic"))
    {
      CCTK_INFO("Evolving with harmonic slicing...");
    }
  if(CCTK_Equals(slicing_type,"hyperbolic"))
    {
      CCTK_INFO("Evolving with hyperbolic slicing...");
    }
  if(CCTK_Equals(slicing_type,"opl"))
    {
      CCTK_INFO("Evolving with one plus log slicing...");
      if(opl_advect_enable == 2) {
	CCTK_INFO("REMEMBER: You must set ghost_zones >= 3 for this choice of upwinding!!");
      }
    }

  return;
}

