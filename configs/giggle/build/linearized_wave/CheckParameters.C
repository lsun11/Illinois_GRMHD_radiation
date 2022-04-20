/*
 *  Subroutines for checking that SOME initial parameters are not crazy.
 *  This is not incredibly useful. -Zach
 */

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

static char *rcsid = "$Header: Exp $";

CCTK_FILEVERSION(ABE_linearizedwave_CheckParameters)

extern "C" void LinearizedWave_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(CCTK_Equals(initial_data,"box"))
  {
    if (!CCTK_Equals(type, "box"))
    {
      CCTK_PARAMWARN("Must have a box grid with box initial data");
    }

  }
  return;
}
