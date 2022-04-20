/*
 *  Subroutine that checks to make sure initial parameters are not crazy.
 */

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

extern "C" void scalarwaveMoL_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    //This function does nothing at the moment.  
    //Given the simplicity of this thorn, not sure it's worth the effort.

  return;
}

