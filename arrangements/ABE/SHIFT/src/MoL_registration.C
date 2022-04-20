//--------------------------------------------------------------------------
// Register with the time stepper 
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <stddef.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

extern "C" void shift_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0, group, rhs, var;

  // Register evolution & RHS gridfunction groups
  group = CCTK_GroupIndex("shift::shift_vars");
  rhs = CCTK_GroupIndex("shift::shift_rhs");

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }

  if (ierr) CCTK_WARN(0,"Problems registering with MoL");

}
