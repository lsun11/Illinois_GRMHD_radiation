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

extern "C" void scalarwaveMoL_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0, group, rhs, var;

  // Register evolution & RHS gridfunction groups
  group = CCTK_GroupIndex("scalarwaveMoL::scalarMoLevolve");
  rhs = CCTK_GroupIndex("scalarwaveMoL::scalarMoLrhs");

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }

  // *** TEST ***
  group = CCTK_GroupIndex("scalarwaveMoL::scalarMoLstagger");
  rhs = CCTK_GroupIndex("scalarwaveMoL::scalarMoLstaggerrhs");

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }
  // ************

  printf("AFTER FIRST REGISTRATIONS!\n");

  // Register analytic gridfunction group
  var = CCTK_VarIndex("scalarwaveMoL::phi_analytic");
  var = CCTK_VarIndex("scalarwaveMoL::phi_analytic_minus_numeric");

  printf("VAR = %d\n",var);

  if (CCTK_IsFunctionAliased("MoLRegisterConstrained"))
  {
    ierr += MoLRegisterConstrained(var);
  }
  else
  {
    CCTK_WARN(0, "MoL function MoLRegisterConstrained not aliased");
    ierr++;
  }

  if (ierr) CCTK_WARN(0,"Problems registering with MoL");

  printf("AFTER SECOND REGISTRATION!\n");

}
