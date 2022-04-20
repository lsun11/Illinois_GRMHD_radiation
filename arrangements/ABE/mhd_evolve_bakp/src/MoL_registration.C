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

extern "C" void mhd_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0, group, rhs, var;

  // Register evolution & RHS gridfunction groups
  if (em_evolve_enable==1 && constrained_transport_scheme==3) {
     group = CCTK_GroupIndex("mhd_evolve::em_Ax");
     rhs = CCTK_GroupIndex("mhd_evolve::em_rhsx");
  } else {
    group = CCTK_GroupIndex("mhd_evolve::em_conservativex");
    rhs = CCTK_GroupIndex("mhd_evolve::em_rhsx");
  }

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }

  if (em_evolve_enable==1 && constrained_transport_scheme==3) {
     group = CCTK_GroupIndex("mhd_evolve::em_Ay");
     rhs = CCTK_GroupIndex("mhd_evolve::em_rhsy");
  } else {
     group = CCTK_GroupIndex("mhd_evolve::em_conservativey");
     rhs = CCTK_GroupIndex("mhd_evolve::em_rhsy");
  }

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }

  if (em_evolve_enable==1 && constrained_transport_scheme==3) {
     group = CCTK_GroupIndex("mhd_evolve::em_Az");
     rhs = CCTK_GroupIndex("mhd_evolve::em_rhsz");
  } else {
     group = CCTK_GroupIndex("mhd_evolve::em_conservativez");
     rhs = CCTK_GroupIndex("mhd_evolve::em_rhsz");
  }

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }

  group = CCTK_GroupIndex("mhd_evolve::em_Phi");
  rhs = CCTK_GroupIndex("mhd_evolve::em_Phi_rhs");

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }

  group = CCTK_GroupIndex("mhd_evolve::em_Blagrangemultiplier");
  rhs = CCTK_GroupIndex("mhd_evolve::em_Blagrangemultiplier_rhs");

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }

  group = CCTK_GroupIndex("mhd_evolve::mhd_conservatives");
  rhs = CCTK_GroupIndex("mhd_evolve::mhd_rhs");


  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }

 
  group = CCTK_GroupIndex("mhd_evolve::rad_conservatives");                                                                
  rhs = CCTK_GroupIndex("mhd_evolve::rad_conservatives_rhs"); 

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
