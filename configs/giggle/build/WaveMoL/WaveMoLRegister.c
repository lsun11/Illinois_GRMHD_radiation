 /*@@
   @file      WaveMoLRegister.c
   @date      Fri Nov  9 13:47:07 2001
   @author    Ian Hawke
   @desc 
   Routine to register the variables with the MoL thorn.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusExamples/WaveMoL/src/WaveMoLRegister.c,v 1.1.1.1 2004/06/26 14:49:02 hawke Exp $";

CCTK_FILEVERSION(CactusExamples_WaveMoL_WaveMoLRegister_c)

/*
#ifndef DEBUG_MOL
#define DEBUG_MOL
#endif
*/

int WaveMoL_RegisterVars(CCTK_ARGUMENTS);

 /*@@
   @routine    WaveMoL_RegisterVars
   @date       Fri Nov  9 13:47:41 2001
   @author     Ian Hawke
   @desc 
   The registration routine.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int WaveMoL_RegisterVars(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr = 0, group, rhs, var;

  group = CCTK_GroupIndex("wavemol::scalarevolvemol_scalar");
  rhs = CCTK_GroupIndex("wavemol::scalarrhsmol_scalar");

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }

  group = CCTK_GroupIndex("wavemol::scalarevolvemol_vector");
  rhs = CCTK_GroupIndex("wavemol::scalarrhsmol_vector");

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function MoLRegisterEvolvedGroup not aliased");
    ierr++;
  }

  var = CCTK_VarIndex("wavemol::energy");

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

#ifdef DEBUG_MOL
  printf("If we've got this far, then we've done with wavetoy registration.\n");
#endif  

  return ierr;
}
