/*
   Set the symmetries for the Scalarwave solver
*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "Symmetry.h"

static char *rcsid="$Header: /home/astro/CVS/Cactus/arrangements/ABE/scalarwave/src/InitSymBound.C,v 1.2 2006/03/10 22:32:24 zetienne Exp $";

CCTK_FILEVERSION(ABE_ScalarWave_InitSymBound)

extern "C" void ABE_ScalarWave_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
    
  int sym[3];
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;

  /* needs to be fixed:
  if(CCTK_Equals(domain,"octant")==1) {
    sym[0] = 1;
    sym[1] = 1;
    sym[2] = 1;
  }
  else {
    sym[0] = 1;
    sym[1] = 1;
    sym[2] = 1;
  }
  */

  SetCartSymVN(cctkGH, sym,"ABE_scalar_wave::phi");
  SetCartSymVN(cctkGH, sym,"ABE_scalar_wave::phidot");

  SetCartSymVN(cctkGH, sym,"ABE_scalar_wave::phi_t");
  SetCartSymVN(cctkGH, sym,"ABE_scalar_wave::phidot_t");

  return;
}
