/*
   Set the symmetries for the Scalarwave solver
*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "Symmetry.h"

static char *rcsid="$Header: /home/astro/CVS/Cactus/arrangements/ABE/scalarwave/src/InitSymBound.C,v 1.2 2006/03/10 22:32:24 zetienne Exp $";

CCTK_FILEVERSION(scalarwaveMoL_InitSymBound)

extern "C" void scalarwaveMoL_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
    
  int sym[3];

  //Both phi and phidot have +1 symmetries across any symmetry axis:
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;

  SetCartSymVN(cctkGH, sym,"scalarwaveMoL::phi");
  SetCartSymVN(cctkGH, sym,"scalarwaveMoL::phidot");

  return;
}
