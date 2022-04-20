/*
   Set the symmetries for the lapse
*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

static char *rcsid="$Header: /peter/piper/picked/a/peck/of/pickled/peppers $";

CCTK_FILEVERSION(Lapse_InitSymBound)

extern "C" void Lapse_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS 
    
  int sym[3];
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;

  SetCartSymVN(cctkGH, sym,"lapse::lapm1");
  SetCartSymVN(cctkGH, sym,"lapse::lapset");

  if(Symmetry == 2) {
    sym[0] = -1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"lapse::lapsex");
    sym[0]=1; sym[1]=-1;
    sym[0] = 1; sym[1] = -1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"lapse::lapsey");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"lapse::lapsez");
  }
  else if(Symmetry == 4) {
    sym[0] = -1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"lapse::lapsex");
    sym[0] = -1; sym[1] = -1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"lapse::lapsey");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"lapse::lapsez");
  }
  else if(Symmetry == 1) {
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"lapse::lapsex");
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"lapse::lapsey");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"lapse::lapsez");
  }
  else if(Symmetry == 0) {
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"lapse::lapsex");
    SetCartSymVN(cctkGH, sym,"lapse::lapsey");
    SetCartSymVN(cctkGH, sym,"lapse::lapsez");
  }


  return;
}
