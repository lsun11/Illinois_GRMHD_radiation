/*
   Set the symmetries for the shift
*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include <stdio.h>

static char *rcsid="$Header: /peter/piper/picked/a/peck/of/pickled/peppers $";

CCTK_FILEVERSION(Shift_InitSymBound)

extern "C" void Shift_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
    
  int sym[3];
  
  if(Symmetry == 2) {
    sym[0] = -1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"Shift::shiftx");
    SetCartSymVN(cctkGH, sym,"Shift::shiftxt");
    
    sym[0]=1; sym[1]=-1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"Shift::shifty");
    SetCartSymVN(cctkGH, sym,"Shift::shiftyt");
    //BUG HERE:
    //sym[1]=1; sym[2]=-1; sym[2] = 1;
    sym[0]=1; sym[1]=1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"Shift::shiftz");
    SetCartSymVN(cctkGH, sym,"Shift::shiftzt");
  }
  else if(Symmetry == 4) {    
    sym[0] = -1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"Shift::shiftx");
    SetCartSymVN(cctkGH, sym,"Shift::shiftxt");
    
    sym[0] = -1; sym[1] = -1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"Shift::shifty");
    SetCartSymVN(cctkGH, sym,"Shift::shiftyt");
    
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"Shift::shiftz");
    SetCartSymVN(cctkGH, sym,"Shift::shiftzt");
  }
  else if(Symmetry == 1) {    
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"Shift::shiftx");
    SetCartSymVN(cctkGH, sym,"Shift::shiftxt");
    
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"Shift::shifty");
    SetCartSymVN(cctkGH, sym,"Shift::shiftyt");
    
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"Shift::shiftz");
    SetCartSymVN(cctkGH, sym,"Shift::shiftzt");
  }
  else if(Symmetry == 0) {
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"Shift::shiftx");
    SetCartSymVN(cctkGH, sym,"Shift::shiftxt");
    SetCartSymVN(cctkGH, sym,"Shift::shifty");
    SetCartSymVN(cctkGH, sym,"Shift::shiftyt");
    SetCartSymVN(cctkGH, sym,"Shift::shiftz");
    SetCartSymVN(cctkGH, sym,"Shift::shiftzt");
  }
  else printf("SYMMETRY = %d NOT SUPPORTED!\n",Symmetry);

  printf("SET SHIFT SYMMETRIES!\n");

  return;
}
