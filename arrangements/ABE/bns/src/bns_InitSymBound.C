/*
  Set the symmetries for the mhd_evolve variables
*/

#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

static char *rcsid="$Header: /peter/piper/picked/a/peck/of/pickled/whatever $";

CCTK_FILEVERSION(bns_InitSymBound)

  extern "C" void bns_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int sym[3];
  printf("SYMMETRY = %d\n",Symmetry);
  /* NEXT EQUATORIAL SYMMETRY: sym[0] AND sym[1] DO NOT MATTER! */
  if(Symmetry==1) {
    sym[0] = 1; 
    sym[1] = 1; 

    sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"bns::ut"); 
    SetCartSymVN(cctkGH, sym,"bns::q"); 

    sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"bns::ux"); 

    sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"bns::uy"); 

    sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"bns::uz"); 

  }
  else {
    printf("BNS_initsymbound: SYMMETRY = %d NOT SUPPORTED (yet or ever)!  Exiting...\n",Symmetry);
    exit(1);
  }
  return;
}

