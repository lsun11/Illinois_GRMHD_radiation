/*
 *  Subroutines for checking that initial parameters are not crazy.
 */

#include "stdio.h" 
#include "stdlib.h" 
#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
//#include "cctk_WarningLevel.h"

static char *rcsid = "$Header: ...pretend like there's some CVS crap here...$";

CCTK_FILEVERSION(BSSN_CheckParameters)

extern "C" void BSSN_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

    if(number_of_mol_ministeps <= 0) {
      printf("YOU FORGOT TO SET THE number_of_mol_ministeps PARAMETER IN YOUR par FILE!\n");
      exit(1);
    }

  return;
}
