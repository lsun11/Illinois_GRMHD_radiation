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

CCTK_FILEVERSION(emfields_ct_InitSymBound)

  extern "C" void emfields_ct_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (em_evolve_enable==1 && Reconstruction==2 && reconstruct_Bitildes_instead_of_Bis==1) {
     printf("PPM does not support Bitildes reconstruction (yet).\n");
     exit(1);
  }

  int sym[3];
  printf("SYMMETRY = %d\n",Symmetry);

  // ZACH SAYS: IGNORE THIS FILE.  IT JUST SETS THE RHS VARIABLES, WHICH DO _NOT_ NEED SYMMETRIES TO BE SET!

  /* FIRST OCTANT SYMMETRY */
  if(Symmetry==2) {

    sym[0] = -1; sym[1] = 1; sym[2] = -Sym_Bz;
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bxtilde_or_Ax_rhs"); 
    sym[0] = 1; sym[1] = -1; sym[2] = -Sym_Bz;
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bytilde_or_Ay_rhs"); 
    sym[0] = 1; sym[1] = 1; sym[2] = Sym_Bz;
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bztilde_or_Az_rhs"); 
  }
  /* NEXT, AXISYMMETRY */
  else if (Symmetry==4) {   
    sym[0] = -1; sym[1] = 1; sym[2] = -Sym_Bz;
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bxtilde_or_Ax_rhs"); 
    sym[0] = -1; sym[1] = -1; sym[2] = -Sym_Bz;
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bytilde_or_Ay_rhs"); 
    sym[0] = 1; sym[1] = 1; sym[2] = Sym_Bz;
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bztilde_or_Az_rhs"); 
  }
  /* NEXT, EQUATORIAL SYMMETRY */
  else if (Symmetry==1) {   
    sym[0] = 1; sym[1] = 1; sym[2] = -Sym_Bz;
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bxtilde_or_Ax_rhs"); 
    sym[0] = 1; sym[1] = 1; sym[2] = -Sym_Bz;
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bytilde_or_Ay_rhs"); 
    sym[0] = 1; sym[1] = 1; sym[2] = Sym_Bz;
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bztilde_or_Az_rhs"); 
  }
  else if (Symmetry==0) {
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bxtilde_or_Ax_rhs");
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bytilde_or_Ay_rhs");
    SetCartSymVN(cctkGH, sym,"mhd_evolve::Bztilde_or_Az_rhs");
  }
  else {
    printf("SYMMETRY = %d NOT SUPPORTED (yet)!  Exiting...\n",Symmetry);
    exit(1);
  }
  return;
}

