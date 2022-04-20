//----------------------------------------------------------------------------
//
// $Id: $
//
//----------------------------------------------------------------------------
//
// 
//
//----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <string.h>
#include <fstream.h>
#include "Symmetry.h"

#define N_INT 1000

static char *rcsid = "$Meow...$";
CCTK_FILEVERSION(mag_bondi_initialdata_part2p5)

extern "C" void mag_bondi_initialdata_part2p5(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
 // Find xmax and zmax 
  int ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                            CCTK_VARIABLE_REAL,zmax_bondi,1,CCTK_VarIndex("grid::Z"));

  ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                            CCTK_VARIABLE_REAL,xmax_bondi,1,CCTK_VarIndex("grid::X"));
  printf("HI XMAX_BONDI, ZMAX_BONDI: %e\t%e\n",xmax_bondi,zmax_bondi);

}
