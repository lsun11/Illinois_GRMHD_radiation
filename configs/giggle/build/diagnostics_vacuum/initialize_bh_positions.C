//--------------------------------------------------------------------------
// Set Xbh1,Ybh1,Zbh1; and Xbh2,Ybh2,Zbh2 if desired.
//--------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <stddef.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

extern "C" void initialize_bh_posns(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(cctk_iteration==0) {
    for(int i=0;i<3;i++) {
      bh_posn_x[i] = 0.0;
      bh_posn_y[i] = 0.0;
      bh_posn_z[i] = 0.0;

      bh_radius_z[i] = -1.0;
    }
  }
}
