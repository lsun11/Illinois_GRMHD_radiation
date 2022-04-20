#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Compute b^2 and store it in Pr 
subroutine mag_bondi_compute_b2(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (puncture_id==1) then 
     call mag_bondi_compute_b2_cpp(cctkGH,cctk_lsh, phi, lapm1, &
                           shiftx,shifty,shiftz,vx,vy,vz,Bx,By,Bz, & 
                           gxx, gxy, gxz, gyy, gyz, gzz, Pr)
  end if

end subroutine mag_bondi_compute_b2
