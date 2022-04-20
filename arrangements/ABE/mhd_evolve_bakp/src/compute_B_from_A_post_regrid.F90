!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine compute_B_from_A_post_regrid(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: dummy

! FOLLOWING CARTSYM'S ARE NEEDED.  Otherwise, we'll end up with memory errors inside ghostzones (even BSSN quantities)!
  call CartSymGN(dummy,cctkGH,'bssn::BSSN_vars')

  call compute_B_from_A(CCTK_PASS_FTOF)

end subroutine compute_B_from_A_post_regrid
