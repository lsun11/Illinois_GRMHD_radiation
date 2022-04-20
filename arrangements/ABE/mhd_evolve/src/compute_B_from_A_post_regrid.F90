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

!  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
!  call CartSymGN(dummy,cctkGH,'shift::shift_vars')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativex')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativey')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativez')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Blagrangemultiplier')
  
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_conservatives')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::microphys_primitives')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::micphys_conservatives')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_pressure')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::temperatures')
    

  call compute_B_from_A(CCTK_PASS_FTOF)

end subroutine compute_B_from_A_post_regrid
