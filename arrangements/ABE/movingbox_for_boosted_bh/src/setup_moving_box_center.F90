#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Setup movingbox center.  
!----------------------------------------------------------------------------------
subroutine setup_movingbox_center(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
!
  xc_prev = x0_bh
  yc_prev = 0.d0
  zc_prev = 0.d0
  position_x(1) = x0_bh
  position_y(1) = 0.d0
  position_z(1) = 0.d0
  
end subroutine setup_movingbox_center
