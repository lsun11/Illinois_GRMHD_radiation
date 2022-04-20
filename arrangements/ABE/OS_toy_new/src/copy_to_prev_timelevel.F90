!-----------------------------------------------------------------------------
!
! $Id: lin_wave.f90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!
!-----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine OS_rad_copy_to_prev_timelevel(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  write(*,*) "Start OS_rad_copy_to_prev_timelevel"
  !particle tracer stuff
  u0_p=u0
  vx_p=vx
  vy_p=vy
  vz_p=vz

end subroutine OS_rad_copy_to_prev_timelevel
