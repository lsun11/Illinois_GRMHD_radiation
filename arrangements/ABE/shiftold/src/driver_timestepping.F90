!------------------------------------------------
! Master timestepping driver for the shift thorn
!------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine shift_timestepping(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !cowling_enable=1 disables bssn+lapse+shift
  if (Spatial_Gauge.eq.0 .or. Spatial_Gauge.eq.3 .or. cowling_enable.eq.1) then
     call zero_shift_timestepping(CCTK_PASS_FTOF)
  else if (Spatial_Gauge.eq.1) then
     call hyperbolic_shift_timestepping(CCTK_PASS_FTOF)
  else if (Spatial_Gauge.eq.5) then
     call hbpuncture_shift_timestepping(CCTK_PASS_FTOF)
  else if (Spatial_Gauge.eq.6) then
     call puncture_firstorder_shift_timestepping(CCTK_PASS_FTOF)
  end if
  
end subroutine shift_timestepping

subroutine zero_shift_timestepping(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  ! Copy shift
  !
  shiftx_rhs = 0.D0
  shifty_rhs = 0.D0
  shiftz_rhs = 0.D0
  
end subroutine zero_shift_timestepping
