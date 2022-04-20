!-----------------------------------------------------------------------------
! A function to initialize the Shift
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine Setup_Initial_Shift(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8   :: ZERO
  parameter (ZERO = 0.D0)
  
  if (Spatial_Gauge.eq.0) then
     shiftx = ZERO
     shifty = ZERO
     shiftz = ZERO

     shiftx_p = ZERO
     shifty_p = ZERO
     shiftz_p = ZERO
  end if

  if (Spatial_Gauge.eq.1.or.Spatial_Gauge.eq.5) then
     shiftxt = ZERO
     shiftyt = ZERO
     shiftzt = ZERO

     shiftxt_p = ZERO
     shiftyt_p = ZERO
     shiftzt_p = ZERO

     Gammax_drive = ZERO
     Gammay_drive = ZERO
     Gammaz_drive = ZERO

     hyper_psi6_init = exp(6.D0*hyper_phi_init)
  end if  
  if (Spatial_Gauge.eq.6) then
     shiftxt_p = ZERO
     shiftyt_p = ZERO
     shiftzt_p = ZERO
     shiftxt_rhs = ZERO
     shiftyt_rhs = ZERO
     shiftzt_rhs = ZERO
     shiftxt = ZERO
     shiftyt = ZERO
     shiftzt = ZERO
  end if
end subroutine Setup_Initial_Shift

