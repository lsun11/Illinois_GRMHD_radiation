!------------------------------------------------
! Master timestepping driver for the shift thorn
!------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine shift_update_bc(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)				:: ext
  integer                                       :: dummy
  real*8                                        :: dX,dY,dZ
!
  integer                                       :: i,j,k
  real*8                                        :: rfromcenterofAH,shiftmag,AH_radius_minimum,horizdirn_x,horizdirn_y,horizdirn_z

  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)

 ext = cctk_lsh

!write(*,*) "Start shift update_bc!!!!!!"

  call CartSymGN(dummy,cctkGH,'shift::shift_vars')

  if(cctk_iteration.gt.0) then
     !cowling_enable=1 disables bssn+lapse+shift
     if (Spatial_Gauge.eq.0 .or. Spatial_Gauge.eq.3 .or. cowling_enable.eq.1) then
        !Do nothing; boundary already updated
     else if (Spatial_Gauge.eq.1) then
        call hyperbolic_shift_bc_driver(CCTK_PASS_FTOF)
     else if (Spatial_Gauge.eq.5 .or. Spatial_Gauge.eq.7) then
        call hbpuncture_shift_bc_driver(CCTK_PASS_FTOF)

     else if (Spatial_Gauge.eq.6) then
        call puncture_firstorder_shift_bc_driver(CCTK_PASS_FTOF)
     end if
  end if

  call CartSymGN(dummy,cctkGH,'shift::shift_vars')

end subroutine shift_update_bc
