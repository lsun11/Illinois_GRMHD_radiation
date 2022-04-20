#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Update center of a moving box.  
!----------------------------------------------------------------------------------
subroutine movingbox_update_centers_part1(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: index,handle,ierr
  real*8 :: delta
!
  find_max_phi1 = 0
  find_max_phi2 = 0
  if (mod(cctk_iteration,out_every)==0 .and. cctk_iteration .gt. 0) then 
     !!delta = abs(sf_origin_x(1)-xc1_prev) + abs(sf_origin_y(1)-yc1_prev) + abs(sf_origin_z(1)-zc1_prev)

     !!if (delta .lt. 1.d-10) find_max_phi1 = 1

     !!delta = abs(sf_origin_x(2)-xc2_prev) + abs(sf_origin_y(2)-yc2_prev) + abs(sf_origin_z(2)-zc2_prev)

     !!if (delta .lt. 1.d-10) find_max_phi2 = 1

     delta = abs(sf_origin_x(3)-xc1_prev) + abs(sf_origin_y(3)-yc1_prev) + abs(sf_origin_z(3)-zc1_prev)

     if (delta .lt. 1.d-10) find_max_phi1 = 1

  end if

end subroutine movingbox_update_centers_part1
