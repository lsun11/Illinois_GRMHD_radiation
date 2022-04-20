#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Update center of a moving box.  
!----------------------------------------------------------------------------------
subroutine movingbox_update_center_part1(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: index,handle,ierr
  real*8 :: delta
!
  find_max_phi = 0
  if (mod(cctk_iteration,out_every)==0 .and. cctk_iteration .gt. 0) then 
     delta = abs(sf_origin_x(1)-xc_prev) + abs(sf_origin_y(1)-yc_prev) + abs(sf_origin_z(1)-zc_prev)

     if (delta .lt. 1.d-10) find_max_phi = 1

     if (find_max_phi == 1) then 
        call CCTK_VarIndex(index,"bssn::phi")
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,max_phi,1,index)
     end if

  end if

end subroutine movingbox_update_center_part1
