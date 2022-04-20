#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Setup movingbox center.  
!----------------------------------------------------------------------------------
subroutine setup_movingbox_centers_part1(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
!
  if (x0_BH1 .lt. -9.99d307) then 
     write(*,*) 'You need to set parameter movingbox_for_bbh::x0_BH1'
     stop
  end if

  if (x0_BH2 .lt. -9.99d307) then
     write(*,*) 'You need to set parameter movingbox_for_bbh::x0_BH2'
     stop
  end if

  if (max_dev1 .gt. 9.99d307) then 
     write(*,*) 'You need to set parameter movingbox_for_bbh::max_dev1'
     stop
  end if

  if (max_dev2 .gt. 9.99d307) then
     write(*,*) 'You need to set parameter movingbox_for_bbh::max_dev2'
     stop
  end if

  xc1_prev = x0_BH1
  yc1_prev = 0.d0
  zc1_prev = 0.d0

  xc2_prev = x0_BH2
  yc2_prev = 0.d0
  zc2_prev = 0.d0

  position_x(1) = x0_BH1
  position_y(1) = 0.d0
  position_z(1) = 0.d0
  position_x(2) = x0_BH2
  position_y(2) = 0.d0
  position_z(2) = 0.d0

end subroutine setup_movingbox_centers_part1
