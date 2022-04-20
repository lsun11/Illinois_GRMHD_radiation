#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine update_movingbox(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (enable_moving_grid==1) then
     position_x(1) = 0.96d0*sin(cctk_time*0.5d0)
     position_y(1) = 1.3d0*(1.d0-cos(cctk_time*0.5d0))
     position_z(1) = 0.d0
     write(*,*) 'Moving box center: ',position_x(1),position_y(1),position_z(1)
  end if

end subroutine update_movingbox
