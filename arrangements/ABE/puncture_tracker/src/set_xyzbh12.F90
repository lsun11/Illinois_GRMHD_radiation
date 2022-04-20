!---------------------------------------------------------------------
! Update xbh1, ybh1, zbh1, xbh2, ybh2, and zbh2 with latest positions
!---------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_xyzbh12(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if(cctk_iteration.gt.0) then
     
     xbh1 = puncture_positions(1,2,1)
     ybh1 = puncture_positions(1,2,2)
     zbh1 = puncture_positions(1,2,3)
     
     xbh2 = puncture_positions(2,2,1)
     ybh2 = puncture_positions(2,2,2)
     zbh2 = puncture_positions(2,2,3)

  else
     ! Initial data should set up xbh1, ybh1, zbh1, xbh2, ybh2, and zbh2
     
     puncture_positions(1,2,1) = xbh1
     puncture_positions(1,2,2) = ybh1
     puncture_positions(1,2,3) = zbh1

     puncture_positions(2,2,1) = xbh2
     puncture_positions(2,2,2) = ybh2
     puncture_positions(2,2,3) = zbh2

     puncture_positions_p = puncture_positions
  end if

end subroutine set_xyzbh12
