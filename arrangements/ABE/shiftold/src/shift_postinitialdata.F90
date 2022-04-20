!----------------------------------------------------------------
! Stuff that needs to be done right after initial data is set up
!----------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine shift_PostInitialData(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 					:: dT, dX,dY,dZ
  integer, dimension(3)		 	        :: ext
  integer					:: dummy,index
!
  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !=================================
  ! HERE WE COPY TO PREV TIMELEVEL:
  shiftx_p = shiftx
  shifty_p = shifty
  shiftz_p = shiftz

  shiftxt_p = shiftxt
  shiftyt_p = shiftyt
  shiftzt_p = shiftzt

  shiftx_p_p = shiftx
  shifty_p_p = shifty
  shiftz_p_p = shiftz

  shiftxt_p_p = shiftxt
  shiftyt_p_p = shiftyt
  shiftzt_p_p = shiftzt


  ! Set rhs's to Zero!
  shiftx_rhs = 0.D0
  shifty_rhs = 0.D0
  shiftz_rhs = 0.D0
  !=================================

end subroutine Shift_PostInitialData
