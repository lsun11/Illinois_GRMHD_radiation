!----------------------------------------------------------------
! Stuff that needs to be done right after initial data is set up
!----------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine lapse_PostInitialData(CCTK_ARGUMENTS)

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
  lapm1_p = lapm1
  lapset_p = lapset

  lapm1_p_p = lapm1
  lapset_p_p = lapset

  ! Set rhs's to Zero!
  lapm1_rhs = 0.D0
  lapset_rhs = 0.D0
  !=================================

  

end subroutine lapse_PostInitialData
