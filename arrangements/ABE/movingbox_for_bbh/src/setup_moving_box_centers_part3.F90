#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Setup movingbox center.  
!----------------------------------------------------------------------------------
subroutine setup_movingbox_centers_part3(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  integer :: handle,vindex,ierr
  real*8 :: levfac_max1,levfac_max2
!
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(vindex,"movingbox_for_bbh::tempx1")
  call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, levfac_max1, 1, vindex)

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(vindex,"movingbox_for_bbh::tempy1")
  call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, levfac_max2, 1, vindex)

  max_levfac1 = int(levfac_max1+0.1d0)
  max_levfac2 = int(levfac_max2+0.1d0)
  write(*,*) 'max_levfac = ',max_levfac1,max_levfac2

end subroutine setup_movingbox_centers_part3
