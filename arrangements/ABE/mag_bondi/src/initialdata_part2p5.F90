#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mag_bondi_initialdata_part2p5(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: ierr,handle,vindex

  ! Find xmax and zmax 
  call CCTK_ReductionHandle(handle,"maximum")

  call CCTK_VarIndex(vindex,"grid::Z")
  call CCTK_Reduce(ierr,cctkGH, -1,handle, 1, &
       CCTK_VARIABLE_REAL,zmax_bondi,1,vindex)

  call CCTK_VarIndex(vindex,"grid::X")
  call CCTK_Reduce(ierr,cctkGH, -1,handle, 1, &
       CCTK_VARIABLE_REAL,xmax_bondi,1,vindex)

  write(*,*), "HI XMAX_BONDI, ZMAX_BONDI:",xmax_bondi,zmax_bondi

end subroutine mag_bondi_initialdata_part2p5
