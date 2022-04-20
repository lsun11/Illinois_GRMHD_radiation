#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for bhns thorn
!-----------------------------------------------------------------------------
subroutine bhns_diagnostics_global0(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                  :: ierr,handle,vindex

  if(MOD(cctk_iteration,out_every)==0) then
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex (vindex, "mhd_evolve::temp8")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,bhns_max_b2,1,vindex)
  end if
end subroutine bhns_diagnostics_global0
