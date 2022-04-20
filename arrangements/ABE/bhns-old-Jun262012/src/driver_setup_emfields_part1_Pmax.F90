#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine BHNS_setup_emfield_part1_Pmax(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                  :: handle,vindex,ierr
  real*8                                   :: reduction_value

  if(CCTK_ITERATION .eq. ITERATION_TO_INSERT_MAGNETIC_FIELDS .or. CCTK_ITERATION .eq. 0) then

     ! Find P_max
     call CCTK_VarIndex(vindex,"mhd_evolve::P")
     call CCTK_ReductionHandle(handle,"maximum")
     if (handle .gt. 0) then
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,vindex)
     else
        call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
     end if
     bhns_P_max = reduction_value

     write(*,*) "INSIDE BHNS_setup_emfield_part1_Pmax.  PMAX = ",bhns_P_max
  end if
end subroutine BHNS_setup_emfield_part1_Pmax
