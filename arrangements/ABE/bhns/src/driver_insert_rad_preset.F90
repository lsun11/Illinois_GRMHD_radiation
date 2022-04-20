#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine BHNS_insert_radiation_preset(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                  :: handle,vindex,ierr
  real*8                                   :: reduction_value

  bhnsinsertRADNOW = 0

if(CCTK_ITERATION .eq. iteration_to_insert_rad .and. iteration_to_insert_rad.gt.0) then
   print *, "CCTK_ITERATION=",CCTK_ITERATION, "iteration_to_insert_rad=",iteration_to_insert_rad
  bhnsinsertRADNOW=1
end if

end subroutine BHNS_insert_radiation_preset
