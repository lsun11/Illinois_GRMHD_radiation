#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_setup_emfields_part4_local_rescaleB_set_conservs(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !write(*,*) "INSIDE driver_setup_emfields_part4_local_rescaleB_set_conservs,dx=",CCTK_DELTA_SPACE(1)

  ! This routine does nothing anymore.

end subroutine driver_setup_emfields_part4_local_rescaleB_set_conservs
