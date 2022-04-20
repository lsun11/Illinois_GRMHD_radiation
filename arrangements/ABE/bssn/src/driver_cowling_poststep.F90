!-----------------------------------------------------------------------------
! Cowling approximation: set rhs's == 0
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine cowling_poststep()

  implicit none

  ! Do nothing here.
  ! Just want Carpet to do prolongation for BSSN_refbd

end subroutine cowling_poststep
