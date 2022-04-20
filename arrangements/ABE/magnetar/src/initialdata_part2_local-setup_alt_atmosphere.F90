!--------------------------------------------------------
! Okay, we've read in the initial data from files.
! Now we set up all other required variables, including:
!  emfields, BSSN variables, primitives, etc.
!--------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine magnetar_setup_alt_atmosphere(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! The alternative atmosphere scheme does not depend on rho_b_max, and so can be defined here.
  if(enable_alt_atmosphere==1) then
     ! When enable_alt_atmosphere==0 (default), the atmosphere will be set up in compute_magnetar!
     call Setup_alt_Atmosphere(CCTK_PASS_FTOF)
  end if

end subroutine magnetar_setup_alt_atmosphere
