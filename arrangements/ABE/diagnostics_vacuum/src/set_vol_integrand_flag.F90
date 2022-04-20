!-----------------------------------------------------------------------------
! Sets volume integration flag (Compute_Vol_Integrands)
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_vol_integrand_flag(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     Compute_Vol_Integrands = 1
  else
     Compute_Vol_Integrands = 0
  end if

end subroutine set_vol_integrand_flag
