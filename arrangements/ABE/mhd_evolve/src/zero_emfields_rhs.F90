!-----------------------------------------------
! Zero out Bi_rhs's (for em_evolve_ct_active=0)
!-----------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine zero_emfields_rhs(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  Bxtilde_or_Ax_rhs = 0.D0
  Bytilde_or_Ay_rhs = 0.D0
  Bztilde_or_Az_rhs = 0.D0

end subroutine zero_emfields_rhs
