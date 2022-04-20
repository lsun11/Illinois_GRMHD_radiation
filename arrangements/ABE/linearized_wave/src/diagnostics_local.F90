#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!------------------------------------------------------
! Linearized wave diagnostics: compute local quantities
!------------------------------------------------------
subroutine lw_diagnostics_local(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext

  ext = cctk_lsh

  if(MOD(cctk_iteration,out_every)==0) then
     call lin_wave_analytic_gxx(ext, X, Y, Z, cctk_time+time_shift, amplitude, width, PsiTau, &
          PhysicalRadius,RadiusDerivative, mode)
     !PsiTau = abs(PsiTau - gxx)
  end if

end subroutine lw_diagnostics_local
