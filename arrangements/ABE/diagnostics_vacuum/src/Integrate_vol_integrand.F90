#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!----------------------------
! Perform volume integration
! Note:
! When a CCTK_Reduce() sum is called in a Carpet simulation, it actually sums the integrand*weight function.
!    The weight function is set to 1/(cctk_levfac(1)*cctk_levfac(2)*cctk_levfac(3)) in non-
!    ghostzone regions of the grid, where e.g., cctk_levfac(1)=1 on the coarsest level and
!    it increases from there by factors of 2 as the grids get finer.
!    This is why we multiply by cctk_delta_space's (the coarsest gridspacing).
!----------------------------
subroutine Integrate_vol_integrand(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8  :: VolIntegral,d3x
  integer :: vindex,handle,ierr

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     d3x = cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)

     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand")
     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)

     if(WhichIntegral==1) then
        M_ADM_VolInt = VolIntegral*d3x
     else if(WhichIntegral==2) then
        J_ADM_VolInt = VolIntegral*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        J_ADM_VolInt_inner = VolIntegral*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand3")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        J_ADM_VolInt_inner2 = VolIntegral*d3x

        print *, "value of whichintegral",WhichIntegral

     else if(WhichIntegral==3) then
        Ham_const_VolIntN = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        Ham_const_VolIntD = sqrt(VolIntegral)*d3x
     else if(WhichIntegral==4) then
        Ham_const_excised_VolIntN = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        Ham_const_excised_VolIntD = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand3")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        Ham_const_excised_innerregion_VolIntN = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand4")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        Ham_const_excised_innerregion_VolIntD = sqrt(VolIntegral)*d3x
     else if(WhichIntegral==5) then
        momx_const_VolIntN = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        momy_const_VolIntN = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand3")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        momz_const_VolIntN = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand4")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        mom_const_VolIntD = sqrt(VolIntegral)*d3x
     else if(WhichIntegral==6) then
        momx_const_excised_VolIntN = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        momy_const_excised_VolIntN = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand3")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        momz_const_excised_VolIntN = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand4")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        mom_const_excised_VolIntD = sqrt(VolIntegral)*d3x
     else if(WhichIntegral==7) then
        Gamx_const_VolInt = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        Gamy_const_VolInt = sqrt(VolIntegral)*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand3")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        Gamz_const_VolInt = sqrt(VolIntegral)*d3x
     else if(WhichIntegral==8) then
        M_constraint = VolIntegral*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        M_constraint_inner = VolIntegral*d3x
     else if(WhichIntegral==9) then
        Jz_constraint = VolIntegral*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        Jz_constraint_inner = VolIntegral*d3x
     else if(WhichIntegral==10) then
        P_constraint = VolIntegral*d3x
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"diagnostics_vacuum::VolIntegrand2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        P_constraint_inner = VolIntegral*d3x
     end if
end if
end subroutine Integrate_vol_integrand
