#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!--------------------------------------------------
! Perform volume integration for matter quantities
! Note:
! When CCTK_Reduce() is called in a Carpet simulation, it actually sums the integrand*weight function.
!    The weight function is set to 1/(cctk_levfac(1)*cctk_levfac(2)*cctk_levfac(3)) in non-
!    ghostzone regions of the grid, where e.g., cctk_levfac(1)=1 on the coarsest level and
!    it increases from there by factors of 2 as the grids get finer.
!    This is why we multiply by cctk_delta_space's (the coarsest gridspacing).
!--------------------------------------------------
subroutine Integrate_vol_integrand_shock(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8  :: VolIntegral,dx
  integer :: vindex,handle,ierr

  if(MOD(cctk_iteration,out_every)==0) then
     dx = cctk_delta_space(1)

     call CCTK_ReductionHandle(handle,"sum")

     if(Shock_Which_Int==1000) then
        call CCTK_VarIndex(vindex,"shocktests::rho_error")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        rho_error_VolInt = VolIntegral*dx

        call CCTK_VarIndex(vindex,"shocktests::P_error")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        P_error_VolInt = VolIntegral*dx

        call CCTK_VarIndex(vindex,"shocktests::vx_error")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        vx_error_VolInt = VolIntegral*dx

        call CCTK_VarIndex(vindex,"shocktests::rho_error_reflection")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        rho_error_VolInt_reflection = VolIntegral*dx

        call CCTK_VarIndex(vindex,"shocktests::P_error_reflection")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        P_error_VolInt_reflection = VolIntegral*dx

        call CCTK_VarIndex(vindex,"shocktests::vx_error_reflection")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        vx_error_VolInt_reflection = VolIntegral*dx
     end if

     if(Shock_Which_Int==1001) then
        call CCTK_VarIndex(vindex,"shocktests::shock_rest_mass_integrand")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        shock_restmass_VolInt = VolIntegral*dx

        call CCTK_VarIndex(vindex,"shocktests::shock_restmass_analytic_integrand")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
        shock_restmass_analytic = VolIntegral*dx

        if(cctk_iteration==0) then
           !Note that this analytic value for the rest mass assumes ZERO initial velocity everywhere!
           call CCTK_VarIndex(vindex,"shocktests::shock_restmass_analytic_initial_integrand")
           call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)
           shock_restmass_analytic_initial = VolIntegral*dx
        end if

     end if
  end if
end subroutine Integrate_vol_integrand_shock
