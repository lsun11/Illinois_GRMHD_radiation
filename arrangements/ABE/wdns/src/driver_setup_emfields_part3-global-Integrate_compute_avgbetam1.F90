#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!--------------------------------------------------
! Perform volume integration for matter quantities
! Note:
! When a CCTK_Reduce() sum is called in a Carpet simulation, it actually sums the integrand*weight function.
!    The weight function is set to 1/(cctk_levfac(1)*cctk_levfac(2)*cctk_levfac(3)) in non-
!    ghostzone regions of the grid, where e.g., cctk_levfac(1)=1 on the coarsest level and
!    it increases from there by factors of 2 as the grids get finer.
!    This is why we multiply by cctk_delta_space's (the coarsest gridspacing).
!--------------------------------------------------
subroutine Integrate_vol_integrand_wdns_emfields_compute_avgbetam1(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8  :: VolIntegral,d3x,wdns_pmV,wdns_pV
  integer :: vindex,handle,ierr

  if(CCTK_ITERATION .eq. ITERATION_TO_RESET_MAGNETIC_FIELDS .or. CCTK_ITERATION .eq. 0) then
     !This routine is only a diagnostic.

     d3x = cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)

     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp1")
     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)

     wdns_pmV = VolIntegral*d3x

     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp2")
     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegral, 1, vindex)

     wdns_pV = VolIntegral*d3x

     wdns_avg_betam1 = wdns_pmV/wdns_pV
     write(*,*) "INSIDE Integrate_vol_integrand_wdns_emfields_compute_avgbetam1."
     write(*,*) "WDNS: AVG BETAM1 [Integral b^2/2 over proper volume] / [Integral P over proper volume]",wdns_avg_betam1
     write(*,*) "WDNS: [Integral b^2/2 over proper volume]",wdns_pmV
     write(*,*) "WDNS: [Integral P over proper volume]",wdns_pV
     ! WDNS: AVG BETAM1    5197.047870508249         22.55140879593439      
     !   4.3392728637169477E-003
  end if

end subroutine Integrate_vol_integrand_wdns_emfields_compute_avgbetam1