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
subroutine Integrate_vol_integrand_movingbox(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8  :: VolIntegralx,VolIntegraly,VolIntegralz,VolIntegraldenom,d3x
  real*8  :: VolIntegralx1,VolIntegraly1,VolIntegralz1,VolIntegraldenom1
  real*8  :: VolIntegralx2,VolIntegraly2,VolIntegralz2,VolIntegraldenom2
  integer :: vindex,handle,ierr

  d3x = cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)

  if(track_bhns.eq.0) then
     
     if(cctk_iteration==0) then
        if(WhichIntegral==(1351)) then
           Box1X_VolInt=0.D0
           Box1Y_VolInt=0.D0
           Box1Z_VolInt=0.D0
           Box1denom_VolInt=1.D0
        end if
     end if
     write(*,*) "box1,radius: ",Box1X_VolInt,Box1Y_VolInt,Box1Z_VolInt,Box1denom_VolInt,radius(1,1)
        
     if(mod(cctk_iteration,Compute_VolIntegrands_Every)==0) then
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"movingbox::box_x")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegralx, 1, vindex)
        
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"movingbox::box_y")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegraly, 1, vindex)
        
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"movingbox::box_z")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegralz, 1, vindex)
        
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"movingbox::box_denom")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegraldenom, 1, vindex)
        
        !WhichIntegral is set to 1351 in vol_integrands-movingbox.F90
        if(WhichIntegral==(1351)) then
           Box1X_VolInt = VolIntegralx*d3x
           Box1Y_VolInt = VolIntegraly*d3x
           Box1Z_VolInt = VolIntegralz*d3x
           Box1denom_VolInt = VolIntegraldenom*d3x
           write(*,*) "Box1 = ",Box1X_VolInt,Box1Y_VolInt,Box1Z_VolInt,VolIntegraldenom
           write(*,*) "Box1 coords = ",Box1X_VolInt/Box1denom_VolInt,Box1Y_VolInt/Box1denom_VolInt,Box1Z_VolInt/Box1denom_VolInt
        end if
     end if
        
  else if(track_bhns.eq.1) then
     
     if (cctk_iteration .eq. 0) then
        Box1denom_VolInt1 = 0.d0
        Box1denom_VolInt2 = 0.d0
     end if
     
     if(mod(cctk_iteration,Compute_VolIntegrands_Every)==0 .and. cctk_iteration .gt. 1) then	
        !  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"movingbox::box_x")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegralx1, 1, vindex)
        
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"movingbox::box_y")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegraly1, 1, vindex)
        
!     call CCTK_ReductionHandle(handle,"sum")
!     call CCTK_VarIndex(vindex,"movingbox::box_z")
!     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegralz1, 1, vindex)

        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"movingbox::box_denom")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegraldenom1, 1, vindex)
        

        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"mhd_evolve::temp1")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegralx2, 1, vindex)
        
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"mhd_evolve::temp2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegraly2, 1, vindex)

!     call CCTK_ReductionHandle(handle,"sum")
!     call CCTK_VarIndex(vindex,"mhd_evolve::temp3")
!     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegralz2, 1, vindex)

        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"mhd_evolve::temp4")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, VolIntegraldenom2, 1, vindex)
        

     !WhichIntegral is set to 1351 in vol_integrands-movingbox.F90
!     if(WhichIntegral==(1351)) then
        Box1X_VolInt1 = VolIntegralx1*d3x
        Box1Y_VolInt1 = VolIntegraly1*d3x
!        Box1Z_VolInt1 = VolIntegralz2*d3x
        Box1denom_VolInt1 = VolIntegraldenom1*d3x

        Box1X_VolInt2 = VolIntegralx2*d3x
        Box1Y_VolInt2 = VolIntegraly2*d3x
!        Box1Z_VolInt2 = VolIntegralz2*d3x
        Box1denom_VolInt2 = VolIntegraldenom2*d3x


        write(*,*) "Box1 = ",Box1X_VolInt1,Box1Y_VolInt1,VolIntegraldenom1
        write(*,*) "Box1 coords = ",Box1X_VolInt1/Box1denom_VolInt1,Box1Y_VolInt1/Box1denom_VolInt1
        
        write(*,*) "Box2 = ",Box1X_VolInt2,Box1Y_VolInt2,VolIntegraldenom2
        write(*,*) "Box2 coords = ",Box1X_VolInt2/Box1denom_VolInt2,Box1Y_VolInt2/Box1denom_VolInt2
     end if
  else
     write(*,*) "ERROR, UNSUPPORTED: System="
  end if
  
end subroutine Integrate_vol_integrand_movingbox
