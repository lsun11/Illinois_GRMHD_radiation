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
  integer :: vindex,handle,ierr

  d3x = cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)

  if(cctk_iteration==0) then
     if(WhichIntegral==(1352)) then
        Box1X_VolInt=0.D0
        Box1Y_VolInt=0.D0
        Box1Z_VolInt=0.D0
        Box1denom_VolInt=1.D0
     end if

     if(WhichIntegral==(1353)) then
        Box2X_VolInt=0.D0
        Box2Y_VolInt=0.D0
        Box2Z_VolInt=0.D0
        Box2denom_VolInt=1.D0
     end if

     if(WhichIntegral==(1354)) then
        Box3X_VolInt=0.D0
        Box3Y_VolInt=0.D0
        Box3Z_VolInt=0.D0
        Box3denom_VolInt=1.D0
     end if

     if(WhichIntegral==(1355)) then
        Box4X_VolInt=0.D0
        Box4Y_VolInt=0.D0
        Box4Z_VolInt=0.D0
        Box4denom_VolInt=1.D0
     end if
  end if

  write(*,*) "box1: ",Box1X_VolInt,Box1Y_VolInt,Box1Z_VolInt,Box1denom_VolInt

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
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

     !WhichIntegral = 1351+current_center_to_track
     if(WhichIntegral==(1352)) then
        Box1X_VolInt = VolIntegralx*d3x
        Box1Y_VolInt = VolIntegraly*d3x
        Box1Z_VolInt = VolIntegralz*d3x
        Box1denom_VolInt = VolIntegraldenom*d3x
        write(*,*) "Box1 = ",Box1X_VolInt,Box1Y_VolInt,Box1Z_VolInt
        write(*,*) "Box1 coords = ",Box1X_VolInt/Box1denom_VolInt,Box1Y_VolInt/Box1denom_VolInt,Box1Z_VolInt/Box1denom_VolInt
     else if(WhichIntegral==1353) then
        Box2X_VolInt = VolIntegralx*d3x
        Box2Y_VolInt = VolIntegraly*d3x
        Box2Z_VolInt = VolIntegralz*d3x
        Box2denom_VolInt = VolIntegraldenom*d3x
        write(*,*) "Box2 coords = ",Box2X_VolInt/Box2denom_VolInt,Box2Y_VolInt/Box2denom_VolInt,Box2Z_VolInt/Box2denom_VolInt
     else if(WhichIntegral==1354) then
        Box3X_VolInt = VolIntegralx*d3x
        Box3Y_VolInt = VolIntegraly*d3x
        Box3Z_VolInt = VolIntegralz*d3x
        Box3denom_VolInt = VolIntegraldenom*d3x
        write(*,*) "Box3 coords = ",Box3X_VolInt/Box3denom_VolInt,Box3Y_VolInt/Box3denom_VolInt,Box3Z_VolInt/Box3denom_VolInt
     else if(WhichIntegral==1355) then
        Box4X_VolInt = VolIntegralx*d3x
        Box4Y_VolInt = VolIntegraly*d3x
        Box4Z_VolInt = VolIntegralz*d3x
        Box4denom_VolInt = VolIntegraldenom*d3x
        write(*,*) "Box4 coords = ",Box4X_VolInt/Box4denom_VolInt,Box4Y_VolInt/Box4denom_VolInt,Box4Z_VolInt/Box4denom_VolInt
     end if

     current_center_to_track=current_center_to_track+1
  end if

end subroutine Integrate_vol_integrand_movingbox
