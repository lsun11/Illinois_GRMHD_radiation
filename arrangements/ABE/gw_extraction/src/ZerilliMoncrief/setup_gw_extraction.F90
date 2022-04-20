!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
!
! Function to compute GW_radius_phys
!
!-----------------------------------------------------------------------------
subroutine setup_gw_extraction(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8, dimension(3) :: pointcoords
  real*8, dimension(1) :: outputblah
  integer :: ierr,handle,vindex

  if(radius_GW .lt. 0.D0) then
     radius_GW = surf_radius
  end if

  ! FOLLOWING FOR GW EXTRACTION CODE:
  if(fisheye_enable==1) then
     if(radius_GW_phys .lt. 0.D0) then

        pointcoords(1) = radius_GW*sin(theta_GW)*cos(phi_GW)
        pointcoords(2) = radius_GW*sin(theta_GW)*sin(phi_GW)
        pointcoords(3) = radius_GW*cos(theta_GW)

        call CCTK_VarIndex(vindex,"fisheye::PhysicalRadius")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,radius_GW_phys)
        call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,dR_GW)
        call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative2")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,ddR_GW)

     end if
  else 
     radius_GW_phys = radius_GW
  end if
end subroutine setup_gw_extraction
