!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
!
! Function to compute EM_radius_phys
!
!-----------------------------------------------------------------------------
subroutine setup_em_extraction(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8, dimension(3) :: pointcoords
  real*8, dimension(1) :: outputblah
  integer :: ierr,handle,vindex

  if(radius_EM .lt. 0.D0) then
     radius_EM = surf_radius
  end if

  ! FOLLOWING FOR EM EXTRACTION CODE:
  if(fisheye_enable==1) then
     if(radius_EM_phys .lt. 0.D0) then

        pointcoords(1) = radius_EM*sin(theta_EM)*cos(phi_EM)
        pointcoords(2) = radius_EM*sin(theta_EM)*sin(phi_EM)
        pointcoords(3) = radius_EM*cos(theta_EM)

        call CCTK_VarIndex(vindex,"fisheye::PhysicalRadius")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,radius_EM_phys)
        call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,dR_EM)
        call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative2")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,ddR_EM)

     end if
  else 
     radius_EM_phys = radius_EM
  end if
end subroutine setup_em_extraction
