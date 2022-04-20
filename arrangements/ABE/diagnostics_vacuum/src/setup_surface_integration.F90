!-----------------------------------------------------------------------------
!
!$Id: setup_surface_integration.F90  $
!
!-----------------------------------------------------------------------------
!
! setup various quantities for surface integrals
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine setup_surface_integration(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  real*8                                   :: dX,dY,dZ
  real*8, dimension(1,3)                   :: pointcoords
  real*8, dimension(1)                     :: interp_output
  real*8                                   :: PI,costheta,sintheta,phiangle
  real*8                                   :: zmax
  integer                                  :: handle,ierr,vindex
  integer                                  :: i,n,j,int_order
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  !If surface radius is set to negative number (as defaulted in param.ccl), set it to following:
  if(surf_radius.lt.0.D0) then
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax,1,vindex)

     surf_radius = -surf_radius * (0.5D0*dZ+zmax)
     scaledsurf=1

     write(*,*) "OISDFO:",ierr,zmax,surf_radius

     if (fisheye_enable==1) then
        pointcoords(1,1)=surf_radius
        pointcoords(1,2)=0.d0
        pointcoords(1,3)=0.d0

        call CCTK_VarIndex(vindex,"fisheye::PhysicalRadius")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,interp_output)
        rbr = interp_output(1)/surf_radius
        call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,interp_output)
        drbdr = interp_output(1)
        call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative2")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,interp_output)
        ddrbddr = interp_output(1)
        write(6,*)'radii:',surf_radius,rbr*surf_radius,rbr,drbdr,ddrbddr
     else
        rbr = 1.d0
        drbdr = 1.d0
        ddrbddr = 0.d0
     end if
  end if

  if(cctk_iteration==0) then

     dX = CCTK_DELTA_SPACE(1)
     dY = CCTK_DELTA_SPACE(2)
     dZ = CCTK_DELTA_SPACE(3)

     ext = cctk_lsh
     global_ext = cctk_gsh

     PI = acos(-1.D0)

     int_order=2

     N_theta = numtheta
     N_phi = numphi
     dphi=2.D0 * PI / numphi
     dcostheta = 2.D0 / numtheta

     if (Symmetry==OCTANT) then 
        N_theta = numtheta/2
        N_phi = numphi/4
        sym_factor=8
     else if (Symmetry==EQUATORIAL) then
        N_theta = numtheta/2
        sym_factor=2
     else if (Symmetry==NO_SYMM) then 
        sym_factor=1
     else if (Symmetry==PI_SYMM) then 
        N_theta = numtheta/2
        N_phi = numphi/2
        sym_factor=4
     else if (Symmetry==AXISYM) then
        N_theta = numtheta/2
        N_phi = 1
        dphi=2.0*PI
        sym_factor=1
     end if

     ntot = N_theta*N_phi

     if (fisheye_enable==1) then
        pointcoords(1,1)=surf_radius
        pointcoords(1,2)=0.d0
        pointcoords(1,3)=0.d0

        call CCTK_VarIndex(vindex,"fisheye::PhysicalRadius")     
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,interp_output)
        rbr = interp_output(1)/surf_radius
        call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,interp_output)
        drbdr = interp_output(1)
        call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative2")
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,interp_output)
        ddrbddr = interp_output(1)
        write(6,*)'radii:',surf_radius,rbr*surf_radius,rbr,drbdr,ddrbddr
     else
        rbr = 1.d0
        drbdr = 1.d0
        ddrbddr = 0.d0
     end if
     write(6,*)'rbr:',rbr,drbdr,ddrbddr

  end if

end subroutine setup_surface_integration
