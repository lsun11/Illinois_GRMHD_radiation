#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!|---------------------------------------------------------------------+
!|
!| Function to compute ADM Mass surface integral centered at (xcenter,ycenter,zcenter)
!|
!|---------------------------------------------------------------------+
subroutine M_surf_integral_offcenter_nofish(cctkGH,xcenter,ycenter,zcenter,output_integral)
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! Variables in the function call
  CCTK_POINTER                             :: cctkGH
  real*8                                   :: xcenter,ycenter,zcenter
  real*8                                   :: output_integral
  ! Variables needed for interpolation
  CCTK_POINTER, dimension(3)               :: interp_coords
  character(60)                            :: options_string
  integer                                  :: interpolation_order,nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle,N_interp_points
  integer                                  :: N_input_arrays,N_output_arrays
  real*8, dimension(N_theta*N_phi)         :: xinterp,yinterp,zinterp
  integer,dimension(13)                    :: input_array_type_codes,input_array_varindices,output_array_type_codes
  CCTK_POINTER,dimension(13)               :: output_array_pointers
  ! Output arrays, dummy indices, parameters
  real*8,dimension(N_theta*N_phi)          :: Psiint,gupxxint,gupxyint,gupxzint,gupyyint,gupyzint,gupzzint
  real*8,dimension(N_theta*N_phi)          :: phixint,phiyint,phizint,Gamxint,Gamyint,Gamzint
  real*8                                   :: phiangle,costheta,sintheta,PI,nxp,nyp,nzp
  integer                                  :: i,j,n

  PI = 3.14159265358979323844D0

  N_dims = 3
  interpolation_order = 2
  N_interp_points = N_theta*N_phi

  interp_handle = -1 
  call CCTK_InterpHandle (interp_handle, "Lagrange polynomial interpolation")
  !  call CCTK_InterpHandle (interp_handle, "uniform cartesian")
  if (interp_handle .lt. 0) then
     call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
  endif

  param_table_handle = -1
  options_string = "order = " // char(ichar('0') + interpolation_order)
  call Util_TableCreateFromString (param_table_handle, options_string)
  if (param_table_handle .lt. 0) then
     call CCTK_WARN(0,"Cannot create parameter table for interpolator")
  endif

  coord_system_handle = -1
  call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
  if (coord_system_handle .lt. 0) then
     call CCTK_WARN(0,"Cannot get handle for cart3d coordinate system ! Forgot to activate an implementation providing coordinates ??")
  endif

  input_array_type_codes = CCTK_VARIABLE_REAL
  output_array_type_codes = CCTK_VARIABLE_REAL

  ! Specify interpolation input arrays, output arrays:
  N_input_arrays = 13
  N_output_arrays = 13

  call CCTK_VarIndex (input_array_varindices(1), "bssn::phi")
  call CCTK_VarIndex (input_array_varindices(2), "bssn::gupxx")
  call CCTK_VarIndex (input_array_varindices(3), "bssn::gupxy")
  call CCTK_VarIndex (input_array_varindices(4), "bssn::gupxz")
  call CCTK_VarIndex (input_array_varindices(5), "bssn::gupyy")
  call CCTK_VarIndex (input_array_varindices(6), "bssn::gupyz")
  call CCTK_VarIndex (input_array_varindices(7), "bssn::gupzz")

  call CCTK_VarIndex (input_array_varindices(8), "bssn::phix")
  call CCTK_VarIndex (input_array_varindices(9), "bssn::phiy")
  call CCTK_VarIndex (input_array_varindices(10), "bssn::phiz")
  call CCTK_VarIndex (input_array_varindices(11), "bssn::Gammax")
  call CCTK_VarIndex (input_array_varindices(12), "bssn::Gammay")
  call CCTK_VarIndex (input_array_varindices(13), "bssn::Gammaz")

  output_array_pointers(1) = CCTK_PointerTo(Psiint)
  output_array_pointers(2) = CCTK_PointerTo(gupxxint)
  output_array_pointers(3) = CCTK_PointerTo(gupxyint)
  output_array_pointers(4) = CCTK_PointerTo(gupxzint)
  output_array_pointers(5) = CCTK_PointerTo(gupyyint)
  output_array_pointers(6) = CCTK_PointerTo(gupyzint)
  output_array_pointers(7) = CCTK_PointerTo(gupzzint)
  output_array_pointers(8) = CCTK_PointerTo(phixint)
  output_array_pointers(9) = CCTK_PointerTo(phiyint)
  output_array_pointers(10) = CCTK_PointerTo(phizint)
  output_array_pointers(11) = CCTK_PointerTo(Gamxint)
  output_array_pointers(12) = CCTK_PointerTo(Gamyint)
  output_array_pointers(13) = CCTK_PointerTo(Gamzint)

  ! Set up interpolation coordinate arrays:
  n = 1
  do i=1,N_theta
     costheta = 1.D0 - (i - 0.5D0)*dcostheta
     sintheta = sqrt(1.D0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5D0)*dphi
        if(N_phi==1) phiangle = 0.D0
        xinterp(n) = surf_radius*sintheta*cos(phiangle)+xcenter
        yinterp(n) = surf_radius*sintheta*sin(phiangle)+ycenter
        zinterp(n) = surf_radius*costheta+zcenter

        n = n + 1
     end do
  end do

!  write(*,*) "surface_integrals-M_surf_NOfish_offcenter.F90: center:",xcenter,ycenter,zcenter,surf_radius

  interp_coords(1) = CCTK_PointerTo(xinterp)
  interp_coords(2) = CCTK_PointerTo(yinterp)
  interp_coords(3) = CCTK_PointerTo(zinterp)

  ! Perform interpolation:
  call CCTK_InterpGridArrays(ierr,cctkGH,N_dims,interp_handle, &
       param_table_handle,coord_system_handle, &
       N_interp_points,input_array_type_codes,interp_coords, &
       N_input_arrays,input_array_varindices, &
       N_output_arrays, output_array_type_codes, output_array_pointers)

  Psiint = exp(Psiint)
  phixint = phixint*Psiint
  phiyint = phiyint*Psiint
  phizint = phizint*Psiint

  !|~~~~~> Integrate the dot product of Dphi with the surface normal
  if(ierr.ge.0) then
     output_integral = 0.D0
     do i=1,ntot
        nxp = (xinterp(i)-xcenter)/surf_radius
        nyp = (yinterp(i)-ycenter)/surf_radius
        nzp = (zinterp(i)-zcenter)/surf_radius
        
        output_integral = output_integral + (nxp*( 0.125D0*Gamxint(i) - (phixint(i)*gupxxint(i) + phiyint(i)*gupxyint(i) + phizint(i)*gupxzint(i)) ) + &
             nyp*( 0.125D0*Gamyint(i) - (phixint(i)*gupxyint(i) + phiyint(i)*gupyyint(i) + phizint(i)*gupyzint(i)) ) + &
             nzp*( 0.125D0*Gamzint(i) - (phixint(i)*gupxzint(i) + phiyint(i)*gupyzint(i) + phizint(i)*gupzzint(i)) ))
     end do

     !|  Multiply by 1/2\pi * surface element.
     output_integral = output_integral * surf_radius * surf_radius * dphi * dcostheta * 0.5D0 / PI*sym_factor
  else
     output_integral = -9999.D0
  end if

end subroutine M_surf_integral_offcenter_nofish
