!-----------------------------------------------------------------------------
! Compute 1/8\pi \int \psi^6 (x A^m_y - yA^m_x) dS_m.
!  Note that we work around a surface centered at xcenter,ycenter,zcenter
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine J_surf_integral_spin(cctkGH,ntheta,nphi,x_cm,y_cm,z_cm,r_local,output_integralx,output_integraly,output_integralz,symm)
  implicit none
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  ! Variables in the function call
  CCTK_POINTER                             :: cctkGH
  real*8                                   :: x_cm,y_cm,z_cm,r_local
  real*8                                   :: output_integralx,output_integraly,output_integralz

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM,symm
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  ! Variables needed for interpolation
  CCTK_POINTER, dimension(3)               :: interp_coords
  character(60)                            :: options_string
  integer                                  :: interpolation_order,nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle,N_interp_points
  integer                                  :: ntheta,nphi
  integer                                  :: N_input_arrays,N_output_arrays
  real*8, dimension(ntheta*nphi)           :: xinterp,yinterp,zinterp
  integer,dimension(13)                    :: input_array_type_codes,input_array_varindices,output_array_type_codes
  CCTK_POINTER,dimension(13)               :: output_array_pointers
  ! Output arrays, dummy indices, parameters
  real*8,dimension(ntheta*nphi)           :: Psi6int,gupxxint,gupxyint,gupxzint,gupyyint,gupyzint,gupzzint
  real*8,dimension(ntheta*nphi)           :: Axxint,Axyint,Axzint,Ayyint,Ayzint,Azzint
  real*8                                   :: phiangle,costheta,sintheta,PI,nxp,nyp,nzp
  integer                                  :: i,j,n

  PI = 3.14159265358979323846D0

  N_dims = 3
  interpolation_order = 2


  sym_factor = 2.d0
  if (Symm==NO_SYMM) sym_factor = 1.d0


  N_interp_points = ntheta*nphi


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

  call CCTK_VarIndex (input_array_varindices(8), "bssn::Axx")
  call CCTK_VarIndex (input_array_varindices(9), "bssn::Axy")
  call CCTK_VarIndex (input_array_varindices(10), "bssn::Axz")
  call CCTK_VarIndex (input_array_varindices(11), "bssn::Ayy")
  call CCTK_VarIndex (input_array_varindices(12), "bssn::Ayz")
  call CCTK_VarIndex (input_array_varindices(13), "bssn::Azz")

  output_array_pointers(1) = CCTK_PointerTo(Psi6int)
  output_array_pointers(2) = CCTK_PointerTo(gupxxint)
  output_array_pointers(3) = CCTK_PointerTo(gupxyint)
  output_array_pointers(4) = CCTK_PointerTo(gupxzint)
  output_array_pointers(5) = CCTK_PointerTo(gupyyint)
  output_array_pointers(6) = CCTK_PointerTo(gupyzint)
  output_array_pointers(7) = CCTK_PointerTo(gupzzint)
  output_array_pointers(8) = CCTK_PointerTo(Axxint)
  output_array_pointers(9) = CCTK_PointerTo(Axyint)
  output_array_pointers(10) = CCTK_PointerTo(Axzint)
  output_array_pointers(11) = CCTK_PointerTo(Ayyint)
  output_array_pointers(12) = CCTK_PointerTo(Ayzint)
  output_array_pointers(13) = CCTK_PointerTo(Azzint)

  ! Set up interpolation coordinate arrays:
  !  n = 1

  ! define d's

  dcostheta =  1.0d0/ntheta
  dphi      =  2.0d0*PI/nphi

  ! if no symmetry then 
  if (Symmetry==NO_SYMM) dcostheta = 2.d0/ntheta


  do i=1,ntheta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,nphi
	n = j + (i-1)*nphi
        phiangle = (j - 0.5d0)*dphi
        xinterp(n) = x_cm + r_local*sintheta*cos(phiangle)
        yinterp(n) = y_cm + r_local*sintheta*sin(phiangle)
        zinterp(n) = z_cm + r_local*costheta
     end do
  end do


  interp_coords(1) = CCTK_PointerTo(xinterp)
  interp_coords(2) = CCTK_PointerTo(yinterp)
  interp_coords(3) = CCTK_PointerTo(zinterp)

  ! Perform interpolation:
  call CCTK_InterpGridArrays(ierr,cctkGH,N_dims,interp_handle, &
       param_table_handle,coord_system_handle, &
       N_interp_points,input_array_type_codes,interp_coords, &
       N_input_arrays,input_array_varindices, &
       N_output_arrays, output_array_type_codes, output_array_pointers)

  Psi6int = exp(6.D0*Psi6int)

  if(ierr.ge.0) then
     !|~~~~~> Integrate the dot product of Dphi with the surface normal:
     output_integralx= 0.D0
     output_integraly= 0.D0
     output_integralz= 0.D0

     do i=1,N_interp_points
        nxp= (xinterp(i)-x_cm)/r_local
        nyp= (yinterp(i)-y_cm)/r_local
        nzp= (zinterp(i)-z_cm)/r_local

        output_integralx = output_integralx + 0.125D0*Psi6int(i)*( &
          nxp*( (yinterp(i)-y_cm)*( Axzint(i)*gupxxint(i) + Ayzint(i)*gupxyint(i) + Azzint(i)*gupxzint(i) )  -&
                (zinterp(i)-z_cm)*( Axyint(i)*gupxxint(i) + Ayyint(i)*gupxyint(i) + Ayzint(i)*gupxzint(i) )) +&
          nyp*( (yinterp(i)-y_cm)*( Axzint(i)*gupxyint(i) + Ayzint(i)*gupyyint(i) + Azzint(i)*gupyzint(i) )  -&
                (zinterp(i)-z_cm)*( Axyint(i)*gupxyint(i) + Ayyint(i)*gupyyint(i) + Ayzint(i)*gupyzint(i) )) +&
          nzp*( (yinterp(i)-y_cm)*( Axzint(i)*gupxzint(i) + Ayzint(i)*gupyzint(i) + Azzint(i)*gupzzint(i) )  -&
                (zinterp(i)-z_cm)*( Axyint(i)*gupxzint(i) + Axyint(i)*gupyzint(i) + Ayzint(i)*gupzzint(i) )))

        output_integraly = output_integraly + 0.125D0*Psi6int(i)*( &
          nxp*( (zinterp(i)-z_cm)*( Axxint(i)*gupxxint(i) + Axyint(i)*gupxyint(i) + Axzint(i)*gupxzint(i) )  -&
                (xinterp(i)-x_cm)*( Axzint(i)*gupxxint(i) + Ayzint(i)*gupxyint(i) + Azzint(i)*gupxzint(i) )) +&
          nyp*( (zinterp(i)-z_cm)*( Axxint(i)*gupxyint(i) + Axyint(i)*gupyyint(i) + Axzint(i)*gupyzint(i) )  -&
                (xinterp(i)-x_cm)*( Axzint(i)*gupxyint(i) + Ayzint(i)*gupyyint(i) + Azzint(i)*gupyzint(i) )) +&
          nzp*( (zinterp(i)-z_cm)*( Axxint(i)*gupxzint(i) + Axyint(i)*gupyzint(i) + Axzint(i)*gupzzint(i) )  -&
                (xinterp(i)-x_cm)*( Axzint(i)*gupxzint(i) + Ayzint(i)*gupyzint(i) + Azzint(i)*gupzzint(i) )))

        output_integralz = output_integralz + 0.125D0*Psi6int(i)*( &
          nxp*( (xinterp(i)-x_cm)*( Axyint(i)*gupxxint(i) + Ayyint(i)*gupxyint(i) + Ayzint(i)*gupxzint(i) )  -&
                (yinterp(i)-y_cm)*( Axxint(i)*gupxxint(i) + Axyint(i)*gupxyint(i) + Axzint(i)*gupxzint(i) )) +&
          nyp*( (xinterp(i)-x_cm)*( Axyint(i)*gupxyint(i) + Ayyint(i)*gupyyint(i) + Ayzint(i)*gupyzint(i) )  -&
                (yinterp(i)-y_cm)*( Axxint(i)*gupxyint(i) + Axyint(i)*gupyyint(i) + Axzint(i)*gupyzint(i) )) +&
          nzp*( (xinterp(i)-x_cm)*( Axyint(i)*gupxzint(i) + Ayyint(i)*gupyzint(i) + Ayzint(i)*gupzzint(i) )  -&
                (yinterp(i)-y_cm)*( Axxint(i)*gupxzint(i) + Axyint(i)*gupyzint(i) + Axzint(i)*gupzzint(i) )))


     end do

     output_integralx = output_integralx * r_local * r_local * dphi * dcostheta / PI*sym_factor 
     output_integraly = output_integraly * r_local * r_local * dphi * dcostheta / PI*sym_factor 
     output_integralz = output_integralz * r_local * r_local * dphi * dcostheta / PI*sym_factor 
  else
     output_integralx = -9999.D0
     output_integraly = -9999.D0
     output_integralz = -9999.D0
  end if

end subroutine J_surf_integral_spin
