!-----------------------------------------------------------------------------
! Compute Komar mass via surface integral
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine Komar_surf_integralv2(cctkGH,output_integral)
  implicit none
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  ! Variables in the function call
  CCTK_POINTER                             :: cctkGH
  real*8                                   :: output_integral
  ! Variables needed for interpolation
  CCTK_POINTER, dimension(3)               :: interp_coords
  character(60)                            :: options_string
  integer                                  :: interpolation_order,nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle,N_interp_points
  integer                                  :: N_input_arrays,N_output_arrays
  real*8, dimension(N_theta*N_phi)         :: xinterp,yinterp,zinterp
  integer,dimension(26)                    :: input_array_type_codes,input_array_varindices,output_array_type_codes
  CCTK_POINTER,dimension(26)               :: output_array_pointers
  ! Output arrays, dummy indices, parameters
  real*8,dimension(N_theta*N_phi)          :: Psi4int,gupxxint,gupxyint,gupxzint,gupyyint,gupyzint,gupzzint
  real*8,dimension(N_theta*N_phi)          :: gxxint,gxyint,gxzint,gyyint,gyzint,gzzint
  real*8,dimension(N_theta*N_phi)          :: Axxint,Axyint,Axzint,Ayyint,Ayzint,Azzint,trKint
  real*8,dimension(N_theta*N_phi)          :: lapsexint,lapseyint,lapsezint
  real*8,dimension(N_theta*N_phi)          :: shiftxint,shiftyint,shiftzint
  real*8                                   :: phiangle,costheta,sintheta,PI,nxp,nyp,nzp
  integer                                  :: i,j,n

  PI = 3.14159265358979323846D0

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
  N_input_arrays = 26
  N_output_arrays = 26

  call CCTK_VarIndex (input_array_varindices(1), "bssn::phi")
  call CCTK_VarIndex (input_array_varindices(2), "bssn::gupxx")
  call CCTK_VarIndex (input_array_varindices(3), "bssn::gupxy")
  call CCTK_VarIndex (input_array_varindices(4), "bssn::gupxz")
  call CCTK_VarIndex (input_array_varindices(5), "bssn::gupyy")
  call CCTK_VarIndex (input_array_varindices(6), "bssn::gupyz")
  call CCTK_VarIndex (input_array_varindices(7), "bssn::gupzz")

  call CCTK_VarIndex (input_array_varindices(8), "bssn::gxx")
  call CCTK_VarIndex (input_array_varindices(9), "bssn::gxy")
  call CCTK_VarIndex (input_array_varindices(10), "bssn::gxz")
  call CCTK_VarIndex (input_array_varindices(11), "bssn::gyy")
  call CCTK_VarIndex (input_array_varindices(12), "bssn::gyz")
  call CCTK_VarIndex (input_array_varindices(13), "bssn::gzz")

  call CCTK_VarIndex (input_array_varindices(14), "bssn::Axx")
  call CCTK_VarIndex (input_array_varindices(15), "bssn::Axy")
  call CCTK_VarIndex (input_array_varindices(16), "bssn::Axz")
  call CCTK_VarIndex (input_array_varindices(17), "bssn::Ayy")
  call CCTK_VarIndex (input_array_varindices(18), "bssn::Ayz")
  call CCTK_VarIndex (input_array_varindices(19), "bssn::Azz")

  call CCTK_VarIndex (input_array_varindices(20), "bssn::trK")

  call CCTK_VarIndex (input_array_varindices(21), "lapse::lapsex")
  call CCTK_VarIndex (input_array_varindices(22), "lapse::lapsey")
  call CCTK_VarIndex (input_array_varindices(23), "lapse::lapsez")

  call CCTK_VarIndex (input_array_varindices(24), "shift::shiftx")
  call CCTK_VarIndex (input_array_varindices(25), "shift::shifty")
  call CCTK_VarIndex (input_array_varindices(26), "shift::shiftz")

! FIXME?!
!  call CCTK_SyncGroup(i,cctkGH,'lapse::lapse_derivatives')

  output_array_pointers(1) = CCTK_PointerTo(Psi4int)
  output_array_pointers(2) = CCTK_PointerTo(gupxxint)
  output_array_pointers(3) = CCTK_PointerTo(gupxyint)
  output_array_pointers(4) = CCTK_PointerTo(gupxzint)
  output_array_pointers(5) = CCTK_PointerTo(gupyyint)
  output_array_pointers(6) = CCTK_PointerTo(gupyzint)
  output_array_pointers(7) = CCTK_PointerTo(gupzzint)
  output_array_pointers(8) = CCTK_PointerTo(gxxint)
  output_array_pointers(9) = CCTK_PointerTo(gxyint)
  output_array_pointers(10) = CCTK_PointerTo(gxzint)
  output_array_pointers(11) = CCTK_PointerTo(gyyint)
  output_array_pointers(12) = CCTK_PointerTo(gyzint)
  output_array_pointers(13) = CCTK_PointerTo(gzzint)
  output_array_pointers(14) = CCTK_PointerTo(Axxint)
  output_array_pointers(15) = CCTK_PointerTo(Axyint)
  output_array_pointers(16) = CCTK_PointerTo(Axzint)
  output_array_pointers(17) = CCTK_PointerTo(Ayyint)
  output_array_pointers(18) = CCTK_PointerTo(Ayzint)
  output_array_pointers(19) = CCTK_PointerTo(Azzint)
  output_array_pointers(20) = CCTK_PointerTo(trKint)
  output_array_pointers(21) = CCTK_PointerTo(lapsexint)
  output_array_pointers(22) = CCTK_PointerTo(lapseyint)
  output_array_pointers(23) = CCTK_PointerTo(lapsezint)
  output_array_pointers(24) = CCTK_PointerTo(shiftxint)
  output_array_pointers(25) = CCTK_PointerTo(shiftyint)
  output_array_pointers(26) = CCTK_PointerTo(shiftzint)

  ! Set up interpolation coordinate arrays:
  n = 1
  do i=1,N_theta
     costheta = 1.D0 - (i - 0.5D0)*dcostheta
     sintheta = sqrt(1.D0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5D0)*dphi
        if(N_phi==1) phiangle = 0.D0
        xinterp(n) = surf_radius*sintheta*cos(phiangle)
        yinterp(n) = surf_radius*sintheta*sin(phiangle)
        zinterp(n) = surf_radius*costheta
        n = n + 1
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

  Psi4int = exp(4.D0*Psi4int)

  !Convert Aijint to Kijint!
  Axxint = (Axxint + (1.D0/3.D0)*trKint*gxxint)*Psi4int
  Axyint = (Axyint + (1.D0/3.D0)*trKint*gxyint)*Psi4int
  Axzint = (Axzint + (1.D0/3.D0)*trKint*gxzint)*Psi4int
  Ayyint = (Ayyint + (1.D0/3.D0)*trKint*gyyint)*Psi4int
  Ayzint = (Ayzint + (1.D0/3.D0)*trKint*gyzint)*Psi4int
  Azzint = (Azzint + (1.D0/3.D0)*trKint*gzzint)*Psi4int

  if(ierr.ge.0) then
     !|~~~~~> Integrate the dot product of Dphi with the surface normal:
     output_integral = 0.D0
     do i=1,ntot
        output_integral = output_integral + sqrt(Psi4int(i))*( &
             (lapsexint(i) - Axxint(i)*shiftxint(i) - Axyint(i)*shiftyint(i) - Axzint(i)*shiftzint(i)) &
             *( xinterp(i)*gupxxint(i) + yinterp(i)*gupxyint(i) + zinterp(i)*gupxzint(i) ) + &
             (lapseyint(i) - Axyint(i)*shiftxint(i) - Ayyint(i)*shiftyint(i) - Ayzint(i)*shiftzint(i)) &
             *( xinterp(i)*gupxyint(i) + yinterp(i)*gupyyint(i) + zinterp(i)*gupyzint(i) ) + &
             (lapsezint(i) - Axzint(i)*shiftxint(i) - Ayzint(i)*shiftyint(i) - Azzint(i)*shiftzint(i)) &
             *( xinterp(i)*gupxzint(i) + yinterp(i)*gupyzint(i) + zinterp(i)*gupzzint(i) ))/surf_radius
     end do
     !  Multiply by -1/2\pi * surface element.
     output_integral = output_integral * surf_radius*surf_radius * dphi * dcostheta * (1.D0/(4.D0*PI))*sym_factor
  else
     output_integral = -9999.D0
  end if

end subroutine Komar_surf_integralv2
