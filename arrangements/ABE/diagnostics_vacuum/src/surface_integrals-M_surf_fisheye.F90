!-----------------------------------------------------------------------------
! Calculate ADM mass via surface integral
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine M_surf_integral(cctkGH,output_integral)
  implicit none
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

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
  real*8                                   :: rphys,fac1,fac2,f4,df4,df4x,df4y,df4z,nx,ny,nz
  real*8                                   :: jxx,jxy,jxz,jyy,jyz,jzz,jxxx,jyyy,jzzz,jxyz
  real*8                                   :: jxxy,jxyy,jxxz,jxzz,jyyz,jyzz,gamx,gamy,gamz
  integer                                  :: i,j,n

  PI = 3.14159265358979323846D0

  N_dims = 3
  interpolation_order = 2
  N_interp_points = N_theta*N_phi

  interp_handle = -1 
!  call CCTK_InterpHandle (interp_handle, "uniform cartesian")
  call CCTK_InterpHandle (interp_handle, "Lagrange polynomial interpolation")
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

  Psiint = exp(Psiint)

  !We're looking at an r=constant surface, and fisheye transfs. are spherically symmetric, so rphys,dR,and ddR are fixed.
  rphys = surf_radius*rbr
  fac1 = (rbr)**(1.D0/3.D0) * drbdr**(1.D0/6.D0)
  fac2 = (drbdr/rphys-1.D0/surf_radius)/3.D0 + ddrbddr/(drbdr*6.D0)

  ! we need to calculate \tilde{Gamma}^i from \bar{\tilde{Gamma}}^i, via jacobian and derivatives
  ! f4 is 'F' in my notes, and df4*nx gives you dF/dx
  f4=fac1**(-4)
  df4=f4/3.0d0*(-4.0d0/rphys*(drbdr-rbr)-2.0d0*ddrbddr/drbdr)

!  n=ntot/50

  if(ierr.ge.0) then
     output_integral = 0.D0
     !|~~~~~> Integrate the dot product of Dphi with the surface normal
     do i=1,ntot
        !Convert phix, phiy, and phiz to physical coordinates:
        nx= xinterp(i)/surf_radius        
        ny= yinterp(i)/surf_radius
        nz= zinterp(i)/surf_radius

        df4x=df4*nx
        df4y=df4*ny
        df4z=df4*nz

        phixint(i) = phixint(i) - nx*fac2
        phiyint(i) = phiyint(i) - ny*fac2
        phizint(i) = phizint(i) - nz*fac2
        !

        Jxx=rbr+nx*nx*(drbdr-rbr)
        Jxy=nx*ny*(drbdr-rbr)
        Jxz=nx*nz*(drbdr-rbr)
        Jyy=rbr+ny*ny*(drbdr-rbr)
        Jyz=ny*nz*(drbdr-rbr)
        Jzz=rbr+nz*nz*(drbdr-rbr)

        jxxx=nx**3*ddrbddr+(drbdr-rbr)/surf_radius*(3.0d0*nx-3.0d0*nx**3)
        jyyy=ny**3*ddrbddr+(drbdr-rbr)/surf_radius*(3.0d0*ny-3.0d0*ny**3)
        jzzz=nz**3*ddrbddr+(drbdr-rbr)/surf_radius*(3.0d0*nz-3.0d0*nz**3)

        jxxy=nx**2*ny*ddrbddr+(drbdr-rbr)/surf_radius*(ny-3.0d0*nx**2*ny)
        jxyy=nx*ny**2*ddrbddr+(drbdr-rbr)/surf_radius*(nx-3.0d0*nx*ny**2)
        jxxz=nx**2*nz*ddrbddr+(drbdr-rbr)/surf_radius*(nz-3.0d0*nx**2*nz)
        jxzz=nx*nz**2*ddrbddr+(drbdr-rbr)/surf_radius*(nx-3.0d0*nx*nz**2)
        jyyz=ny**2*nz*ddrbddr+(drbdr-rbr)/surf_radius*(nz-3.0d0*ny**2*nz)
        jyzz=ny*nz**2*ddrbddr+(drbdr-rbr)/surf_radius*(ny-3.0d0*ny*nz**2)

        jxyz=nx*ny*nz*(ddrbddr-3.0d0*(drbdr-rbr)/surf_radius)

        Gamx=f4*(jxx*gamxint(i)+jxy*gamyint(i)+jxz*gamzint(i))+ &
             gupxxint(i)*(0.5d0*df4x*jxx-f4*jxxx)+ &
             gupyyint(i)*(0.5d0*df4y*jxy-f4*jxyy)+ &
             gupzzint(i)*(0.5d0*df4z*jxz-f4*jxzz)+ &
             gupxyint(i)*(0.5d0*(df4x*jxy+df4y*jxx)-2.0d0*f4*jxxy)+ &
             gupxzint(i)*(0.5d0*(df4x*jxz+df4z*jxx)-2.0d0*f4*jxxz)+ &
             gupyzint(i)*(0.5d0*(df4y*jxz+df4z*jxy)-2.0d0*f4*jxyz)
        Gamy=f4*(jxy*gamxint(i)+jyy*gamyint(i)+jyz*gamzint(i))+ &
             gupxxint(i)*(0.5d0*df4x*jxy-f4*jxxy)+ &
             gupyyint(i)*(0.5d0*df4y*jyy-f4*jyyy)+ &
             gupzzint(i)*(0.5d0*df4z*jyz-f4*jyzz)+ &
             gupxyint(i)*(0.5d0*(df4x*jyy+df4y*jxy)-2.0d0*f4*jxyy)+ &
             gupxzint(i)*(0.5d0*(df4x*jyz+df4z*jxy)-2.0d0*f4*jxyz)+ &
             gupyzint(i)*(0.5d0*(df4y*jyz+df4z*jyy)-2.0d0*f4*jyyz) 
        Gamz=f4*(jxz*gamxint(i)+jyz*gamyint(i)+jzz*gamzint(i))+ &
             gupxxint(i)*(0.5d0*df4x*jxz-f4*jxxz)+ &
             gupyyint(i)*(0.5d0*df4y*jyz-f4*jyyz)+ &
             gupzzint(i)*(0.5d0*df4z*jzz-f4*jzzz)+ &
             gupxyint(i)*(0.5d0*(df4x*jyz+df4y*jxz)-2.0d0*f4*jxyz)+ &
             gupxzint(i)*(0.5d0*(df4x*jzz+df4z*jxz)-2.0d0*f4*jxzz)+ &
             gupyzint(i)*(0.5d0*(df4y*jzz+df4z*jyz)-2.0d0*f4*jyzz) 

        output_integral = output_integral + Psiint(i) * fac1 * &
             (nx*(gupxxint(i)*phixint(i)+gupxyint(i)*phiyint(i)+gupxzint(i)*phizint(i))+ &
             ny*(gupxyint(i)*phixint(i)+gupyyint(i)*phiyint(i)+gupyzint(i)*phizint(i))+ & 
             nz*(gupxzint(i)*phixint(i)+gupyzint(i)*phiyint(i)+gupzzint(i)*phizint(i)))- &
             rbr**2/8.0d0*(nx*gamx+ny*Gamy+nz*Gamz)

     end do

     !|  Multiply by -1/2\pi * surface element.
     output_integral = output_integral * -0.5D0/PI * surf_radius * surf_radius * dphi * dcostheta*sym_factor
  else
     output_integral = -9999.D0
  end if

end subroutine M_surf_integral
