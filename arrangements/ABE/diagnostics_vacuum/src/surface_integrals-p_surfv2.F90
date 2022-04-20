!-----------------------------------------------------------------------------
! Compute 1/8\pi \int \psi^6 (x A^m_y - yA^m_x) dS_m.
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine P_surf_integralv2(cctkGH,output_integralx,output_integraly,output_integralz)
  implicit none
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  ! Variables in the function call
  CCTK_POINTER                             :: cctkGH
  real*8                                   :: output_integralx,output_integraly,output_integralz
  ! Variables needed for interpolation
  CCTK_POINTER, dimension(3)               :: interp_coords
  character(60)                            :: options_string
  integer                                  :: interpolation_order,nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle,N_interp_points
  integer                                  :: N_input_arrays,N_output_arrays
  real*8, dimension(N_theta*N_phi)         :: xinterp,yinterp,zinterp
  integer,dimension(14)                    :: input_array_type_codes,input_array_varindices,output_array_type_codes
  CCTK_POINTER,dimension(14)               :: output_array_pointers
  ! Output arrays, dummy indices, parameters
  real*8,dimension(N_theta*N_phi)          :: Psi6int,gupxxint,gupxyint,gupxzint,gupyyint,gupyzint,gupzzint
  real*8,dimension(N_theta*N_phi)          :: Axxint,Axyint,Axzint,Ayyint,Ayzint,Azzint,trKint
  real*8                                   :: phiangle,costheta,sintheta,PI,nxp,nyp,nzp
  real*8                                   :: AupxxL2,AupxyL2,AupxzL2,AupyxL2,AupyyL2,AupyzL2,AupzxL2,AupzyL2,AupzzL2
  real*8                                   :: AupxxL,AupxyL,AupxzL,AupyxL,AupyyL,AupyzL,AupzxL,AupzyL,AupzzL
  real*8                                   :: xhat,yhat,zhat
  real*8                                   :: out_x,out_y,out_z
  real*8                                   :: juxx,juyy,juzz,juxy,juxz,juyz,rr,drr
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
  N_input_arrays = 14
  N_output_arrays = 14

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

  call CCTK_VarIndex (input_array_varindices(14), "bssn::trK")

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
  output_array_pointers(14) = CCTK_PointerTo(trKint)

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

  Psi6int = exp(6.D0*Psi6int)

  if(ierr.ge.0) then
     !|~~~~~> Integrate the dot product of Dphi with the surface normal
     output_integralx = 0.D0
     output_integraly = 0.D0
     output_integralz = 0.D0
     rr=1.0/rbr
     drr=1.0/drbdr
     
     do i=1,ntot
        AupxxL=Axxint(i)*gupxxint(i)+Axyint(i)*gupxyint(i)+Axzint(i)*gupxzint(i)
        AupxyL=Axyint(i)*gupxxint(i)+Ayyint(i)*gupxyint(i)+Ayzint(i)*gupxzint(i)
        AupxzL=Axzint(i)*gupxxint(i)+Ayzint(i)*gupxyint(i)+Azzint(i)*gupxzint(i)
        AupyxL=Axxint(i)*gupxyint(i)+Axyint(i)*gupyyint(i)+Axzint(i)*gupyzint(i)
        AupyyL=Axyint(i)*gupxyint(i)+Ayyint(i)*gupyyint(i)+Ayzint(i)*gupyzint(i)
        AupyzL=Axzint(i)*gupxyint(i)+Ayzint(i)*gupyyint(i)+Azzint(i)*gupyzint(i)
        AupzxL=Axxint(i)*gupxzint(i)+Axyint(i)*gupyzint(i)+Axzint(i)*gupzzint(i)
        AupzyL=Axyint(i)*gupxzint(i)+Ayyint(i)*gupyzint(i)+Ayzint(i)*gupzzint(i)
        AupzzL=Axzint(i)*gupxzint(i)+Ayzint(i)*gupyzint(i)+Azzint(i)*gupzzint(i)
        
        xhat=xinterp(i)/surf_radius
        yhat=yinterp(i)/surf_radius
        zhat=zinterp(i)/surf_radius
        
        juxx=xhat**2*(drr-rr)+rr
        juxy=xhat*yhat*(drr-rr)
        juxz=xhat*zhat*(drr-rr)
        juyy=yhat**2*(drr-rr)+rr
        juyz=yhat*zhat*(drr-rr)
        juzz=zhat**2*(drr-rr)+rr
        
        AupxxL2=AupxxL*juxx+AupxyL*juxy+AupxzL*juxz
        AupxyL2=AupxxL*juxy+AupxyL*juyy+AupxzL*juyz
        AupxzL2=AupxxL*juxz+AupxyL*juyz+AupxzL*juzz
        AupyxL2=AupyxL*juxx+AupyyL*juxy+AupyzL*juxz
        AupyyL2=AupyxL*juxy+AupyyL*juyy+AupyzL*juyz
        AupyzL2=AupyxL*juxz+AupyyL*juyz+AupyzL*juzz
        AupzxL2=AupzxL*juxx+AupzyL*juxy+AupzzL*juxz
        AupzyL2=AupzxL*juxy+AupzyL*juyy+AupzzL*juyz
        AupzzL2=AupzxL*juxz+AupzyL*juyz+AupzzL*juzz
        
        out_x=xhat*AupxxL2+yhat*AupyxL2+zhat*AupzxL2-2./3.*xhat*trkint(i)*drr
        out_y=xhat*AupxyL2+yhat*AupyyL2+zhat*AupzyL2-2./3.*yhat*trkint(i)*drr
        out_z=xhat*AupxzL2+yhat*AupyzL2+zhat*AupzzL2-2./3.*zhat*trkint(i)*drr
        
        output_integralx=output_integralx+0.125D0*Psi6int(i)*out_x
        output_integraly=output_integraly+0.125D0*Psi6int(i)*out_y
        output_integralz=output_integralz+0.125D0*Psi6int(i)*out_z
     end do

     output_integralx = output_integralx * surf_radius * surf_radius * dphi * dcostheta / PI *sym_factor
     output_integraly = output_integraly * surf_radius * surf_radius * dphi * dcostheta / PI *sym_factor
     output_integralz = output_integralz * surf_radius * surf_radius * dphi * dcostheta / PI *sym_factor
  else
     output_integralx = -9999.D0
     output_integraly = -9999.D0
     output_integralz = -9999.D0
  end if

end subroutine P_surf_integralv2
