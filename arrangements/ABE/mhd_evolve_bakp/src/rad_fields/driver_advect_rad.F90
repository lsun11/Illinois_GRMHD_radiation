!-------------------------------------------------------
!    :: Driver routine for MHD timestepping, v2.0 ::
! (i.e., computing RHS's of all conservative variables)
!-------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_advect_rad(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dX,dY,dZ,b2bt
  integer               :: index,ierr,handle,dummy
  CCTK_REAL             :: reduction_value
  integer               :: AXISYM,i,j,k
  parameter(AXISYM = 4)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(rad_evolve_enable==1) then


        call advect_rad_tau_ct_cpp(m, cctkGH,ext,X,Y,Z, tau_rad_rhs, &
	     tau_rad_flux, Symmetry)


	call advect_Srad_cpp(m, cctkGH, cctk_lsh, cctk_nghostzones, Symmetry, &
             S_rad_x_rhs, S_rad_y_rhs, S_rad_z_rhs, &
	     S_radx_flux, S_rady_flux, S_radz_flux, &
  	     dX,dY,dZ)


	 

        call rad_source_cpp(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry, &
       	     enable_HARM_energyvariable, &
       	     dX, dY,dZ, Z, &
       	     S_rad_x_rhs, S_rad_y_rhs, S_rad_z_rhs,  &
	     tau_rad_rhs, rho_star, &
       	     P, h, u0,  &
       	     vx, vy, vz,  &
       	     sbt, sbx, sby, sbz, &
       	     lapm1, shiftx, shifty, shiftz,  &
       	     phi,  &
       	     gxx, gxy, gxz, gyy, gyz, gzz,  &
       	     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,  &
       	     lapm1_f, shiftx_f, shifty_f, shiftz_f,  &
       	     phi_f,  &
       	     gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f, &
       	     temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8, &
       	     E_rad,F_rad0,F_radx,F_rady,F_radz,P_rad, &
       	     rad_closure_scheme,rad_opacity_abs,rad_opacity_sct, rad_const, &
	     T_fluid)

	     

     if(excision_enable == 1) then
        call remove_interior(ext,X,Y,Z,tau_rad_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,S_rad_x_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,S_rad_y_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,S_rad_z_rhs,excision_zone_gf,Symmetry)      
     end if

 end if 

end subroutine driver_advect_rad
