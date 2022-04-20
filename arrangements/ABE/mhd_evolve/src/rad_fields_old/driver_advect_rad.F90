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
  real*8                :: dX,dY,dZ,b2bt,pow_axi
  integer               :: index,ierr,handle,dummy
  CCTK_REAL             :: reduction_value
  integer               :: AXISYM,i,j,k
  parameter(AXISYM = 4)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  
  do k=1, cctk_lsh(3)
     do j=1, cctk_lsh(2)
          do i=1, cctk_lsh(1)
if (i==78.and.j==70.and.k==2) then
    write(*,*) " In driver_advect_rad.F90 checkpoint <<0>>  tau(i,j,k) is ", tau(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<0>>  tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<0>>  tau_rad(i,j,k) is ", tau_rad(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<0>>  tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<0>>  T_fluid(i,j,k) is ", T_fluid(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<0>>  rho_b(i,j,k), P(i,j,k) are ", rho_b(i,j,k), P(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<0>>  S_rad_x(i,j,k) = ",S_rad_x(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<0>>  S_rad_x_rhs(i,j,k) = ",S_rad_x_rhs(i,j,k)
end if
        end do
    end do
 end do



       call flux_rad_cpp (m, cctkGH, ext, X, Y,&
                   tau_rad_flux, &
                   S_radx_flux, S_rady_flux, S_radz_flux, &
                   E_radr, E_radl, &
                   F_radxr,F_radxl, &
                   F_radyr, F_radyl, &
                   F_radzr, F_radzl, &
                   Pr, Pl,&
                   rho_br, rho_bl,&
                   Bxr, Bxl, Byr, Byl, Bzr, Bzl,&
                   v02r, v02l,&
                   vxr,vxl,vyr,vyl,vzr,vzl, &
                   gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f, &
                   gupxx_f, gupyy_f, gupzz_f,gupxy_f, gupxz_f, gupyz_f, &
                   cmax,cmin, &
                   shiftx_f, shifty_f, shiftz_f, &
                   lapm1_f, phi_f, &
                   pow_axi,Symmetry, rad_closure_scheme,&
                   enable_OS_collapse,&
                   neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th)




  do k=1, cctk_lsh(3)
     do j=1, cctk_lsh(2)
          do i=1, cctk_lsh(1)
if (i==78.and.j==70.and.k==2) then
    write(*,*) " In driver_advect_rad.F90 checkpoint <<1>>  tau(i,j,k) is ", tau(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<1>>  tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<1>>  tau_rad(i,j,k) is ", tau_rad(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<1>>  tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<1>> T_fluid(i,j,k) is ", T_fluid(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<1>> rho_b(i,j,k), P(i,j,k) are ", rho_b(i,j,k), P(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<1>>  S_rad_x(i,j,k) = ",S_rad_x(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<1>>  S_rad_x_rhs(i,j,k) = ",S_rad_x_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<1>>  S_radx_flux(i,j,k) = ",S_radx_flux(i,j,k)
end if
      	end do
    end do
  end do
       

        call advect_rad_tau_ct_cpp(m, cctkGH,cctk_lsh, tau_rad_rhs, &
	     tau_rad_flux, Symmetry,dX,dY,dZ, tau_rad_flux_x, tau_rad_flux_xp1)



  
  call advect_Srad_cpp(m, cctkGH, cctk_lsh, cctk_nghostzones, Symmetry, &
             S_rad_x_rhs, S_rad_y_rhs, S_rad_z_rhs, &
	     S_radx_flux, S_rady_flux, S_radz_flux, &
  	     dX,dY,dZ, S_radx_flux_x, S_radx_flux_xp1)

    do k=1, cctk_lsh(3)
     do j=1, cctk_lsh(2)
          do i=1, cctk_lsh(1)
if (i==78.and.j==70.and.k==2) then
    write(*,*) " In driver_advect_rad.F90 checkpoint <<2>>  tau(i,j,k) is ", tau(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<2>>  tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<2>>  tau_rad(i,j,k) is ", tau_rad(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<2>>  tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<2>> T_fluid(i,j,k) is ", T_fluid(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<2>> rho_b(i,j,k), P(i,j,k) are ", rho_b(i,j,k), P(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<2>>  S_rad_x(i,j,k) = ",S_rad_x(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<2>>  S_rad_x_rhs(i,j,k) = ",S_rad_x_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<2>>  S_radx_flux(i,j,k) = ",S_radx_flux(i,j,k)
end if
end do
    end do
  end do
  
        write(*,*) "Before rad_source_cpp, rad_closure_scheme is", rad_closure_scheme
        write(*,*) "rad_opacity_abs is", rad_opacity_abs
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
       	     E_rad,F_radx,F_rady,F_radz,&
       	     rad_closure_scheme,rad_opacity_abs,rad_opacity_sct, rad_const, &
	     T_fluid)


     if(excision_enable == 1) then
        call remove_interior(ext,X,Y,Z,tau_rad_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,S_rad_x_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,S_rad_y_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,S_rad_z_rhs,excision_zone_gf,Symmetry)      
     end if



 call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_conservatives')
 call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives')

  do k=1, cctk_lsh(3)
     do j=1, cctk_lsh(2)
          do i=1, cctk_lsh(1)
if (i==78.and.j==70.and.k==2) then
    write(*,*) " In driver_advect_rad.F90 checkpoint <<3>>  tau(i,j,k) is ", tau(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<3>>  tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<3>>  tau_rad(i,j,k) is ", tau_rad(i,j,k)
    write(*,*) " In driver_advect_rad.F90 checkpoint <<3>>  tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<3>> T_fluid(i,j,k) is ", T_fluid(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<3>> rho_b(i,j,k), P(i,j,k) are ", rho_b(i,j,k), P(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<3>>  S_rad_x(i,j,k) = ",S_rad_x(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<3>>  S_rad_x_rhs(i,j,k) = ",S_rad_x_rhs(i,j,k)
    write(*,*) " In driver_advect_rad.F90,checkpoint <<3>>  S_radx_flux(i,j,k) = ",S_radx_flux(i,j,k)
end if
        end do
    end do
  end do
 



end subroutine driver_advect_rad
