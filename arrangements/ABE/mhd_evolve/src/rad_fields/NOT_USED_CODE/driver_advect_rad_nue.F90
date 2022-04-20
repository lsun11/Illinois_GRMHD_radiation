!-------------------------------------------------------
!    :: Driver routine for MHD timestepping, v2.0 ::
! (i.e., computing RHS's of all conservative variables)
!-------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_advect_rad_nue(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dX,dY,dZ,b2bt,pow_axi
  integer               :: index,ierr,handle,dummy
  CCTK_REAL             :: reduction_value
  integer               :: i,j,k

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)


        do k=1,ext(3)
           do j=1,ext(2)
              do i=1,ext(1)
                 if (isnan(tau_rad(i,j,k)) .or. isnan(S_rad_x(i,j,k))) then
                    write(*,*) "inside driver_advect_rad.F90[[1]], tau_rad and S_rad_x=", tau_rad(i,j,k), S_rad_x(i,j,k)
                    write(*,*) "i,j,k =", i,j,k
                 end if
              end do
           end do
        end do


!  write(*,*) "Inside driver_advect_rad.F90, flux_direction = ", m      

  call flux_rad_nue_cpp (m, cctkGH, ext, X, Y, Z,&
       tau_rad_nue_flux,&
       S_radx_nue_flux, S_rady_nue_flux, S_radz_nue_flux, &
       E_rad_nuer, E_rad_nuel, &
       F_radx_nuer, F_radx_nuele, &
       F_rady_nuer, F_rady_nuele, &
       F_radz_nuer, F_radz_nuele, &
       FaFar_nue, FaFal_nue,&
       Pr, Pl,&
       rho_br, rho_bl,&
       Bxr, Bxl, Byr, Byl, Bzr, Bzl,&
       v02_rad_nuer, v02_rad_nuel,&
       vxr,vxl,vyr,vyl,vzr,vzl, &
       gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f, &
       gupxx_f, gupyy_f, gupzz_f,gupxy_f, gupxz_f, gupyz_f, &
       cmax_rad_nue,cmin_rad_nue, &
       shiftx_f, shifty_f, shiftz_f, &
       lapm1_f, phi_f, &
       pow_axi,Symmetry, rad_closure_scheme,&
       enable_OS_collapse,&
       neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th,&
       Erad_atm_cut)

  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_fijs')


  call advect_rad_tau_ct_nue_cpp(m, cctkGH,cctk_lsh, tau_rad_nue_rhs, &
       tau_rad_nue_flux, Symmetry,dX,dY,dZ)


  call advect_Srad_nue_cpp(m, cctkGH, cctk_lsh, cctk_nghostzones, Symmetry, &
       S_rad_x_nue_rhs, S_rad_y_nue_rhs, S_rad_z_nue_rhs, &
       S_radx_nue_flux, S_rady_nue_flux, S_radz_nue_flux, &
       dX,dY,dZ)


   call rad_source_nue_cpp(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry, &
       enable_HARM_energyvariable, &
       dX, dY,dZ,&
       S_rad_x_nue_rhs, S_rad_y_nue_rhs, S_rad_z_nue_rhs,  &
       tau_rad_nue_rhs, rho_star, &
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
       E_rad_nue,F_radx_nue,F_rady_nue,F_radz_nue,Y_e, optd,&
       rad_closure_scheme,rad_opacity_abs,rad_opacity_sct, rad_const, &
       ka_gf_nue, ks_gf_nue, emission_gf_nue, chi_rad_nue,&
       T_fluid, Erad_atm_cut, enable_OS_collapse, compute_microphysics, rad_fourforce_enable, T_fluid_cgs_atm)


        do k=1,ext(3)
           do j=1,ext(2)
              do i=1,ext(1)
                 if (isnan(tau_rad(i,j,k)) .or. isnan(S_rad_x(i,j,k))) then
                    write(*,*) "inside driver_advect_rad.F90[[2]], tau_rad and S_rad_x=", tau_rad(i,j,k), S_rad_x(i,j,k)
                    write(*,*) "i,j,k =", i,j,k
                 end if
              end do
           end do
        end do


  if(excision_enable == 1) then
     call remove_interior(ext,X,Y,Z,tau_rad_rhs,excision_zone_gf,Symmetry)
     call remove_interior(ext,X,Y,Z,S_rad_x_rhs,excision_zone_gf,Symmetry)
     call remove_interior(ext,X,Y,Z,S_rad_y_rhs,excision_zone_gf,Symmetry)
     call remove_interior(ext,X,Y,Z,S_rad_z_rhs,excision_zone_gf,Symmetry)      
  end if
  
  
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')  
  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives') 
  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_pressure')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::temperatures')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::micphys_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::microphys_primitives')
  
  
end subroutine driver_advect_rad_nue
