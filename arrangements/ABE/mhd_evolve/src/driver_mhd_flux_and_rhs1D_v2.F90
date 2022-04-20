!-----------------------------------------------------
! Compute MHD fluxes and rhs's in one direction: v2.0
!-----------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mhd_flux_and_rhs1D_v2(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS: 
  real*8                :: dX,dY,dZ
  integer, dimension(3) :: ext

  CCTK_REAL :: reduction_value
  integer :: dummy,handle,index,ierr
  integer :: i,j,k
  integer :: PPM_PLUS, PPM, CENO,MC,SPPM
  parameter(PPM_PLUS = 1, PPM = 2,CENO = 3,MC = 4, SPPM = 5)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

!  write(*,*) "Inside driver_mhd_flux_and_rhs1D_v2.F90, flux_direction = ",m
     
  call mhdfluxnew_hybrid(m,cctkGH,cctk_lsh,cctk_nghostzones, &
       enable_HARM_energyvariable, &
       X,dX,dY,dZ, &
       drho_b_m, drhoYe_m, dP_m, &
       dvx_m,dvy_m,dvz_m, &
       Pr, Pl,  &
       rho_br, rho_bl, Y_er, Y_el,&
       T_fluidle, T_fluidr, &
       vxr, vxl, vyr, vyl, vzr, vzl,  &
       Bxr, Bxl, Byr, Byl, Bzr, Bzl,  &
       v02r, v02l,  &
       gupxx_f, gupyy_f, gupzz_f,  &
       cmax, cmin,  &
       lapm1_f,  &
       shiftx_f, shifty_f, shiftz_f, & 
       gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f,  &
       phi_f, Symmetry,  &
       enable_OS_collapse, rho_b_atm,&
       neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab,gamma_th, & 
       use_central_scheme_instead_of_hll, compute_microphysics)

  
     

  ! now do advection
  !Advect rho and tau:
  call advect_rho_e_cpp(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry, &
       rho_star_rhs,tau_rhs,drho_b_m,dP_m, &
       dX,dY,dZ,X,Y,Z,drho_b_m_x, drho_b_m_xp1)


  !Advect S_i's
  call advect_z_cpp(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry, &
       X, mhd_st_y,vy,P,lapm1,phi, &
       shiftx,shifty,shiftz, &
       gxx,gxy,gxz,gyy,gyz,gzz,  &
       sbt,sbx,sby,sbz,  &
       mhd_st_x_rhs, mhd_st_y_rhs, mhd_st_z_rhs,  &
       dvx_m,dvy_m,dvz_m,  &
       dX,dY,dZ,dvx_m_x,dvx_m_xp1)

!  write(*,*) "Inside driver_mhd_flux_and_rhs1D_v2.F90, check1"  

  !Advect rho*Ye for M1 scheme: 
  if(rad_evolve_enable.eq.1) then
     call advect_rhoYe_cpp(m, cctkGH, cctk_lsh, cctk_nghostzones, Symmetry, &
       rhoYe_rhs, drhoYe_m, &
       dX,dY,dZ,X,Y,Z, drhoYe_m_x, drhoYe_m_xp1)
  end if


  !Add third-order accurate curvature terms to S_i_rhs
  
  
  call mhd_source_z_tau_cpp(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry, &
       enable_HARM_energyvariable, &
       dX, dY,dZ, Z, &
       mhd_st_x_rhs, mhd_st_y_rhs, mhd_st_z_rhs,  &
       tau_rhs, rho_star, & 
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
       E_rad,F_radx,F_rady,F_radz, Y_e, optd, eta_nue,&
       E_rad_nue,F_radx_nue,F_rady_nue,F_radz_nue, &
       E_rad_nux,F_radx_nux,F_rady_nux,F_radz_nux,&
       rad_fourforce_enable,rad_closure_scheme,rad_opacity_abs,rad_opacity_sct, rad_const, &
       chi_rad,chi_rad_nue,&
       T_fluid, enable_OS_collapse, compute_microphysics, microphysics_scheme,&
       T_fluid_cgs_atm, rad_fix)


!  write(*,*) "Inside driver_mhd_flux_and_rhs1D_v2.F90, check2" 

  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_pressure')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::micphys_conservatives') 
  call CartSymGN(dummy,cctkGH,'mhd_evolve::microphys_primitives')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::temperatures')
end subroutine mhd_flux_and_rhs1D_v2
