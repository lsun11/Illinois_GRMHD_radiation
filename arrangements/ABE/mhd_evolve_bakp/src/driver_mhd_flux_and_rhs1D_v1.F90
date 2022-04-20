!-----------------------------------------------------
! Compute MHD fluxes and rhs's in one direction: v1.0
!-----------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mhd_flux_and_rhs1D_v1(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS: 
  real*8                :: dX,dY,dZ
  integer, dimension(3) :: ext

  CCTK_REAL :: reduction_value
  integer :: dummy,handle,index,ierr
  integer :: PPM_PLUS, PPM, CENO,MC,SPPM
  parameter(PPM_PLUS = 1, PPM = 2,CENO = 3,MC = 4, SPPM = 5)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(m==1 .and. cctk_iteration.eq.1 .and. iter_count.eq.1) then
     write(*,*) "IF THE CODE CRASHES HERE, YOU PROBABLY COMPILED WITH A NEW INTEL COMPILER,"
     write(*,*) "WHICH IS INCOMPATIBLE WITH use_new_code==0 !"
  end if
  
  if(cctk_nghostzones(2) .gt. 1) then
     write(*,*) "ERROR: use_new_code = 0 (i.e., the OLD mhd_evolve) does NOT support arbitrary number of ghostzones!"
     stop
  end if

  call mhdflux_hybrid(ext,X,drho_b_m,dP_m,dvx_m,dvy_m, &
       dvz_m,Pr,Pl,rho_br,rho_bl, &
       vxr,vxl,vyr,vyl,vzr,vzl, &
       Bxr,Bxl,Byr,Byl,Bzr,Bzl, &
       v02r,v02l,gupxx_f, &
       gupyy_f,gupzz_f,cmax,cmin,m, &
       lapm1_f, shiftx_f, shifty_f, shiftz_f, &
       gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f, &
       phi_f, Symmetry,neos, ergo_star,ergo_sigma,&
       rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th,enable_HARM_energyvariable)

  if(excision_enable == 1) then
     call remove_interior2(ext,X,Y,Z,drho_b_m,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,dP_m,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,dvx_m,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,dvy_m,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,dvz_m,excision_zone_gf,Symmetry)
  end if

  ! Do higher-order flux derivatives?
  !  if(Integration_Order > 2) then
  !     fs_to_fhats(m)
  !  end if

  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_nablas')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_nablas')

  ! now do advection
  call advect_rho_e(ext,rho_star_rhs,tau_rhs, &
       drho_b_m,dP_m,m, &
       X,Y,Z,Symmetry)

  call advect_z(ext,X,mhd_st_y,vy,P,lapm1,phi, & 
       shiftx,shifty,shiftz, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       sbt,sbx,sby,sbz, &
       mhd_st_x_rhs,mhd_st_y_rhs,mhd_st_z_rhs, &
       dvx_m,dvy_m,dvz_m,m,dX,dY,dZ,Symmetry)

  call mhd_source_z_tau(ext,dX,dY,dZ,mhd_st_x_rhs,mhd_st_y_rhs,mhd_st_z_rhs, &
       tau_rhs, &
       rho_star,h,u0,vx,vy,vz,P, &
       sbt,sbx,sby,sbz, lapm1, &
       shiftx,shifty,shiftz,phi, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Symmetry, &
       lapm1_f,shiftx_f,shifty_f,shiftz_f, &
       gxx_f,gxy_f,gxz_f,gyy_f,gyz_f,gzz_f, &
       phi_f,m,Z,enable_HARM_energyvariable)

  if(excision_enable == 1) then
     call scalar_excision_bc(ext,X,Y,Z, &
          rho_star_rhs,Symmetry,excision_zone_gf)
     call scalar_excision_bc(ext,X,Y,Z, &
          tau_rhs,Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z, &
          mhd_st_x_rhs,mhd_st_y_rhs,mhd_st_z_rhs, &
          Symmetry,excision_zone_gf)
  end if

  !  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_temp')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')

end subroutine mhd_flux_and_rhs1D_v1
