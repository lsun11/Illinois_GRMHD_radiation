!-----------------------------------------------------------------------------
!
! $Id: lin_wave.f90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!
!-----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bns_copy_to_prev_timelevel(CCTK_ARGUMENTS)
  implicit none
  integer :: dummy,index
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  lapm1_p = lapm1
  lapset_p = lapset
  shiftx_p = shiftx
  shifty_p = shifty
  shiftz_p = shiftz

  phi_p = phi 

  gxx_p = gxx
  gyy_p = gyy
  gzz_p = gzz
  gxy_p = gxy
  gxz_p = gxz
  gyz_p = gyz

  Gammax_p = Gammax
  Gammay_p = Gammay
  Gammaz_p = Gammaz
  !======================================
  ! Set other matter sources to prev timelevel values
  !======================================
  Sx_p = Sx
  Sy_p = Sy
  Sz_p = Sz
  !======================================
  ! Set K_ij to prev timelevel value
  !======================================
  trK_p = trK
  Axx_p = Axx
  Axy_p = Axy
  Axz_p = Axz
  Ayy_p = Ayy
  Ayz_p = Ayz
  Azz_p = Azz

  !======================================
  ! Set rhs's to Zero!
  !======================================
  shiftx_rhs = 0.D0
  shifty_rhs = 0.D0
  shiftz_rhs = 0.D0

  !===========================================
  ! Finally, copy remaining matter variables to _p level
  !===========================================
  rho_p = rho
  S_p = S
  Sxx_p = Sxx
  Sxy_p = Sxy
  Sxz_p = Sxz
  Syy_p = Syy
  Syz_p = Syz
  Szz_p = Szz

  rho_star_p = rho_star
  w_p = w
  h_p = h
  tau_p = tau
  st_x_p = st_x
  st_y_p = st_y
  st_z_p = st_z
  mhd_st_x_p = mhd_st_x
  mhd_st_y_p = mhd_st_y
  mhd_st_z_p = mhd_st_z

  Bx_p = Bx
  By_p = By
  Bz_p = Bz


  !-----------------------------------------
  !Set aux. MHD gridfunctions to zero, so that valgrind is happy.
  !Should probably comment these out if you're debugging the code
  gupxx_f = 0.D0
  gupyy_f = 0.D0
  gupzz_f = 0.D0 
  cmax = 0.D0
  cmin = 0.D0
  lapm1_f = 0.D0
  shiftx_f = 0.D0
  shifty_f = 0.D0
  shiftz_f = 0.D0 
  gxx_f = 0.D0
  gxy_f = 0.D0
  gxz_f = 0.D0
  gyy_f = 0.D0
  gyz_f = 0.D0 
  gzz_f = 0.D0 
  phi_f = 0.D0

  dP_m = 0.D0
  Pr = 0.D0
  Pl = 0.D0
  drho_b_m = 0.D0
  rho_br = 0.D0
  rho_bl = 0.D0
  dvx_m = 0.D0
  dvy_m = 0.D0
  dvz_m = 0.D0
  vxr = 0.D0
  vxl = 0.D0
  vyr = 0.D0
  vyl = 0.D0
  vzr = 0.D0
  vzl = 0.D0
  cmax = 0.D0
  cmin = 0.D0
  Bxr = 0.D0
  Bxl = 0.D0
  Byr = 0.D0
  Byl = 0.D0
  Bzr = 0.D0
  Bzl = 0.D0
  v02r = 0.D0
  v02l = 0.D0
  dBx_m = 0.D0
  dBy_m = 0.D0
  dBz_m = 0.D0

  rho_star_rhs = 0.D0
  tau_rhs = 0.D0 
  st_x_rhs = 0.D0 
  st_y_rhs = 0.D0
  st_z_rhs = 0.D0 
  v02R = 0.D0
  v02l = 0.D0
  cmin = 0.D0
  cmax = 0.D0
  
  rhob_floor = 0.D0
  P_floor = 0.D0

  d0rho_b_m = 0.D0 
  d0P_m = 0.D0
  d0vx_m = 0.D0 
  d0vy_m = 0.D0
  d0vz_m = 0.D0 
  d02rho_b_m = 0.D0
  d02P_m = 0.D0
  d02vx_m = 0.D0
  d02vy_m = 0.D0
  d02vz_m = 0.D0
  d2rho_b_bck = 0.D0
  d2P_bck = 0.D0
  d2vx_bck = 0.D0
  d2vy_bck = 0.D0
  d2vz_bck = 0.D0 
  d2rho_b_fwd = 0.D0
  d2P_fwd = 0.D0
  d2vx_fwd = 0.D0
  d2vy_fwd = 0.D0
  d2vz_fwd = 0.D0
  
  !-----------------------------------------

  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_vs')

  !======================================
  ! Solve for lapse in maximal slicing case
  !======================================
!  if (Slicing == MAX_SLICING) {
!    Compute_Ricci(t0);
!    Lapse->Guess(t0);
!    Lapse->Predict_Lapse(t0,t0);
!  }
end subroutine bns_copy_to_prev_timelevel
