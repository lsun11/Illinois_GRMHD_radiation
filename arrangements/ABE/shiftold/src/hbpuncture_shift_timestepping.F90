!-----------------------------------------------------------------------------
! Timestepping routine for 2nd order Gamma-driving shift (Spatial_Gauge == 5)
!-----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine hbpuncture_shift_timestepping(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  interface
     subroutine hbpuncture_shiftRHS(cctkGH,cctk_nghostzones,ext, &
          dT,dx,dy,dz,hbpunc_advect_enable,eta, &
          RadiusDerivative, &
          shiftx_old,shifty_old,shiftz_old, &
          shiftx_rhs,shifty_rhs,shiftz_rhs, &
          dtshiftx_old,dtshifty_old,dtshiftz_old, &
          dtshiftx_rhs,dtshifty_rhs,dtshiftz_rhs, &
          Gammax,Gammay,Gammaz, &
          Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
          shiftxt_timederiv,shiftyt_timederiv,shiftzt_timederiv,Symmetry,  &
          hbpuncture_shift_convert_Gammai_fisheye_to_physical, &
          r,eta_final_value,eta_falloff_enable,eta_falloff_radius,eta_falloff_dr, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,phi)
       implicit none
       CCTK_POINTER                            :: cctkGH
       integer, dimension(3)                   :: ext,cctk_nghostzones
       real*8                                  :: dT,dx,dy,dz,eta,eta_final_value,eta_falloff_radius,eta_falloff_dr
       real*8, dimension(ext(1),ext(2),ext(3)) :: RadiusDerivative,r
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftx_old,shifty_old,shiftz_old
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftx_rhs,shifty_rhs,shiftz_rhs
       real*8, dimension(ext(1),ext(2),ext(3)) :: dtshiftx_old,dtshifty_old,dtshiftz_old
       real*8, dimension(ext(1),ext(2),ext(3)) :: dtshiftx_rhs,dtshifty_rhs,dtshiftz_rhs
       real*8, dimension(ext(1),ext(2),ext(3)) :: Gammax,Gammay,Gammaz
       real*8, dimension(ext(1),ext(2),ext(3)) :: Gammax_rhs,Gammay_rhs,Gammaz_rhs
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftxt_timederiv,shiftyt_timederiv,shiftzt_timederiv
       real*8, dimension(ext(1),ext(2),ext(3)) :: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,phi
       integer                                 :: Symmetry,hbpunc_advect_enable,eta_falloff_enable
       integer                                 :: hbpuncture_shift_convert_Gammai_fisheye_to_physical
     end subroutine hbpuncture_shiftRHS
  end interface
  integer                                   :: n1,n2,n3,m,i,dummy,index
  integer, dimension(3)                     :: ext,fake_ext
  real*8                                    :: dT,dX,dY,dZ
  real*8                                    :: HALF,ld_eps,ld_c
!!!!  real*8, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))   :: tempx,tempy,tempz
  integer                    :: PI_SYMM, AXISYM
  parameter(PI_SYMM = 3, AXISYM = 4)
  parameter ( HALF = 0.5D0 )
  !
  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)
  dT = CCTK_DELTA_TIME

  ! WE USE shiftxt_timederiv,shiftyt_timederiv,shiftzt_timederiv AS TEMPORARY STORAGE HERE!

  if(CCTK_TIME.ge.time_to_switch_to_eta_final_value) then
     eta = eta_final_value
  end if

  write(*,*) "Inside hbpuncture_shift_timestepping: eta=",eta

  call hbpuncture_shiftRHS(cctkGH,cctk_nghostzones,cctk_lsh, &
       dT,dx,dy,dz,hbpunc_advect_enable,eta, &
       RadiusDerivative, &
       shiftx,shifty,shiftz, &
       shiftx_rhs,shifty_rhs,shiftz_rhs, &
       shiftxt,shiftyt,shiftzt, &
       shiftxt_rhs,shiftyt_rhs,shiftzt_rhs, &
       Gammax,Gammay,Gammaz, &
       Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
       shiftxt_timederiv,shiftyt_timederiv,shiftzt_timederiv,Symmetry,  &
       hbpuncture_shift_convert_Gammai_fisheye_to_physical, &
       r,eta_final_value,eta_falloff_enable,eta_falloff_radius,eta_falloff_dr, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,phi)

  if(enable_lower_order_at_boundaries==1) then
     call hbpuncture_shiftRHS_2(cctkGH,cctk_nghostzones,cctk_lsh, &
          dT,dx,dy,dz,hbpunc_advect_enable,eta, &
          RadiusDerivative, &
          shiftx,shifty,shiftz, &
          shiftx_rhs,shifty_rhs,shiftz_rhs, &
          shiftxt,shiftyt,shiftzt, &
          shiftxt_rhs,shiftyt_rhs,shiftzt_rhs, &
          Gammax,Gammay,Gammaz, &
          Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
          shiftxt_timederiv,shiftyt_timederiv,shiftzt_timederiv,Symmetry,  &
          hbpuncture_shift_convert_Gammai_fisheye_to_physical)
  end if

  if(hbpunc_advect_enable .ne. 0) then
     call hbpuncture_upwind(cctkGH,cctk_nghostzones,cctk_lsh, dT,  dx,  dy,  dz, &
          RadiusDerivative, hbpunc_advect_enable, bssn_enable_shift_upwind, &
          shiftx,shifty,shiftz, &
          shiftx_rhs,shifty_rhs,shiftz_rhs, &
          shiftxt,shiftyt,shiftzt, &
          shiftxt_rhs,shiftyt_rhs,shiftzt_rhs, &
          Gammax,Gammay,Gammaz, &
          shiftxt_timederiv,shiftyt_timederiv,shiftzt_timederiv,Symmetry)
     if(enable_lower_order_at_boundaries==1) then
        call hbpuncture_upwind_4(cctkGH,cctk_nghostzones,cctk_lsh, dT,  dx,  dy,  dz, &
             RadiusDerivative, hbpunc_advect_enable, bssn_enable_shift_upwind, &
             shiftx,shifty,shiftz, &
             shiftx_rhs,shifty_rhs,shiftz_rhs, &
             shiftxt,shiftyt,shiftzt, &
             shiftxt_rhs,shiftyt_rhs,shiftzt_rhs, &
             Gammax,Gammay,Gammaz, &
             shiftxt_timederiv,shiftyt_timederiv,shiftzt_timederiv,Symmetry)
        call hbpuncture_upwind_2(cctkGH,cctk_nghostzones,cctk_lsh, dT,  dx,  dy,  dz, &
             RadiusDerivative, hbpunc_advect_enable, bssn_enable_shift_upwind, &
             shiftx,shifty,shiftz, &
             shiftx_rhs,shifty_rhs,shiftz_rhs, &
             shiftxt,shiftyt,shiftzt, &
             shiftxt_rhs,shiftyt_rhs,shiftzt_rhs, &
             Gammax,Gammay,Gammaz, &
             shiftxt_timederiv,shiftyt_timederiv,shiftzt_timederiv,Symmetry)
     end if
  end if
end subroutine hbpuncture_shift_timestepping
