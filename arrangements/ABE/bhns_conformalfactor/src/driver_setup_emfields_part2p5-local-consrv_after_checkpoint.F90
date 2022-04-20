#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine BHNS_setup_emfield_part2p5_local_consrv_aftercp(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if(CCTK_ITERATION .eq. ITERATION_TO_INSERT_MAGNETIC_FIELDS .and. ITERATION_TO_INSERT_MAGNETIC_FIELDS.ne.0) then

     ! We are faced with a tough choice. If we add magnetic fields at a checkpoint, we need to fill in all the timelevels 
     !  appropriately. Magnetic fields depend on the pressure on those timelevels, but pressure is not stored on all 
     !  timelevels, so wee need to call the primitives solver (using conservatives on all the timelevels) and then
     !  calculate A's & then B's. The B's are then correct on all points except the stencil size, which is:
     !  1 (Pressure averaging to set A_phi) +
     !  1 (A_phi's -> Ai's) +
     !  1 (Ai's -> Bi's)
     !  -----------------
     !  3 = stencil size = # of ghostzones.
     ! So we need at least 3 ghostzones to properly set up the magnetic fields.

     if(1==0) then
        ! First do the _p_p timelevel.
        call BSSN_compute_gupij(cctkGH,cctk_lsh,             &
             gxx_p_p,gxy_p_p,gxz_p_p,gyy_p_p,gyz_p_p,gzz_p_p,&
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)

        call recompute_conservatives_fast_standalone_gf(cctkGH,cctk_lsh, &
             rho_b,P,vx,vy,vz, &
             phi_p_p,lapm1_p_p, &
             shiftx_p_p,shifty_p_p,shiftz_p_p, &
             gxx_p_p,gxy_p_p,gxz_p_p,gyy_p_p,gyz_p_p,gzz_p_p, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             Bx,By,Bz, &
             rho_star_p_p,tau_p_p,mhd_st_x_p_p,mhd_st_y_p_p,mhd_st_z_p_p)

        ! Then do the _p timelevel.
        call BSSN_compute_gupij(cctkGH,cctk_lsh, &
             gxx_p,gxy_p,gxz_p,gyy_p,gyz_p,gzz_p,&
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)

        call recompute_conservatives_fast_standalone_gf(cctkGH,cctk_lsh, &
             rho_b,P,vx,vy,vz, &
             phi_p,lapm1_p, &
             shiftx_p,shifty_p,shiftz_p, &
             gxx_p,gxy_p,gxz_p,gyy_p,gyz_p,gzz_p, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             Bx,By,Bz, &
             rho_star_p,tau_p,mhd_st_x_p,mhd_st_y_p,mhd_st_z_p)

        ! Finally, do the current timelevel. This way gupij are correctly set for this timelevel.
        call BSSN_compute_gupij(cctkGH,cctk_lsh, &
             gxx,gxy,gxz,gyy,gyz,gzz,            &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)

        call recompute_conservatives_fast_standalone_gf(cctkGH,cctk_lsh, &
             rho_b,P,vx,vy,vz, &
             phi,lapm1, &
             shiftx,shifty,shiftz, &
             gxx,gxy,gxz,gyy,gyz,gzz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             Bx,By,Bz, &
             rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z)
     end if
  end if

end subroutine BHNS_setup_emfield_part2p5_local_consrv_aftercp
