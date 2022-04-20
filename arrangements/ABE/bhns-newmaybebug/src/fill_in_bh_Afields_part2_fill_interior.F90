#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Now that we have values of A^i on two spheres, we can smoothly extrapolate A^i into the interior of the BHs and compute B^i.
!-----------------------------------------------------------------------------
subroutine fill_in_bh_Afields_part2_fill_interior(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                  :: i
  real*8	  		           :: bhx,bhy,bhz
  real*8,dimension(100)                    :: r_array

  if(MOD(cctk_iteration,refill_horizons_magfields_every)==0 .and. cctk_iteration.gt.0) then

     bhx = bh_posn_x(1)
     bhy = bh_posn_y(1)
     bhz = bh_posn_z(1)

     do i=1,RADIAL_INTERP_ORDER
        r_array(i)=horiz_radius*bhns_bh_filling_radius(i)
     end do
     if(NUM_ZERO_PTS.eq.0) r_array(100)=horiz_radius*bhns_bh_filling_radius(100)

     call bhns_fill_Ai_BH_interior(cctkGH,cctk_lsh,Symmetry, &
          Ax, X,Y,Z, &
          bhx,bhy,bhz, &
          NUM_ZERO_PTS,RADIAL_INTERP_ORDER, &
          INPUTARRAY_THETASIZE,INPUTARRAY_PHISIZE, &
          r_array, &
          Axintr)
     call bhns_fill_Ai_BH_interior(cctkGH,cctk_lsh,Symmetry, &
          Ay, X,Y,Z, &
          bhx,bhy,bhz, &
          NUM_ZERO_PTS,RADIAL_INTERP_ORDER, &
          INPUTARRAY_THETASIZE,INPUTARRAY_PHISIZE, &
          r_array, &
          Ayintr)

     call bhns_fill_Ai_BH_interior(cctkGH,cctk_lsh,Symmetry, &
          Az, X,Y,Z, &
          bhx,bhy,bhz, &
          NUM_ZERO_PTS,RADIAL_INTERP_ORDER, &
          INPUTARRAY_THETASIZE,INPUTARRAY_PHISIZE, &
          r_array, &
          Azintr)

     call bhns_fill_Ai_BH_interior(cctkGH,cctk_lsh,Symmetry, &
          rho_b, X,Y,Z, &
          bhx,bhy,bhz, &
          NUM_ZERO_PTS,RADIAL_INTERP_ORDER, &
          INPUTARRAY_THETASIZE,INPUTARRAY_PHISIZE, &
          r_array, &
          rhobintr)

     call bhns_fill_Ai_BH_interior(cctkGH,cctk_lsh,Symmetry, &
          P, X,Y,Z, &
          bhx,bhy,bhz, &
          NUM_ZERO_PTS,RADIAL_INTERP_ORDER, &
          INPUTARRAY_THETASIZE,INPUTARRAY_PHISIZE, &
          r_array, &
          Pintr)

     where(rho_b .lt. rho_b_atm) 
        rho_b = rho_b_atm
     end where

     where(P .lt. 0.D0)
        P = 0.D0
     end where

     !The following subroutine applies symmetry BCs on A^i and computes B^i from A^i:
     call bhns_compute_B_from_A(CCTK_PASS_FTOF)

  end if
end subroutine fill_in_bh_Afields_part2_fill_interior
