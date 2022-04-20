!-----------------------------------------------------------------------------
! Cowling approximation: set rhs's == 0
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine BSSN_timestepping_Cowling(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8                  :: dT,dX,dY,dZ, i,j,k
  integer, dimension(3)   :: ext
  integer                 :: dummy,index
  integer, parameter      :: AXISYM = 4

!
  dT = CCTK_DELTA_TIME
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)
  
  ext = cctk_lsh


!The following call might not be necessary, but I've added it here just in case...
  call BSSN_compute_gupij(cctkGH,cctk_lsh, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)

  phi_rhs = 0.D0
  chi_rhs = 0.D0
  trK_rhs = 0.D0
  gxx_rhs = 0.D0
  gxy_rhs = 0.D0
  gxz_rhs = 0.D0
  gyy_rhs = 0.D0
  gyz_rhs = 0.D0
  gzz_rhs = 0.D0
  Axx_rhs = 0.D0
  Axy_rhs = 0.D0
  Axz_rhs = 0.D0
  Ayy_rhs = 0.D0
  Ayz_rhs = 0.D0
  Azz_rhs = 0.D0
  Gammax_rhs = 0.D0
  Gammay_rhs = 0.D0
  Gammaz_rhs = 0.D0

!  lapsex = lapsex_p
!  lapsey = lapsey_p
!  lapsez = lapsez_p
! FOLLOWING HANDLED BY SHIFT THORN (sets shifti_rhs=0)
!  shiftx = shiftx_p
!  shifty = shifty_p
!  shiftz = shiftz_p

  !needed for Komar mass calculation:
  !  S = S_p

!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
!  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry)
!  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  
end subroutine BSSN_timestepping_Cowling
