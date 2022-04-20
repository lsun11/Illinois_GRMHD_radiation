!--------------------------------------------------------------------------------
! Timestepping routines for first-order Gamma-driving shift (Spatial_Gauge == 6)
!--------------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine puncture_firstorder_shift_timestepping(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  interface
     subroutine puncture_firstorder_shiftRHS(cctkGH, &
          nghostzones,ext,Symmetry, &
          dx,dy,dz,x,y,z,r, &
          PhysicalRadius,RadiusDerivative,RadiusDerivative2, &
          Gammax,Gammay,Gammaz, &
          gupxx,gupxy,gupxz, &
          gupyy,gupyz,gupzz,  &
          betax,betay,betaz, &
          betax_rhs,betay_rhs,betaz_rhs, &
          eta,firstorder_shift_convert_Gammai_fisheye_to_physical)
       implicit none
       CCTK_POINTER                            :: cctkGH
       integer, dimension(3)                   :: ext,nghostzones
       real*8, dimension(ext(1),ext(2),ext(3)) :: x,y,z,r,PhysicalRadius
       real*8, dimension(ext(1),ext(2),ext(3)) :: RadiusDerivative,RadiusDerivative2
       real*8, dimension(ext(1),ext(2),ext(3)) :: Gammax,Gammay,Gammaz
       real*8, dimension(ext(1),ext(2),ext(3)) :: gupxx,gupxy,gupxz
       real*8, dimension(ext(1),ext(2),ext(3)) :: gupyy,gupyz,gupzz
       real*8, dimension(ext(1),ext(2),ext(3)) :: betax,betay,betaz
       real*8, dimension(ext(1),ext(2),ext(3)) :: betax_rhs,betay_rhs,betaz_rhs
       real*8                                  :: dx,dy,dz,eta
       integer                                 :: firstorder_shift_convert_Gammai_fisheye_to_physical
       integer :: Symmetry
     end subroutine puncture_firstorder_shiftRHS
  end interface
  integer                                   :: n1,n2,n3,m,i,dummy,index
  integer, dimension(3)                     :: ext,fake_ext
  real*8                                    :: dT,dX,dY,dZ
  real*8                                    :: HALF,ld_eps,ld_c
!  real*8, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))   :: tempx,tempy,tempz
  integer                    :: PI_SYMM, AXISYM
  parameter(PI_SYMM = 3, AXISYM = 4)
  parameter ( HALF = 0.5D0 )
  !
  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)
  dT = CCTK_DELTA_TIME

  call puncture_firstorder_shiftRHS(cctkGH, &
       cctk_nghostzones,cctk_lsh,Symmetry, &
       dx,dy,dz,x,y,z,r, &
       PhysicalRadius,RadiusDerivative,RadiusDerivative2, &
       Gammax,Gammay,Gammaz, &
       gupxx,gupxy,gupxz, &
       gupyy,gupyz,gupzz,  &
       shiftx,shifty,shiftz, &
       shiftx_rhs,shifty_rhs,shiftz_rhs, &
       eta,firstorder_shift_convert_Gammai_fisheye_to_physical)

  if(hbpunc_advect_enable .eq. 2) then
     call puncture_firstorder_shift_upwind(cctkGH, &
          cctk_nghostzones,cctk_lsh,Symmetry, &
          dx,dy,dz, &
          shiftx,shifty,shiftz, &
          shiftx_rhs,shifty_rhs,shiftz_rhs)
  end if

end subroutine puncture_firstorder_shift_timestepping
