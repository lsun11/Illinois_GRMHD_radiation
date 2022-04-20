!------------------------------------
! Timestepping routine for opl lapse
!------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine opl_lapse_timestepping(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS   

  interface
     subroutine opl_lapse_rhs(cctkGH, cctk_nghostzones,ext,Symmetry,opl_alap,f_of_alpha, &
          r,PhysicalRadius,RadiusDerivative, &
          lapse_old, lapse_rhs, &
          phi,trK,dX,dY,dZ)
       implicit none
       CCTK_POINTER                            :: cctkGH
       integer, dimension(3)                   :: cctk_nghostzones,ext
       integer                                 :: opl_alap,Symmetry,f_of_alpha
       real*8, dimension(ext(1),ext(2),ext(3)) :: r, PhysicalRadius, RadiusDerivative
       real*8, dimension(ext(1),ext(2),ext(3)) :: lapse_old,lapse_rhs
       real*8, dimension(ext(1),ext(2),ext(3)) :: phi,trK
       real*8                :: dX,dY,dZ
     end subroutine opl_lapse_rhs
  end interface
  integer               :: dummy,i
  integer, dimension(3) :: ext,fake_ext
  real*8                :: dT,dX,dY,dZ
  integer               :: PI_SYMM, AXISYM
  parameter(PI_SYMM = 3, AXISYM = 4)
  !
  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)
  dT = CCTK_DELTA_TIME

  call opl_lapse_rhs(cctkGH,cctk_nghostzones,cctk_lsh,Symmetry, opl_alap,f_of_alpha, &
       r,PhysicalRadius,RadiusDerivative, &
       lapm1, lapm1_rhs, phi, trK,dx,dy,dz)

  if(opl_advect_enable.ne.0) then
     call opl_lapse_upwind(cctkGH,cctk_nghostzones,cctk_lsh,Symmetry,opl_advect_enable, &
          dT,dX,dY,dZ, &
          lapsex,lapsey,lapsez, &
          shiftx,shifty,shiftz, &
          lapm1, lapm1_rhs)
     if(enable_lower_order_at_boundaries==1) then
        call opl_lapse_upwind_4(cctkGH,cctk_nghostzones,cctk_lsh,Symmetry,opl_advect_enable, &
             dT,dX,dY,dZ, &
             lapsex,lapsey,lapsez, &
             shiftx,shifty,shiftz, &
             lapm1, lapm1_rhs)
        call opl_lapse_upwind_2(cctkGH,cctk_nghostzones,cctk_lsh,Symmetry,opl_advect_enable, &
             dT,dX,dY,dZ, &
             lapsex,lapsey,lapsez, &
             shiftx,shifty,shiftz, &
             lapm1, lapm1_rhs)
     end if
  end if

  !  call CCTK_SyncGroup(dummy,cctkGH,'lapse::lapse_vars')
  !  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')

end subroutine opl_lapse_timestepping
