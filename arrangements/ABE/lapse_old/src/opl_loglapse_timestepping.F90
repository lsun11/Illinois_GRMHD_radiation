!------------------------------------
! Timestepping routine for opl lapse
!------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine opl_loglapse_timestepping(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS   

  interface
     subroutine opl_loglapse_rhs(cctkGH, cctk_nghostzones,ext,Symmetry, &
          lapse_rhs, trK)
       implicit none
       CCTK_POINTER                            :: cctkGH
       integer, dimension(3)                   :: cctk_nghostzones,ext
       integer                                 :: Symmetry
       real*8, dimension(ext(1),ext(2),ext(3)) :: lapse_rhs,trK
     end subroutine opl_loglapse_rhs
  end interface
  integer               :: dummy,i,j,k
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

  !  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  !  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_aux_restrict')

  !$omp parallel do
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           lapm1(i,j,k) = log(1.D0 + lapm1(i,j,k))
        end do
     end do
  end do
  !$omp end parallel do

  !  lapm1 = log(1.D0 + lapm1)

  call opl_loglapse_rhs(cctkGH,cctk_nghostzones,cctk_lsh,Symmetry, &
       lapm1_rhs, trK)

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

  if(opl_advect_enable==1) then
     write(*,*) "opl_advect_enable==1 NOT SUPPORTED with opl_loglapse! (need to fix lapsex,lapsey,lapsez)"
     stop
  end if

  !  call CCTK_SyncGroup(dummy,cctkGH,'lapse::lapse_vars')
  !  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')

  !$omp parallel do
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           lapm1(i,j,k) = exp(lapm1(i,j,k)) - 1.D0
           lapm1_rhs(i,j,k) = (lapm1(i,j,k) + 1.D0) * lapm1_rhs(i,j,k)
        end do
     end do
  end do
  !$omp end parallel do

!  call hbpuncture_add_dissipation(cctk_nghostzones,cctk_lsh,dx,dy,dz,lapm1_rhs,lapm1)

  !  lapm1 = exp(lapm1) - 1.D0
  !  lapm1_rhs = (lapm1 + 1.D0) * lapm1_rhs

end subroutine opl_loglapse_timestepping
