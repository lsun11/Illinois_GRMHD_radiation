#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!---------------------------------------------------------------------------------------
! gradually increase dt to lower mom. constraint violation near beginning of simulation
!---------------------------------------------------------------------------------------
subroutine change_dt(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8                                   :: dT,dXglob,dYglob,dZglob

  dXglob = cctk_delta_space(1)
  dYglob = cctk_delta_space(2)
  dZglob = cctk_delta_space(3)

  dT = CCTK_DELTA_TIME

  if(CCTK_ITERATION.ge.0) then
     dtfac = 0.1D0
  end if

  if(CCTK_ITERATION.ge.4) then
     dtfac = 0.125D0
  end if

  if(CCTK_ITERATION.ge.1024) then
     dtfac = 0.15D0
  end if

  if(CCTK_ITERATION.ge.1536) then
     dtfac = 0.175D0
  end if

  if(CCTK_ITERATION.ge.2048) then
     dtfac = 0.2D0
  end if

  if(CCTK_ITERATION.ge.2560) then
     dtfac = 0.25D0
  end if

  if(CCTK_ITERATION.ge.3072) then
     dtfac = 0.3D0
  end if

  if(CCTK_ITERATION.ge.3584) then
     dtfac = 0.4D0
  end if

  if(CCTK_ITERATION.ge.4096) then
     dtfac = 0.45D0
  end if

!  if(dtfac.ne.dt) then
!     CCTK_ITERATION = 0
!     CCTK_TIME = 0
!  end if

!!$  if(CCTK_ITERATION.eq.0) then
!!$     dtfac = 0.45*min(1.D0*CCTK_ITERATION/3072.D0,1.D0)
!!$  else 
!!$     dtfac = 0.45*min(1.D0/3072.D0,1.D0)
!!$  end if

  cctk_delta_time = max(dXglob,dYglob,dZglob)*dtfac

  write(*,*) "iter,dt_orig = ",CCTK_ITERATION,dT
  write(*,*) "iter,dtfac,dX/dT = ",CCTK_ITERATION,dtfac,CCTK_DELTA_TIME/dXglob
  write(*,*) "iter,dt,dtfac = ",CCTK_ITERATION,CCTK_DELTA_TIME,dtfac


end subroutine change_dt
