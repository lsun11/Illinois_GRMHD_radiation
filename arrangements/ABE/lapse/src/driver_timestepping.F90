!------------------------------------------
! Master timestepping driver for the lapse
!------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine lapse_timestepping(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                       :: dummy
  integer, dimension(3)                         :: ext
  real*8                                        :: dX,dY,dZ
!
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  if (CCTK_Equals(slicing_type,"frozen").eq.1 .or. cowling_enable==1) then
     call frozen_lapse_timestepping(CCTK_PASS_FTOF)
  else if (CCTK_Equals(slicing_type,"geodesic").eq.1) then
     call frozen_lapse_timestepping(CCTK_PASS_FTOF)
  else if (CCTK_Equals(slicing_type,"harmonic").eq.1) then
     call harmonic_lapse_timestepping(CCTK_PASS_FTOF)
  else if (CCTK_Equals(slicing_type,"hyperbolic").eq.1) then
     call hyperbolic_lapse_timestepping(CCTK_PASS_FTOF)
  else if (CCTK_Equals(slicing_type,"opl").eq.1) then
     call opl_lapse_timestepping(CCTK_PASS_FTOF)
  else if (CCTK_Equals(slicing_type,"opl_loglapse").eq.1) then
     call opl_loglapse_timestepping(CCTK_PASS_FTOF)
  end if
!  call CCTK_SyncGroup(dummy,cctkGH,'lapse::lapse_vars')
end subroutine Lapse_timestepping

subroutine frozen_lapse_timestepping(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  ! Set lapse rhs = zero
  !
  lapm1_rhs = 0.D0

end subroutine frozen_lapse_timestepping

subroutine harmonic_lapse_timestepping(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  real*8                        :: ONE,SIX
  parameter (ONE = 1.D0, SIX = 6.D0)
  
  ! Calculate lapse
  !
  lapm1_rhs = (lapm1 + ONE) * exp(SIX*(phi-phi_p)) - ONE
  
end subroutine harmonic_lapse_timestepping

