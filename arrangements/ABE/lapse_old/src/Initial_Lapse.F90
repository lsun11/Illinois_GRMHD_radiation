!-----------------------------------------------------------------------------
!
!$Id: Initial_Lapse.F90  $
!
!-----------------------------------------------------------------------------
!
! A function to initialize the Lapse
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
!
! Predictor
!
!-----------------------------------------------------------------------------
subroutine Setup_Initial_Lapse(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  write (*,*) "INSIDE SETUP INITIAL LAPSE!!"
  
  if (CCTK_Equals(slicing_type,"geodesic") .eq. 1) then
     write(*,*) "setting GEODESIC LAPSE!"
     call initial_lapse_geodesic(CCTK_PASS_FTOF)
  else if (CCTK_Equals(slicing_type,"frozen") .eq. 1) then
     write (*,*) "frozen lapse"
  else if (CCTK_Equals(slicing_type,"harmonic") .eq. 1) then
     call initial_lapse_harmonic(CCTK_PASS_FTOF)
  else if (CCTK_Equals(slicing_type,"opl") .eq. 1) then
     call initial_lapse_opl(CCTK_PASS_FTOF)
  else if (CCTK_Equals(slicing_type,"opl_loglapse") .eq. 1) then
     call initial_lapse_opl(CCTK_PASS_FTOF)
  end if
  
end subroutine Setup_Initial_Lapse

subroutine initial_lapse_geodesic(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !
  ! Set lapse to zero
  !
  lapm1 = 0.d0
  lapm1_p = 0.d0
  write (*,*) "SETTING LAPSE TO ZERO!"
end subroutine initial_lapse_geodesic

subroutine initial_lapse_harmonic(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
!  integer			:: Symmetry
  real*8                        :: ONE,SIX
  parameter (ONE = 1.D0, SIX = 6.D0)
  
  ! Calculate lapse
  !
  lapm1 = exp(SIX*(phi)) - ONE 
  lapm1_p = exp(SIX*(phi)) - ONE 
  
end subroutine initial_lapse_harmonic

subroutine initial_lapse_opl(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
!  integer			:: Symmetry
  real*8                        :: ONE
  parameter (ONE = 1.D0)
  
  ! Calculate lapse
  !
  if(opl_a0lap.ge.0) then
     write(*,*)'A0LAP:',opl_a0lap
     lapm1 = exp(-1*opl_a0lap*(phi)) - ONE 
     lapm1_p = exp(-1*opl_a0lap*(phi)) - ONE 
  end if


end subroutine initial_lapse_opl
