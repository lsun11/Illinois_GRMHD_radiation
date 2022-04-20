!/*@@
!  @file      InitSym_Om.F90
!  @date      August 2013
!  @author    Vasileios Paschalidis
!  @desc 
!  
!  @enddesc 
!@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


!/*@@
!  @routine    ang_freq_ofBfield_InitSym
!  @date       August 2013
!  @author     Vasileios Paschalidis
!  @desc 
!  
!  @enddesc 
!  @calls     
!  @calledby   
!  @history 
!
!  @endhistory 
!@@*/

subroutine ang_freq_ofBfield_InitSym(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  integer ones(3)
  integer ierr

  ones(1)=1
  ones(2)=1
  ones(3)=1

  call SetCartSymVN(ierr,cctkGH,ones,'ang_freq_ofBfield::Bfreq1')

  call SetCartSymVN(ierr,cctkGH,ones,'ang_freq_ofBfield::Bfreq2')

  ! Vasilis says: The symmetries above should be OK. May need to 
  ! revisit this though. 

  return
end subroutine ang_freq_ofBfield_InitSym


