!/*@@
!  @file      StartUp.F
!  @date      November 1999
!  @author    Gabrielle Allen
!  @desc 
!  
!  @enddesc 
!@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

!/*@@
!  @routine    InitSym.F
!  @date       November 1999
!  @author     Gabrielle Allen
!  @desc 
!  
!  @enddesc 
!  @calls     
!  @calledby   
!  @history 
!
!  @endhistory 
!@@*/

subroutine PsiKadelia_InitSym(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS

  integer ones(3)
  integer ierr

  ones(1)=1
  ones(2)=1
  ones(3)=1

  call SetCartSymVN(ierr,cctkGH,ones,'gw_extraction::psi0re')
  ones(3) = -1
  call SetCartSymVN(ierr,cctkGH,ones,'gw_extraction::psi0im')
  ones(3) = 1

  return
end subroutine PsiKadelia_InitSym


