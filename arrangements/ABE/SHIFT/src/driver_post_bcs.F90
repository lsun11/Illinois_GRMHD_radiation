!--------------------------------
! Post-bc routines for the shift
!--------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine shift_postbc(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                       :: dummy
  integer, dimension(3)                         :: ext
  real*8                                        :: dX,dY,dZ, i,j,k
!
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)

!write(*,*) "Start shift_postbc!!!!!!"

  ext = cctk_lsh
  if(cctk_iteration.gt.0) then
     call fill_shift_symmetry_gz(ext,X,Y,Z,Symmetry,shiftx,shifty,shiftz,shiftxt,shiftyt,shiftzt)
  end if


end subroutine shift_postbc
