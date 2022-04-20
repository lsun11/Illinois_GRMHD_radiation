!-----------------------------------
! Print the current iteration, time
!-----------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine scalarwaveMoL_print_time(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: k
  real*8 :: levelnumber

  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0

  write(*,*) "=========================================================="
  write(*,"(A9,I6,A6,I2,A22,F13.6)") ' Iter. # ',cctk_iteration, ', Lev: ',int(levelnumber),', Integrating to time: ',cctk_delta_time/cctk_levfac(1)+cctk_time
  write(*,*) "=========================================================="

  iter_count = 1

  do k=1,cctk_lsh(3)
     !write(*,*) "extents:",k,Z(1,1,k)
  end do

end subroutine scalarwaveMoL_print_time
