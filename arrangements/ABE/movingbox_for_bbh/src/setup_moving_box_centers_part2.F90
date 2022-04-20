#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Setup movingbox center.  
!----------------------------------------------------------------------------------
subroutine setup_movingbox_centers_part2(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  real*8 :: levfac, rad, dist2_bh1, dist2_bh2, rad2
  integer :: i,j,k
!
  levfac = dble(cctk_levfac(1))

  rad = min( 0.25d0*abs(x0_BH1-x0_BH2), max_dev1, max_dev2)
  rad2 = rad*rad
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
	do i=1,cctk_lsh(1) 
           dist2_bh1 = (X(i,j,k)-x0_BH1)**2 + Y(i,j,k)**2 + Z(i,j,k)**2
	   dist2_bh2 = (X(i,j,k)-x0_BH2)**2 + Y(i,j,k)**2 + Z(i,j,k)**2
	   if (dist2_bh1 .lt. rad2) then 
	      tempx1(i,j,k) = levfac
  	   else 
	      tempx1(i,j,k) = 0.d0
	   end if
           if (dist2_bh2 .lt. rad2) then
              tempy1(i,j,k) = levfac
           else
              tempy1(i,j,k) = 0.d0
           end if
	end do
     end do
  end do

end subroutine setup_movingbox_centers_part2
