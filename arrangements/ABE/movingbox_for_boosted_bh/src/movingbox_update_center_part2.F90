#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Update center of a moving box.  
!----------------------------------------------------------------------------------
subroutine movingbox_update_center_part2(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i,j,k,istart,iend,jstart,jend,kstart,kend
!
  if (mod(cctk_iteration,out_every)==0 .and. cctk_iteration .gt. 0) then
     if (find_max_phi==1) then 

        istart = cctk_nghostzones(1)
        jstart = cctk_nghostzones(2)
        kstart = cctk_nghostzones(3)
        iend = cctk_lsh(1) - istart
        jend = cctk_lsh(2) - jstart
        kend = cctk_lsh(3) - kstart

        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 temp0(i,j,k) = 0.d0
                 tempx(i,j,k) = 0.d0
                 tempy(i,j,k) = 0.d0
              end do
           end do
        end do
        !$omp end parallel do

        ! Search for maximum phi at the equatorial plane
        !$omp parallel do
        do k=kstart,kend
           do j=jstart,jend
              do i=istart,iend
                 if (abs(phi(i,j,k) - max_phi) .lt. 1.d-10) then 
	    	    temp0(i,j,k) = 1.d0
	            tempx(i,j,k) = X(i,j,k)
		    tempy(i,j,k) = Y(i,j,k)
                 end if
              end do
           end do
        end do
        !$omp end parallel do

     end if
   end if

end subroutine movingbox_update_center_part2
