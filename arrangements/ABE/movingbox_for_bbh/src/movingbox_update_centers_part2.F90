#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Update center of a moving box.  
!----------------------------------------------------------------------------------
subroutine movingbox_update_centers_part2(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i,j,k,istart,iend,jstart,jend,kstart,kend,levfac
  real*8 :: max_phi, imax,jmax,kmax
  real*8 :: dist2, rad2
  logical :: found_max
!
  if (mod(cctk_iteration,out_every)==0 .and. cctk_iteration .gt. 0) then

     levfac = cctk_levfac(1)

     if (find_max_phi1==1) then 

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
                 tempx1(i,j,k) = -1.d100
                 tempy1(i,j,k) = -1.d100
              end do
           end do
        end do
        !$omp end parallel do

        ! Search for maxima phi at the equatorial plane
        max_phi = -1.d308
        found_max = .false.
        rad2 = max_dev1*max_dev1
        do k=kstart,kend
           do j=jstart,jend
              do i=istart,iend
	         dist2 = (X(i,j,k)-position_x(1))**2 + (Y(i,j,k)-position_y(1))**2 + Z(i,j,k)**2
		 if (dist2 .lt. rad2 .and. levfac==max_levfac1) then 
		    if (phi(i,j,k) .gt. max_phi) then 
		       found_max = .true.
		       max_phi = phi(i,j,k)
		       imax = i
		       jmax = j
		       kmax = k
		    end if
		 end if
              end do
           end do
        end do

        if (found_max) then 
	   tempx1(imax,jmax,kmax) = X(imax,jmax,kmax)
	   tempy1(imax,jmax,kmax) = Y(imax,jmax,kmax)
	end if

     end if

     if (find_max_phi2==1) then 

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
                 tempx2(i,j,k) = -1.d100
                 tempy2(i,j,k) = -1.d100
              end do
           end do
        end do
        !$omp end parallel do

        ! Search for maxima phi at the equatorial plane
        max_phi = -1.d308
        found_max = .false.
        rad2 = max_dev2*max_dev2
        do k=kstart,kend
           do j=jstart,jend
              do i=istart,iend
	         dist2 = (X(i,j,k)-position_x(2))**2 + (Y(i,j,k)-position_y(2))**2 + Z(i,j,k)**2
		 if (dist2 .lt. rad2 .and. levfac==max_levfac2) then 
		    if (phi(i,j,k) .gt. max_phi) then 
		       found_max = .true.
		       max_phi = phi(i,j,k)
		       imax = i
		       jmax = j
	               kmax = k
		    end if
		 end if
              end do
           end do
        end do

        if (found_max) then 
	   tempx2(imax,jmax,kmax) = X(imax,jmax,kmax)
	   tempy2(imax,jmax,kmax) = Y(imax,jmax,kmax)
	end if

     end if

   end if

end subroutine movingbox_update_centers_part2
