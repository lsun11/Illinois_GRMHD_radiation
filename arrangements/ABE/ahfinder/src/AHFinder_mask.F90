/*@@
  @file      AHFinder_mask.F
  @date      November 1998
  @author    Miguel Alcubierre
  @desc 
             Set up horizon mask for black hole excision.
  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine AHFinder_mask(CCTK_ARGUMENTS,rhor)

!     Set up the horizon mask.  All this routine does is
!     set up the value of the grid function "ahmask" to zero
!     for points inside the horizon.
!
!     It also finds out the radius of the largest sphere that
!     fits inside the horizon and it centered in (xc,yc,zc).
!     This is useful for simple excision.

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k
      integer reduce_handle

      CCTK_REAL zero,one
      CCTK_REAL rhor,rhortemp
      CCTK_REAL buffer,dd,aux
      CCTK_REAL xa,ya,za


!     *******************
!     ***   NUMBERS   ***
!     *******************

      zero = 0.0D0
      one  = 1.0D0


!     *********************************
!     ***   FIND HORIZON FUNCTION   ***
!     *********************************

      call AHFinder_fun(CCTK_ARGUMENTS)


!     ***********************************
!     ***   FIND BUFFER ZONE FACTOR   ***
!     ***********************************

!     The mask will always be set inside the horizon minus
!     a safety buffer zone.  To find out how big is the buffer
!     zone I check the value of the parameter "ahf_maskshrink".
!     This parameter says by which factor of the monopole term
!     to shrink the mask with respect to the real horizon.

!     Here I use the fact that the function "ahfgrid" is zero
!     at the horizon, and is linear on the monopole term
!     (which measures the overall scale).  What I then do is
!     set up the value of "buffer" to a negative number that
!     roughly measures the desired coordinate distance from
!     the mask to the horizon.  The reason why it is negative
!     is because the horizon function is negative inside the
!     horizon, and I just want to compare its value with the
!     value of "buffer" to decide when a given point is masked.

      buffer = - (one-ahf_maskshrink)*c0(0)

!     Now I make sure that there are at least "ahf_maskbuffer"
!     grid points in the buffer zone. The final buffer zone will
!     be whatever is largest between the number of grid points
!     specified by "ahf_maskbuffer" and the shrink factor
!     specified by "ahf_maskshrink".

      dd  = max(dx,dy,dz)
      aux = ahf_maskbuffer*dd

      if (abs(buffer).lt.aux) buffer = - aux


!     *****************************************
!     ***   FIND RADIUS OF LARGEST SPHERE   ***
!     *****************************************

!     Find the radius of the largest sphere that
!     fits inside the horizon.

      rhor = 1.0D10

      do k=1,nz
         do j=1,ny
            do i=1,nx

               if (ahfgrid(i,j,k).ge.buffer) then
                  aux = sqrt((x(i,j,k)-xc)**2 + (y(i,j,k)-yc)**2 &
                       + (z(i,j,k)-zc)**2)
                  if (rhor.gt.aux) rhor = aux
               end if

            end do
         end do
      end do

!     Now find the minimum of rhor across processors.

      reduce_handle = -1
      call CCTK_ReductionArrayHandle(reduce_handle,"minimum")
      if (reduce_handle.lt.0) then
        call CCTK_WARN(1,"Cannot get handle for minimum reduction ! Forgot to activate an implementation providing reduction operators ??")
      end if

      call CCTK_ReduceLocalScalar(ierr ,cctkGH,-1,reduce_handle, &
           rhor,rhortemp,CCTK_VARIABLE_REAL)
      rhor = rhortemp

      if (ierr.ne.0) then
         call CCTK_WARN(1,"AHFinder_mask:Reduction of rhor failed!")
      end if


!     ***********************
!     ***   LEGO SPHERE   ***
!     ***********************

!     If we want a lego shere mask life is simple.
!     Just set the mask inside of horizon (minus a
!     buffer zone) equal to 0.

      if (CCTK_EQUALS(ahf_masktype,'lego')) then
         do k=1,nz
            do j=1,ny
               do i=1,nx

                  if (ahfgrid(i,j,k).lt.buffer) then
                     ahmask(i,j,k) = zero
                  end if

               end do
            end do
         end do


!     ************************
!     ***   CUBICAL MASK   ***
!     ************************

!     If we want a cubical mask, set the mask based on the value
!     of rhor.  Notice that if we want the corners of the cube to
!     touch the sphere of radius rhor, we must make the side of
!     the cube equal to rhor/sqrt(3), which is about 0.57*rhor.

      else if (CCTK_EQUALS(ahf_masktype,'cube')) then

         do k=1,nz
            do j=1,ny
               do i=1,nx

                  aux = 0.57D0*rhor

                  xa = abs(x(i,j,k)-xc)
                  ya = abs(y(i,j,k)-yc)
                  za = abs(z(i,j,k)-zc)

                  if ((xa.lt.aux).and.(ya.lt.aux).and.(za.lt.aux)) then
                     ahmask(i,j,k) = zero
                  end if

               end do
            end do
         end do


!     **************************
!     ***   POLYHEDRA MASK   ***
!     **************************

!     Here we want to excise a polyhedra.  But not just any
!     polyhedra, a polyhedra defined by taking a cube and
!     cutting away the corners all the way up to the middle
!     of the edges.  This results in a figure made up of
!     6 squares and 8 equilateral triangles.

!     Notice that for such a figure, if we want the corners
!     to touch a sphere of radius rhor, then the sides have
!     to be rhor/sqrt(2) which is roughly 0.7*rhor.

      else if (CCTK_EQUALS(ahf_masktype,'poly')) then

         do k=1,nz
            do j=1,ny
               do i=1,nx

                  aux = 0.7D0*rhor

                  xa = abs(x(i,j,k)-xc)
                  ya = abs(y(i,j,k)-yc)
                  za = abs(z(i,j,k)-zc)

                  if ((xa.lt.aux).and.(ya.lt.aux).and.(za.lt.aux).and. &
                       ((xa+ya+za).lt.2.0D0*aux)) then
                     ahmask(i,j,k) = zero
                  end if

               end do
            end do
         end do

      end if


!     *****************************************************
!     ***   APPLY SYMMETRY BOUNDARIES AND SYNCHRONIZE   ***
!     *****************************************************

      call CCTK_SyncGroup(ierr,cctkGH,"ahfinder::ahfmask")
      call CartSymGN(ierr,cctkGH,"ahfinder::ahfmask")

!     If desired, copy ahf_mask to spacemask mask (emask).

      if (use_mask.eq.1) then
         emask = ahmask
      end if


!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder_mask

