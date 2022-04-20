/*@@
  @file      AHFinder_initguess.F
  @date      July 1999
  @author    Miguel Alcubierre
  @desc 
             Find initial guess for minimization algorithm
  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      subroutine AHFinder_initguess(CCTK_ARGUMENTS,rmx)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j
      integer nn0,nn2,min_i,min_j

      CCTK_REAL zero,one,cc0,cc2,rmx,aux,min_aa
      CCTK_REAL, allocatable, dimension(:,:) :: aa,a0,a2
      
!     Description of variables
!
!     i,j       Counters.     
!
!     aa        Area array for initial guess.
!     a0        c0(0) array for initial guess.
!     a2        c0(2) array for initial guess.
!
!     rmx       Radius of largest sphere that fits in the grid.
!
!     nn0       Number of subdivisions of c0(0) for initial guess.
!     nn2       Number of subdivisions of c0(2) for initial guess.


!     **********************
!     ***   INITIALIZE   ***
!     **********************

!     define numbers

      zero = 0.0D0
      one  = 1.0D0

      aux = rmx

!     Initialize arrays.

      c0 = zero
      cc = zero
      cs = zero

      cc0 = zero
      cc2 = zero

!     Get {nn0,nn2}.

      nn0 = ahf_nn0
      nn2 = ahf_nn0

      if (cartoon) nn2 = 0

!     Allocate storage.

      allocate(aa(0:nn0,0:nn2))
      allocate(a0(0:nn0,0:nn2))
      allocate(a2(0:nn0,0:nn2))

      aa = zero; a0 = zero; a2 = zero


!     ****************************
!     ***  FIND INITIAL GUESS  ***
!     ****************************

!     IMPORTANT:  In order to choose a good initial guess,
!     I map a section of the space {c0(0),c0(2)} looking
!     for a local minimum.  This is VERY expensive, so the
!     parameter space is not very well convered.

      if (veryver.and.(myproc.eq.0)) then
         write(*,*)
         write(*,*) 'Choosing initial guess ...'
      end if

!     lmax >= 2 and good guess.

      if ((lmax.ge.2).and.(.not.sloppy)) then
         if (veryver.and.(myproc.eq.0)) then
            write(*,*)
            write(*,*) 'Choosing GOOD initial guess ...'
         end if

         do j=0,nn2

            c0(2) = (0.5D0*dble(j)/dble(nn2) - 0.25D0)*aux

            call AHFinder_fun(CCTK_ARGUMENTS)
            call AHFinder_exp(CCTK_ARGUMENTS)

            do i=0,nn0

               c0(0) = 4.0D0*dx + (aux-4.0D0*dx)*dble(i)/dble(nn0)

               a0(i,j) = c0(0)
               a2(i,j) = c0(2)

!     sqrt(5) is a normalization for P_lm's
               if(c0(0).gt.abs(c0(2)*sqrt(5.0))) then
                  call AHFinder_int(CCTK_ARGUMENTS)
                  
                  if (minarea) then
                     aa(i,j) = intarea
                  else
                     aa(i,j) = intexp2
                  end if
                  
                  if (guessverbose) then
                     write(*,*)
                     write(*,"(A12,2ES14.6)") ' a0 and a2: ',a0(i,0),a2(i,0)
                     write(*,"(A6,ES14.6)")   ' H  : ',intexp
                     write(*,"(A6,ES14.6)")   ' H^2: ',intexp2
                  end if
               endif   
            end do
         end do
         
         if ((.not.inner).or.(minarea)) then
            
            min_i = 1
            min_j = 1
            min_aa = 1.0d15
            
            do i=1,nn0-1
               do j=1,nn2-1
                  if ((aa(i,j).gt.zero).and. &
                       (aa(i,j).lt.aa(i+1,j)).and. &
                       (aa(i,j).lt.aa(i-1,j)).and. &
                       (aa(i,j).lt.aa(i,j+1)).and. &
                       (aa(i,j).lt.aa(i,j-1))) then
                     cc0 = a0(i,j)
                     cc2 = a2(i,j)
                  end if
                  if (aa(i,j) .le. min_aa) then
                     min_i = i
                     min_j = j
                     min_aa = aa(i,j)
                  end if
               end do
            end do
            
         else
            
            min_i = nn0-1
            min_j = 1
            min_aa = 1.0d15
            
            do i=nn0-1,1,-1
               do j=1,nn2-1
                  if ((aa(i,j).gt.zero).and. &
                       (aa(i,j).lt.aa(i+1,j)).and. &
                       (aa(i,j).lt.aa(i-1,j)).and. &
                       (aa(i,j).lt.aa(i,j+1)).and. &
                       (aa(i,j).lt.aa(i,j-1))) then
                     cc0 = a0(i,j)
                     cc2 = a2(i,j)
                  end if
                  if (aa(i,j).le.min_aa) then
                     min_i = i
                     min_j = j
                     min_aa = aa(i,j)
                  end if
               end do
            end do
            
         end if
         
         if (guess_absmin) then
            cc0 = a0(min_i,min_j)
            cc2 = a2(min_i,min_j)
         endif
         
!     Set up initial guess.
         
         c0(0) = cc0
         c0(2) = cc2
         
!     lmax < 2 or sloppy guess.
         
      else
         if (veryver.and.(myproc.eq.0)) then
            write(*,*)
            write(*,*) 'Choosing SLOPPY initial guess ...'
         end if
         
         call AHFinder_fun(CCTK_ARGUMENTS)
         call AHFinder_exp(CCTK_ARGUMENTS)
         
         do i=0,nn0
            
            c0(0) = 4.0D0*dx + (aux-4.0D0*dx)*dble(i)/dble(nn0)
            
            a0(i,0) = c0(0)
            if (veryver.and.(myproc.eq.0)) then
               write(*,*)
               write(*,*) 'About to call AHFinder_int...',c0(0)
            end if
           
            call AHFinder_int(CCTK_ARGUMENTS)

            if (minarea) then
               aa(i,0) = intarea
            else
               aa(i,0) = intexp2
            end if

            if (guessverbose) then
               write(*,*)
               write(*,"(A6,ES14.6)") ' a0 : ',a0(i,0)
               write(*,"(A6,ES14.6)") ' H  : ',intexp
               write(*,"(A6,ES14.6)") ' H^2: ',intexp2
            end if

         end do

         if ((.not.inner).or.(minarea)) then

            do i=1,nn0-1
               if ((aa(i,0).gt.zero).and. &
                    (aa(i,0).lt.0.9D0*aa(i+1,0)).and. &
                    (aa(i,0).lt.0.9D0*aa(i-1,0))) then
                  cc0 = a0(i,0)
               end if
            end do

         else

            do i=nn0-1,1,-1
               if ((aa(i,0).gt.zero).and. &
                    (aa(i,0).lt.0.9D0*aa(i+1,0)).and. &
          (aa(i,0).lt.0.9D0*aa(i-1,0))) then
                  cc0 = a0(i,0)
               end if
            end do

         end if

!        Set up initial guess.

         c0(0) = cc0

      end if

!     If no local minimum was found, start with
!     a large sphere anyway.

      if (c0(0).eq.zero) then
         c0    = zero
         c0(0) = rmx
      end if


!     ***************
!     ***   END   ***
!     ***************

!     Print message.

      if (veryver.and.(myproc.eq.0)) then
         write(*,*)
         write(*,*) 'Initial guess found'
         if (lmax.lt.2) then
            write(*,"(A8,ES14.6)") ' c0(0): ',c0(0)
         else
            write(*,"(A14,2ES14.6)") ' c0(0),c0(2): ',c0(0),c0(2)
         end if
      end if

!     Deallocate storage.

      deallocate(aa,a0,a2)


      end subroutine AHFinder_initguess

