/*@@
  @file      AHFinder_min.F
  @date      April 1998
  @author    Miguel Alcubierre
  @desc 
             Minimization routine.
  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      subroutine AHFinder_min(CCTK_ARGUMENTS,NN,status,logf)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      logical status,found

      integer i,NN
      integer l,m
      integer iter,itmax

      CCTK_REAL FTOL,FRET
      CCTK_REAL zero,half,one

      CCTK_REAL, dimension (NN)    :: PP
      CCTK_REAL, dimension (NN,NN) :: XI

      character(len=200) :: logf

!     Description of variables:
!
!     i,j,l,m   Counters.
!
!     ITER      Number of iterations taken in minimization.
!     ITMAX     Maximum number of iterations allowed.
!
!     FTOL      Fractional tolerance for return criterion.
!
!     FRET      Value of fucntion on return from minimization.
!
!     PP         Vector of coefficients.
!
!     XI        Matrix of initial directions.


!     ***********************************
!     ***   DEFINE  {zero,half,one}   ***
!     ***********************************

      zero = 0.0D0
      half = 0.5D0
      one  = 1.0D0


!     *********************************
!     ***   FIND INITIAL VECTOR PP   ***
!     *********************************

      PP = zero
      i = 0
 
!     Copy {xc,yc,zc}.

      if (wander) then
         if (.not.refx) then
            i = i+1
            PP(1) = xc
         end if
         if (.not.refy) then
            i = i+1
            PP(2) = yc
         end if
         if (.not.refz) then
            i = i+1
            PP(3) = zc
         end if
      end if

!     Copy c0(l).

      do l=0,lmax,1+stepz
         i = i+1
         PP(i) = c0(l)
      end do

!     Copy {cc,cs}.

      if (nonaxi) then

!        Copy cc(l,m).

         do l=1,lmax
            do m=1+stepx,l,1+stepx
               if (stepz*mod(l-m,2).eq.0) then
                  i = i+1
                  PP(i) = cc(l,m)
               end if
            end do
         end do

!        Copy cs(l,m).

         if (.not.refy) then
            do l=1,lmax
               do m=1,l,1+stepx
                  if (stepz*mod(l-m,2).eq.0) then
                     i = i+1
                     PP(i) = cs(l,m)
                  end if
               end do
            end do
         end if

      end if


!     *********************************************
!     ***   FIND MATRIX OF INITIAL DIRECTIONS   ***
!     *********************************************

!     The initial set of directions are taken to be
!     the unit vectors (the first index identifies the
!     components, and second index the vectors).

      XI = zero

      do i=1,NN
         XI(i,i) = one
      end do


!     ****************************************
!     ***   CALL MINIMIZATION SUBROUTINE   ***
!     ****************************************

!     Find ITMAX.

      ITMAX = ahf_maxiter

!     Find FTOL.  FTOL is the fractional tolerance in the
!     function value so that failure to decrease by more
!     than this amount on one interation signals that we
!     are done.

      FTOL = ahf_tol

!     Here I call the routine POWELL from Numerical Recipies.

      found = .false.
      firstfun = .true.

      call POWELL(CCTK_ARGUMENTS,PP,XI,NN,FTOL,ITER,ITMAX,FRET,found)


!     ****************************
!     ***   FIND  {c0,cc,cs}   ***
!     ****************************

      i = 0

!     Find {xc,yc,zc}.

      if (wander) then
         if (.not.refx) then
            i = i+1
            xc = PP(1)
         end if
         if (.not.refy) then
            i = i+1
            yc = PP(2)
         end if
         if (.not.refz) then
            i = i+1
            zc = PP(3)
         end if
      end if

!     Find c0(l).

      do l=0,lmax,1+stepz
         i = i+1
         c0(l) = PP(i)
      end do

!     Find {cc,cs}.

      if (nonaxi) then

!        Find cc(l,m).

         do l=1,lmax
            do m=1+stepx,l,1+stepx
               if (stepz*mod(l-m,2).eq.0) then
                  i = i+1
                  cc(l,m) = PP(i)
               end if
            end do
         end do

!        Find cs(l,m).

         if (.not.refy) then
            do l=1,lmax
               do m=1,l,1+stepx
                  if (stepz*mod(l-m,2).eq.0) then
                     i = i+1
                     cs(l,m) = PP(i)
                  end if
               end do
            end do
         end if

      end if


!     ***************************************
!     ***   DECIDE ON THE RETURN STATUS   ***
!     ***************************************

!     If a minimum was found, set status to true.

      if (found) then

         status = .true.

!     If no minimum was found, it might still be because
!     there were to many iterations, in which case there
!     is a minimum, we just did not converge to it in the
!     maximum number of iterations allowed.

      else

         if (ITER.eq.ITMAX) then

            status = .true.

            if ((verbose.or.logfile).and.(myproc.eq.0)) then

               write(*,*)
               write(*,*) 'AHFinder:  Too many iterations.'

               if (logfile) then
                  open(11,file=logf,form='formatted',
     .            status='old',position='append')
                  write(11,*)
                  write(11,*) 'AHFinder:  Too many iterations.'
                  close(11)
               end if

            end if

         else

            status = .false.

         end if

      end if


!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder_min
