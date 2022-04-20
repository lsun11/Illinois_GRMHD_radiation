/*@@
  @file      AHFinder_pow.F
  @date      April 1998
  @author    Miguel Alcubierre
  @desc 
             Multidimensional minimization routines taken
             from Numerical Recipes.  The routines have some
             modifications to adapt them to the specific problem
             we are trying to solve.
  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


      module F1COM

!     Data to be shared by minimization routines.

      implicit none

      logical error

      integer NCOM

      CCTK_REAL XMIN
      CCTK_REAL, allocatable, dimension(:) :: PCOM,XICOM


!     ***************
!     ***   END   ***
!     ***************

      end module F1COM




!***********************************************************************
!***********************************************************************




      subroutine POWELL(CCTK_ARGUMENTS,PP,XI,N,FTOL,ITER,ITMAX,FRET,found)

!     Minimization of a function FUNC on N variables (FUNC is
!     not an argument, it is a fixed function name).  Input consists
!     of an initial starting point PP that is a vector of length N;
!     an initial matrix XI of dimensions NxN and whose columns contain
!     the initial set of directions (usually the N unit vectors); and
!     FTOL, the fractional tolerance in the function value such that
!     failure to decrease by more than this ammount on one iteration
!     signals doneness.  On output, PP is set to the best point found,
!     XI is the then-current direction set, FRET is the returned
!     function value at PP, and ITER is the number of iterations taken.
!
!     Extra output added by myself (Miguel Alcubierre) is a logical
!     variable "found" that is set to true is a minimum was indeed
!     found, and to false if it was not.

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS

      logical found

      integer i,j,N
      integer ITER,ITMAX,IBIG

      CCTK_REAL FUNC,DEL,FP,FPTT,FRET,FTOL,T,TOL1,ZEPS

      CCTK_REAL, dimension(1:N) :: PP,PT,PTT,XIT
      CCTK_REAL, dimension(1:N,1:N) :: XI

      character(len=200) :: logf

      parameter (ZEPS=1.0D-10)


!     *************************
!     ***   START ROUTINE   ***
!     *************************

      FRET = FUNC(CCTK_ARGUMENTS,PP,N)

!     Save initial point.

      do j=1,N
         PT(j) = PP(j)
      end do

      ITER = 0

   10 ITER = ITER + 1

      if ((myproc.eq.0).and.veryver) then
         write(*,*)
         write(*,"(A20,I2)") 'POWELL ITERATION  = ',ITER
      end if

      if ((myproc.eq.0).and.logfile) then
         logf = filestr(1:nfile)//"/ahf_logfile"
         open(11,file=logf,form='formatted',status='old', &
              position='append')
         write(11,*)
         write(11,"(A20,I2)") 'POWELL ITERATION  = ',ITER
         close(11)
      end if

      FP = FRET
      IBIG = 0

!     Initialize DEL.  It will be the biggest function decrease.

      DEL = 0.0D0

!     In each iteration, loop over all directions in the set.

      do i=1,N

!        Copy the direction.

         do j=1,N
            XIT(j) = XI(j,i)
         end do

         FPTT = FRET

!        Minimize along it.

         call LINMIN(CCTK_ARGUMENTS,PP,XIT,N,FRET,FTOL)

!        And record it if it is the largest decrease so far.

         if (abs(FPTT-FRET).gt.DEL) then
            DEL  = abs(FPTT-FRET)
            IBIG = I
         end if

   20   continue

      end do

!     Termination criterion.

      TOL1 = 0.5D0*FTOL*(abs(FP)+abs(FRET)) + ZEPS

      if (abs(FP-FRET).le.TOL1) then
         found = .true.
         return
      else if (abs(FRET).lt.ZEPS) then
         found = .true.
         return
      end if

!     Too many iterations?

      if (ITER.ge.ITMAX) return

!     Construct the extrapolated point and the average direction.

      do j=1,N
         PTT(j) = 2.0D0*PP(j) - PT(j)
         XIT(j) = PP(j) - PT(j)
         PT(j)  = PP(j)
      end do

!     Function value at extrapolated point.

      FPTT = FUNC(CCTK_ARGUMENTS,PTT,N)

!     First reason not to use new direction.

      if (FPTT.ge.FP) goto 10

!     Second reason not to use new direction.

      T = 2.0D0*(FP-2.0D0*FRET+FPTT)*(FP-FRET-DEL)**2 &
           - DEL*(FP-FPTT)**2

      if (T.gt.0.0D0) goto 10

!     Move to the minimum of the new direction.

      call LINMIN(CCTK_ARGUMENTS,PP,XIT,N,FRET,FTOL)

!     Save new direction.

      do j=1,N
         XI(j,IBIG) = XIT(j)
      end do

!     Back for another iteration.

      goto 10


!     ***************
!     ***   END   ***
!     ***************

      end subroutine POWELL




!***********************************************************************
!***********************************************************************



      subroutine LINMIN(CCTK_ARGUMENTS,PP,XI,N,FRET,TOL)

!     Given an N dimensional point PP and an N dimensional direction XI,
!     moves and resets PP to where the function FUNC takes on a
!     minimum along the direction XI from PP, and replaces XI by the
!     actual vector displacement that PP was moved.  Also returns as
!     FRET the value of FUNC at the returned location PP.  This is
!     actually all acomplished by calling the routines MNBRAK and
!     PARABOLA.

      use F1COM

      implicit none

      DECLARE_CCTK_ARGUMENTS

      integer i,j,N
      integer l,m

      CCTK_REAL XA,XB,XX,FA,FB,FX
      CCTK_REAL FRET,PARABOLA,TOL
      CCTK_REAL zero
      CCTK_REAL dx,dy,dz

      CCTK_REAL, dimension(1:N) :: PP,XI

      CCTK_REAL F1DIM

      EXTERNAL  F1DIM


!     **********************************************
!     ***   ALLOCATE STORAGE FOR SHARED ARRAYS   ***
!     **********************************************

      zero = 0.0d0

      allocate(PCOM(1:N),XICOM(1:N))


!     *************************
!     ***   START ROUTINE   ***
!     *************************

      dx = cctk_delta_space(1)
      dy = cctk_delta_space(2)
      dz = cctk_delta_space(3)

      NCOM = N

      do j=1,N
         PCOM(j)  = PP(j)
         XICOM(j) = XI(j)
      end do

!     Initial guess for brackets.

      XA = 0.0D0
      XX = 0.1D0*min(dx,dy,dz)

!     Bracket the minimum.

      call MNBRAK(CCTK_ARGUMENTS,XA,XX,XB,FA,FX,FB,F1DIM)

!     Check if function is constant.

      if ((FA.ne.FX).or.(FB.ne.FX)) then

!        Call 1D minimization.

         FRET = PARABOLA(CCTK_ARGUMENTS,XA,XX,XB,FA,FX,FB,F1DIM,TOL)

!        Construct the vector results to return.

         XI = XMIN*XI
         PP  = PP + XI

      end if

      
!     ******************************
!     ***   DEALLOCATE STORAGE   ***
!     ******************************

      deallocate(PCOM,XICOM)


!     ***************
!     ***   END   ***
!     ***************

      end subroutine LINMIN




!***********************************************************************
!***********************************************************************




      CCTK_REAL function F1DIM(CCTK_ARGUMENTS,XX)

!     Must accompany LINMIN

      use F1COM

      implicit none

      DECLARE_CCTK_ARGUMENTS

      integer j

      CCTK_REAL XX,FUNC

      CCTK_REAL, dimension(1:NCOM) :: XT


!     *************************
!     ***   START ROUTINE   ***
!     *************************

      do j=1,NCOM
         XT(j) = PCOM(j) + XX*XICOM(j)
      end do

      F1DIM = FUNC(CCTK_ARGUMENTS,XT,NCOM)


!     ***************
!     ***   END   ***
!     ***************

      end function F1DIM




!***********************************************************************
!***********************************************************************



      subroutine MNBRAK(CCTK_ARGUMENTS,XA,XB,XC,FA,FB,FC,FUNC)

!     Given a function FUNC, and given distinct initial points
!     XA and XB, this routine searches in the downhill direction
!     (defined by the function as evaluated at the initial points)
!     and returns new points XA, XB, XC which bracket a minimum
!     of the function.  Also returned are the function values
!     at the three points, FA, FB and FC.

      implicit none

      DECLARE_CCTK_ARGUMENTS

      CCTK_REAL XA,XB,XC
      CCTK_REAL FA,FB,FC,FU,FUNC
      CCTK_REAL DUM,Q,SS,U,ULIM
      CCTK_REAL GOLD,GLIMIT,TINY

      parameter (GOLD=1.618034D0,GLIMIT=4.0D0,TINY=1.0D-20)


!     *************************
!     ***   START ROUTINE   ***
!     *************************

      FA = FUNC(CCTK_ARGUMENTS,XA)
      FB = FUNC(CCTK_ARGUMENTS,XB)

!     Switch roles of A and B so that we can go downhill
!     in the direction from A to B.

      if (FB.gt.FA) then

         DUM = XA
         XA  = XB
         XB  = DUM

         DUM = FB
         FB  = FA
         FA  = DUM

      end if

!     First guess for C.

      XC = XB + GOLD*(XB-XA)

      FC = FUNC(CCTK_ARGUMENTS,XC)

!     Check if the function is constant.

      if ((FA.eq.FB).and.(FB.eq.FC)) then
         return
      end if

!     Keep returning here until we bracket.

      do while(FB.gt.FC)

!        Compute U by parabollic extrapolation from A, B, C.
!        TINY is used to prevent any possible division by zero.

         SS = (XB-XA)*(FB-FC)
         Q = (XB-XC)*(FB-FA)
         U = XB - ((XB-XC)*Q - (XB-XA)*SS) &
              / (2.0D0*SIGN(MAX(ABS(Q-SS),TINY),Q-SS))

!        We will not go further than this.

         ULIM = XB + GLIMIT*(XC-XB)

!        Now to test various possibilities.

         if ((XB-U)*(U-XC).gt.0.0D0) then

!           Parabolic U is between B and C, try it.

            FU = FUNC(CCTK_ARGUMENTS,U)

            if (FU.lt.FC) then

!              Got a minimum between B and C.

               XA = XB
               FA = FB
               XB = U
               FB = FU
               return

            else if (FU.gt.FB) then

!              Got a minimum between A and U.

               XC = U
               FC = FU
               return

            end if

!           Parabolic fit was no use.  Use default magnification.

            U = XC + GOLD*(XC-XB)

            FU = FUNC(CCTK_ARGUMENTS,U)

         else if ((XC-U)*(U-ULIM).gt.0.0D0) then

!           Parabolic fit is between C and its allowed limit.

            FU = FUNC(CCTK_ARGUMENTS,U)

            if (FU.lt.FC) then

               XB = XC
               XC = U
               U  = XC + GOLD*(XC-XB)
               FB = FC
               FC = FU

               FU = FUNC(CCTK_ARGUMENTS,U)

            end if

         else if ((U-ULIM)*(ULIM-XC).ge.0.0D0) then

!           Limit parabolic U to maximum allowed value.

            U = ULIM

            FU = FUNC(CCTK_ARGUMENTS,U)

         else

!           Reject parabolic U, use default magnification.

            U = XC + GOLD*(XC-XB)
            FU = FUNC(CCTK_ARGUMENTS,U)

         end if

!        Eliminate oldest point and continue.

         XA = XB
         XB = XC
         XC = U
         FA = FB
         FB = FC
         FC = FU

      end do


!     ***************
!     ***   END   ***
!     ***************

      end subroutine MNBRAK




!***********************************************************************
!***********************************************************************




      CCTK_REAL function PARABOLA(CCTK_ARGUMENTS,XA,XB,XC,FXA,FXB,FXC, &
           FUNC,TOL)

!     Given a function F, and given a bracketing triplet of abscissas
!     XA, XB, XC (such that XB is between XA and XC, and F(XB) is less
!     than both F(XA) and F(XC)), this routine isolates the minimum
!     to a fractional precision of about TOL using inverse parabolic
!     interpolation. The abscissa of the minimum is returned as XMIN,
!     and the minimum function value is returned as PARABOLA, the
!     returned function value.

      use F1COM

      implicit none

      DECLARE_CCTK_ARGUMENTS

      integer ITER,ITMAX

      CCTK_REAL XA,XB,XC,FXA,FXB,FXC,FUNC,TOL
      CCTK_REAL A,B,D,PP,Q,SS,U,XX,XM,TOL1
      CCTK_REAL FA,FB,FU,FX,FP,ZEPS

      parameter (ITMAX=100,ZEPS=1.0D-10)


!     *************************
!     ***   START ROUTINE   ***
!     *************************

      error = .false.

!     A and B must be in ascending order, though the input
!     abscissas need not be.

      if (XC.lt.XA) then
         A = XC
         B = XA
         FA = FXC
         FB = FXA
      else
         A = XA
         B = XC
         FA = FXA
         FB = FXC
      end if

      XX = XB
      FX = FXB
      FP = FA

!     Main program loop.

      do ITER=1,ITMAX

         XM = 0.5D0*(A + B)

!        Test for done here.

         TOL1 = 0.5D0*TOL*(abs(FP)+abs(FX)) + ZEPS
         if (abs(FP-FX).le.TOL1) goto 3

!        Construct a trial parabolic fit.

         SS = (XX-A)*(FX-FB)
         Q = (XX-B)*(FX-FA)
         PP = (XX-B)*Q - (XX-A)*SS
         Q = 2.0D0*(Q-SS)

!        Check if Q is zero.  This can only happen if the
!        three points are colinear (in which case they where
!        not bracketing a minimum in the first place), or
!        if two of them are in fact the same point.

         if (Q.eq.0.0D0) then
            write(*,*) 'AHFinder: problem in PARABOLA'
            error = .true.
            return            
         end if

!        Find new point.

         D = - PP/Q
         U = XX + D

!        Check if X was already the minimum.

         if (abs(U-XX).lt.ZEPS) goto 3

!        This is the one function evaluation per iteration.

         FU = FUNC(CCTK_ARGUMENTS,U)

!        And now we have to decide what to do with our function
!        evaluation.

         if (FU.le.FX) then

            if (U.ge.XX) then
               A = XX
               FA = FX
            else
               B = XX
               FB = FX
            end if

            XX = U

            FP = FX
            FX = FU

         else

            if (U.lt.XX) then
               A = U
               FA = FU
            else
               B = U
               FB = FU
            end if

         end if

      end do

      write(*,*) 'PARABOLA exceeded maximum number of iterations.'

    3 XMIN = XX
      PARABOLA = FX


!     ***************
!     ***   END   ***
!     ***************

      end function PARABOLA




!***********************************************************************
!***********************************************************************



      CCTK_REAL function FUNC(CCTK_ARGUMENTS,PP,N)

!     This is a shell to find the surface integral of the expansion.

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      logical flag

      integer i,l,m,N
      integer iter

      CCTK_REAL zero,half,one
      CCTK_REAL aux

      CCTK_REAL, dimension(1:N) :: PP

      character(len=200) :: logf

      save iter,logf


!     **************************
!     ***   DEFINE NUMBERS   ***
!     **************************

      zero = 0.0D0
      half = 0.5D0
      one  = 1.0D0


!     ************************
!     ***   OPEN LOGFILE   ***
!     ************************

      logf = filestr(1:nfile)//"/ahf_logfile"

      if (logfile.and.(myproc.eq.0)) then
         open(12,file=logf,form='formatted',status='old', &
              position='append')
      end if


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


!     *********************************************
!     ***   CALCULATE EXPANSION IF WE NEED TO   ***
!     *********************************************

!     Firstcall?

      if (firstfun) then

         iter = 1

!        Find horizon function and expansion.

         call AHFinder_fun(CCTK_ARGUMENTS)
         call AHFinder_exp(CCTK_ARGUMENTS)

!        Save old values of coefficients.

         xc_old = xc
         yc_old = yc
         zc_old = zc

         c0_old = c0
         cc_old = cc
         cs_old = cs

         firstfun = .false.

!     Not first call.

      else

         iter = iter + 1

!        Check if coefficients changed from last call.
!        I check this because the minimization algorithm
!        sometimes reevaluates the function at the same
!        point twice.  Also, if only c0(0) changed, the
!        expansion will not be affected.

         flag = .false.

         if (wander) then
            if (xc.ne.xc_old) flag = .true.
            if (yc.ne.yc_old) flag = .true.
            if (zc.ne.zc_old) flag = .true.
         end if

         do l=1+stepz,lmax,1+stepz
            if (c0(l).ne.c0_old(l)) flag = .true.
         end do

         if (nonaxi) then
            do l=1,lmax
               do m=1+stepx,l,1+stepx
                  if (stepz*mod(l-m,2).eq.0) then
                     if (cc(l,m).ne.cc_old(l,m)) flag = .true.
                  end if
               end do
               if (.not.refy) then
                  do m=1,l,1+stepx
                     if (stepz*mod(l-m,2).eq.0) then
                        if (cs(l,m).ne.cs_old(l,m)) flag = .true.
                     end if
                  end do
               end if
            end do
         end if

!        If necessary, find horizon function and expansion.

         if (flag) then
            call AHFinder_fun(CCTK_ARGUMENTS)
            call AHFinder_exp(CCTK_ARGUMENTS)
         end if

!        Save old values of coefficients.

         c0_old = c0
         cc_old = cc
         cs_old = cs

      end if

      if (myproc.eq.0) then

         if (veryver) then
            write(*,*)
            write(*,"(A13,I5)") 'FUNC CALL  =',iter
         end if

         if (logfile) then
            write(12,*)
            write(12,"(A13,I5)") 'FUNC CALL  =',iter
         end if

      end if


!     **********************************
!     ***   FIND SURFACE INTEGRALS   ***
!     **********************************

      call AHFinder_int(CCTK_ARGUMENTS)

!     Messages to screen and logfile.

      if (myproc.eq.0) then

!        Find inverse of area.  Be carefull not to divide by zero!
!        Also, a value of 1.0D10 indicates an error in the integration
!        routine, and all integrals end up with the same value.
!        So dont divide by it or I will be confused later.

         if ((intarea.ne.zero).and.(intarea.ne.1.0D10)) then
            aux = one/intarea
         else
            aux = one
         end if

!        Screen.

         if (veryver) then
            write(*,*)
            write(*,"(A21,ES14.6)") ' Surface area       =',intarea
            write(*,"(A21,ES14.6)") ' Mean value of H    =',aux*intexp
            write(*,"(A21,ES14.6)") ' Mean value of H^2  =',aux*intexp2
            write(*,*) 'Number of interpolated points:     ', &
                 inside_min_count
            write(*,*) 'Number of those that are negative: ', &
                 inside_min_neg_count
         end if

!        Logfile.

         if (logfile) then
            write(12,*)
            write(12,"(A21,ES14.6)") ' Surface area       =',intarea
            write(12,"(A21,ES14.6)") ' Mean value of H    =',aux*intexp
            write(12,"(A21,ES14.6)") ' Mean value of H^2  =',aux*intexp2
            write(12,*) 'Number of interpolated points:     ', &
                 inside_min_count
            write(12,*) 'Number of those that are negative: ', &
                 inside_min_neg_count
         end if

      end if


!     ******************************
!     ***   WRITE COEFFICIENTS   ***
!     ******************************

      if (myproc.eq.0) then

!        Headers to screen.

         if (veryver) then

            write(*,*)
            write(*,"(A20)") ' Shape coefficients:'

            if (offset) then
               write(*,*)
               write(*,"(A8,ES14.6)") '   xc  =',xc
               write(*,"(A8,ES14.6)") '   yc  =',yc
               write(*,"(A8,ES14.6)") '   zc  =',zc
            end if

            write(*,*)
            write(*,*) '  l            c0_l'
            write(*,*)
            write(*,"(I4,A6,ES14.6)") 0,'      ',c0(0)

         end if

!        Headers to logfile.

         if (logfile) then

            write(12,*)
            write(12,"(A20)") ' Shape coefficients:'

            if (offset) then
               write(12,*)
               write(12,"(A8,ES14.6)") '   xc  =',xc
               write(12,"(A8,ES14.6)") '   yc  =',yc
               write(12,"(A8,ES14.6)") '   zc  =',zc
            end if

            write(12,*)
            write(12,*) '  l            c0_l'
            write(12,*)
            write(12,"(I4,A6,ES14.6)") 0,'      ',c0(0)

         end if

!        c0 coefficients.

         do l=1+stepz,lmax,1+stepz
            if (veryver) then
               write(*,"(I4,A6,ES14.6)") l,'      ',c0(l)
            end if
            if (logfile) then
               write(12,"(I4,A6,ES14.6)") l,'      ',c0(l)
            end if
         end do

!        cc and cs coefficients (only for nonaxisymmetric case).

         if (nonaxi) then

!           Reflection symmetry on y:  some coefficients are not there.

            if (.not.refy) then

               if (veryver) then
                  write(*,*)
                  write(*,*) '  l   m       cc_lm          cs_lm'
                  write(*,*)
               end if

               if (logfile) then
                  write(12,*)
                  write(12,*) '  l   m       cc_lm          cs_lm'
                  write(12,*)
               end if

               do l=1,lmax
                  do m=1,l
                     if (stepz*mod(l-m,2).eq.0) then
                        if (veryver) then
                           write(*,10) l,m,'  ',cc(l,m),' ',cs(l,m)
                        end if
                        if (logfile) then
                           write(12,10) l,m,'  ',cc(l,m),' ',cs(l,m)
                        end if
   10                   format(I4,I4,A2,ES14.6,A1,ES14.6)
                     end if
                  end do
               end do

!           No reflection symmetry on y.

            else

               if (veryver) then
                  write(*,*)
                  write(*,*) '  l   m       cc_lm'
                  write(*,*)
               end if

               if (logfile) then
                  write(12,*)
                  write(12,*) '  l   m       cc_lm'
                  write(12,*)
               end if

               do l=1,lmax
                  do m=1+stepx,l,1+stepx
                     if (stepz*mod(l-m,2).eq.0) then
                        if (veryver) then
                           write(*,20) l,m,'  ',cc(l,m)
                        end if
                        if (logfile) then
                           write(12,20) l,m,'  ',cc(l,m)
                        end if
   20                   format(I4,I4,A2,ES14.6)
                     end if
                  end do
               end do

            end if
         end if
      end if


!     *********************
!     ***   FIND FUNC   ***
!     *********************

!     The value of FUNC will be equal to the surface integral
!     of the square of the expansion or to the surface area,
!     depending on what we want to minimize.

      if (minarea) then
         FUNC = intarea
      else
         if(find_trapped_surface) then
            FUNC = intexpdel2
         else
            FUNC = intexp2
         endif
      end if


!     *************************
!     ***   CLOSE LOGFILE   ***
!     *************************

      if (logfile.and.(myproc.eq.0)) then
         close(12)
      end if


!     ***************
!     ***   END   ***
!     ***************

      end function FUNC
