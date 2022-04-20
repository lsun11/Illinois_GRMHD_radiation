/*@@
  @file      AHFinder_leg.F
  @date      April 1998
  @author    Miguel Alcubierre
  @desc 
             Associated Legendre functions.
  @enddesc 
@@*/

#include "cctk.h"

      CCTK_REAL function LEGEN(L,M,X)

!     This function computes the associated Legendre polynomial
!     P^M_L(X).  Here M and L are integers satisfying 0 <= M <= L,
!     and X lies in the range -1 <= X <= 1.
!
!     The general function is taken from Numerical Recipes.
!     However, I hardwired the first few polynomials for M=0,
!     and I added a section that tries to use as much information
!     as possible from previous calls to the routine to avoid
!     having to start from scratch every time.

      use AHFinder_dat

      implicit none

      integer L,M
      integer LOLD1,LOLD2,MOLD1,MOLD2
      integer I,LL

      CCTK_REAL X,Y
      CCTK_REAL PMM,PLL,PMMP1,SOMX2,FACT
      CCTK_REAL factor,aux
      CCTK_REAL zero,one,two
      CCTK_REAL XOLD1,XOLD2
      CCTK_REAL LEGOLD1,LEGOLD2

      save LOLD1,LOLD2,MOLD1,MOLD2,XOLD1,XOLD2
      save LEGOLD1,LEGOLD2


!     **************************
!     ***   DEFINE NUMBERS   ***
!     **************************

      zero = 0.0D0
      one  = 1.0D0
      two  = 2.0D0


!     *************************************************
!     ***   HARDWIRE FIRST FEW POLYNOMIALS FOR M=0  ***
!     *************************************************

!     Notice that it does not make sense to go much higher
!     than L=12 since the explicit expressions for the
!     polynomials have large terms with alternating signs,
!     so for high L we rapidly loose accuracy in the sum.
!     The method from Numerical Recipes that comes below
!     is slow, but does not have this problem.

      Y = X**2

!     Find polynomials for M=0.

      if ((M.eq.0).and.(L.le.12)) then

         if (L.eq.0) then
            LEGEN = one
            goto 10
         else if (L.eq.1) then
            LEGEN = X
            goto 10
         else if (L.eq.2) then
            LEGEN = 5.0D-1*(3.0D0*Y - one)
            goto 10
         else if (L.eq.3) then
            LEGEN = 5.0D-1*(5.0D0*Y - 3.0D0)*X
            goto 10
         else if (L.eq.4) then
            LEGEN = 1.25D-1*(35.0D0*Y**2 - 30.0D0*Y + 3.0D0)
            goto 10
         else if (L.eq.5) then
            LEGEN = 1.25D-1*(63.0D0*Y**2 - 70.0D0*Y + 15.0D0)*X
            goto 10
         else if (L.eq.6) then
            LEGEN = 6.25D-2*(231.0D0*Y**3 - 315.0D0*Y**2 + 105.0D0*Y &
                 - 5.0D0)
            goto 10
         else if (L.eq.7) then
            LEGEN = 6.25D-2*(429.0D0*Y**3 - 693.0D0*Y**2 + 315.0D0*Y &
                 - 35.0D0)*X
            goto 10
         else if (L.eq.8) then
            LEGEN = 7.8125D-3*(6435.0D0*Y**4 - 12012.0D0*Y**3 &
                 + 6930.0D0*Y**2 - 1260.0D0*Y + 35.0D0)
            goto 10
         else if (L.eq.9) then
            LEGEN = 7.8125D-3*(12155.0D0*Y**4 - 25740.0D0*Y**3 &
                 + 18018.0D0*Y**2 - 4620.0D0*Y + 315.0D0)*X
            goto 10
         else if (L.eq.10) then
            LEGEN = 3.90625D-3*(46189.0D0*Y**5 - 109395.0D0*Y**4 &
                 + 90090.0D0*Y**3 - 30030.0D0*Y**2 + 3465.0D0*Y &
                 - 63.0D0)
            goto 10
         else if (L.eq.11) then
            LEGEN = 3.90625D-3*(88179.0D0*Y**5 - 230945.0D0*Y**4 &
                 + 218790.0D0*Y**3  - 90090.0D0*Y**2 + 15015.0D0*Y &
                 - 693.0D0)*X
            goto 10
         else if (L.eq.12) then
            LEGEN = 9.76563D-4*(676039.0D0*Y**6 - 1939938.0D0*Y**5 &
                 + 2078505.0D0*Y**4 - 1021020.0D0*Y**3 &
                 + 225225.0D0*Y**2 - 18018.0D0*Y + 231.0D0)
            goto 10
         end if

      end if


!     ***********************************************
!     ***   USE RECURSIVE RELATIONS IF POSSIBLE   ***
!     ***********************************************

!     Here I try to optimize the calculation of the Legendre
!     polynomials (which is ordinarily very expensive) by using
!     the recurrence relations and information from the last
!     calls to the routine.

      if (firstleg) then

         LOLD1   = -2
         MOLD1   = -2
         XOLD1   = zero
         LEGOLD1 = zero

         firstleg = .false.

      else

!        See if we have the same values of M and X as the last call.

         if ((M.eq.MOLD1).and.(X.eq.XOLD1).and.(L.ne.M)) then

!           If L is the same as in the last call, then
!           LEGEN should also be the same. This will
!           almost certainly never happen.

            if (L.eq.LOLD1) then

               LEGEN = LEGOLD1
               goto 10

!           This is the most interesting case L = LOLD1+1.

            else if (L.eq.LOLD1+1) then

!              LOLD1 = MOLD1.

               if (LOLD1.eq.MOLD1) then

                  LEGEN = X*dble(2*M+1)*LEGOLD1
                  goto 10

!              LOLD1 = LOLD2 + 1.

               else if ((L.eq.LOLD2+2).and.(X.eq.XOLD2)) then

                  LEGEN = (X*dble(2*L-1)*LEGOLD1 - dble(L+M-1)*LEGOLD2) &
                       / dble(L-M)
                  goto 10

               end if

            end if

         end if

      end if


!     *************************************
!     ***   NUMERICAL RECIPES ROUTINE   ***
!     *************************************

!     Compute P^M_M.

      PMM = one

      if (M.gt.0) then

         SOMX2 = sqrt((one-X)*(one+X))
         FACT = one

         do I=1,M
            PMM = - PMM*FACT*SOMX2
            FACT = FACT + two
         end do

      end if

      if (L.eq.M) then

         LEGEN = PMM

!     Compute P^M_M+1.

      else

         PMMP1 = X*dble(2*M+1)*PMM

         if (L.eq.M+1) then

            LEGEN = PMMP1

!        Compute P^M_L, L > M+1.

         else

            do LL=M+2,L
               PLL = (X*dble(2*LL-1)*PMMP1 - dble(LL+M-1)*PMM) &
                    / dble(LL-M)
               PMM = PMMP1
               PMMP1 = PLL
            end do

            LEGEN = PLL

          end if

      end if


!     ***************************
!     ***   SAVE OLD VALUES   ***
!     ***************************

   10 continue

      LOLD2   = LOLD1
      MOLD2   = MOLD1
      XOLD2   = XOLD1
      LEGOLD2 = LEGOLD1

      LOLD1   = L
      MOLD1   = M
      XOLD1   = X
      LEGOLD1 = LEGEN


!     **********************************************
!     ***   MULTIPLY WITH NORMALIZATION FACTOR   ***
!     **********************************************

!     For M=0 I use the normalization factor:
!
!     sqrt(2*L + 1)
!
!     For M non-zero I use the normalization factor:
!
!     sqrt( 2 (2*L + 1) (L-M)! / (L+M)! )
!
!     The extra factor of sqrt(2) is there because I
!     use a basis of sines and cosines and not the
!     standard complex exponential.
!
!     With this normalization, I ensure that my basis
!     functions f _lm are such that:
!
!     /
!     | f_lm f_l;mp sin(theta) dtheta dphi  =  4 pi delta_mm' delta_ll'
!     /
!
!     Notice the extra factor of 4 pi.  This is there because
!     for a sphere, I want the (0,0) coefficient to correspond
!     to the radius.

      factor = dble(2*L+1)

      if (M.ne.0) then
 
         aux = one
 
         do i=L-M+1,L+M
            aux = aux*dble(i)
         end do
 
         factor = two*factor/aux
 
      end if

      LEGEN = sqrt(factor)*LEGEN


!     ***************
!     ***   END   ***
!     ***************

      end function LEGEN



