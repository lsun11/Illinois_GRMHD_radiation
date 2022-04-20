/*@@
  @file      AHFinder_fun.F
  @date      April 1998
  @author    Miguel Alcubierre
  @desc 
             Find horizon function.
  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine AHFinder_fun(CCTK_ARGUMENTS)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k
      integer l,m,ll

      CCTK_REAL LEGEN
      CCTK_REAL xp,yp,zp,rp
      CCTK_REAL phiangle,cost,cosa,sina
      CCTK_REAL zero,half,one,two
      CCTK_REAL pi,halfpi,twopi
      CCTK_REAL aux1,aux2
      CCTK_REAL red_tmp
      CCTK_REAL local_legen,local_legen_old,local_legen_old_old 
      CCTK_REAL factor
      CCTK_REAL cc_cosa(lmax,lmax),cs_sina(lmax,lmax)
      CCTK_REAL dsqrt_factor(lmax,lmax),dsqrt_factor2(lmax)
      CCTK_REAL lm_factor1(0:lmax,0:lmax),lm_factor2(0:lmax,0:lmax)
      CCTK_REAL lm_factor3(lmax)
      CCTK_REAL sinphi,cosphi,sinphi1,cosphi1,sinphi2,cosphi2
      CCTK_REAL legen_base


!     **************************
!     ***   DEFINE NUMBERS   ***
!     **************************

      zero = 0.0D0
      half = 0.5D0
      one  = 1.0D0
      two  = 2.0D0

      pi     = acos(-one)
      halfpi = half*pi
      twopi  = two*pi


!     *******************************************
!     ***   FACTORS FOR SPHERICAL HARMONICS   ***
!     *******************************************

!     Compute the normalization factors only once for all (l,m)
!
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
!     | f_lm f_l;mp sin(theta) dtheta dphiangle  =  4 pi delta_mm' delta_ll'
!     /
!
!     Notice the extra factor of 4 pi.  This is there because
!     for a sphere, I want the (0,0) coefficient to correspond
!     to the radius.
!
!     Also compute the multiplying factors in the recursion relations

!     Initialize.

      lm_factor1 = zero
      lm_factor2 = zero
      lm_factor3 = zero

      dsqrt_factor  = zero
      dsqrt_factor2 = zero

!     M=0.

      do l=1,lmax

         dsqrt_factor2(l) = c0(l)*sqrt(dble(2*l+1))

         lm_factor1(l,0) = dble(2*l-1)/dble(l)
         lm_factor2(l,0) = dble(l-1)/dble(l)
         lm_factor3(l)   = dble(2*l-1)

      end do

!     M > 0.

      if (nonaxi) then

         do l=1,lmax

            dsqrt_factor(l,1) = sqrt(2.0D0*dble(2*l+1)/dble(l*(l+1)))

            do m=2,l
               dsqrt_factor(l,m) = dsqrt_factor(l,m-1)
     .                           / sqrt(dble((l-m+1)*(l+m)))
            end do

         end do

         do l=2,lmax

            lm_factor1(l,1) = dble(2*l-1)/dble(l-1)
            lm_factor2(l,1) = dble(l)/dble(l-1)

            do m=2,l-1
               lm_factor1(l,m) = dble(2*l-1)/dble(l-m)
               lm_factor2(l,m) = dble(l+m-1)/dble(l-m)
            end do

         end do

      end if


!     *********************************
!     ***   FIND HORIZON FUNCTION   ***
!     *********************************

!     Loop over j.

      do j=1,ny

!        Find yp.

         yp = y(1,j,1) - yc

!        Loop over i.

         do i=1,nx

!           Find xp.

            xp = x(i,1,1) - xc

!           Find sines and cosines of phiangle.

            if (nonaxi) then

!              Find phiangle.

               if (yp.ne.zero) then
                  phiangle = atan2(yp,xp)
               else
                  phiangle = zero
               end if

!              Compute the expensive sines and cosines only once
!              for all (l,m) using the recursion relations:
!
!              cos(m*phiangle) = 2*cos((m-1)*phiangle)*cos(phiangle)-cos((m-2)*phiangle)
!              sin(m*phiangle) = 2*sin((m-1)*phiangle)*cos(phiangle)-sin((m-2)*phiangle)

               cosphi  = cos(phiangle)
               cosphi1 = one
               cosphi2 = cosphi

               if (.not.refy) then
                  sinphi  = sin(phiangle)
                  sinphi1 = zero
                  sinphi2 = sinphi
               end if

!              M=1.

               do l=1,lmax

                  cc_cosa(l,1) = cc(l,1)*cosphi*dsqrt_factor(l,1)

                  if (.not.refy) then
                     cs_sina(l,1) = cs(l,1)*sinphi*dsqrt_factor(l,1)
                  end if

               end do

!              M>1.

               do m=2,lmax

                  cosa = 2.0D0*cosphi2*cosphi - cosphi1

                  cosphi1 = cosphi2
                  cosphi2 = cosa

                  do l=m,lmax
                     cc_cosa(l,m) = cc(l,m)*cosa*dsqrt_factor(l,m)
                  end do 

                  if (.not.refy) then

                     sina = 2.0D0*sinphi2*cosphi - sinphi1

                     sinphi1 = sinphi2
                     sinphi2 = sina

                     do l=m,lmax
                        cs_sina(l,m) = cs(l,m)*sina*dsqrt_factor(l,m)
                     end do

                  end if

              end do

            end if

!           Loop over k.

            do k=1,nz

!              Find zp.

               zp = z(i,j,k) - zc

!              Find rp.

               rp = sqrt(xp**2 + yp**2 + zp**2)

!              Monopole term.

               aux1 = c0(0)

!              Axisymmetric terms.

               if (rp.ne.zero) then
                  cost = zp/rp
               else
                  cost = one
               end if

!              Compute the contribution from M=0, L=1,Lmax using
!              the recursion relations from Numerical Recipes.

               if (lmax.gt.0) then

!                 L=1

                  legen_base = one
                  local_legen_old_old = one
                  local_legen_old = cost

                  aux1 = aux1 + dsqrt_factor2(1)*local_legen_old

!                 L=2,lmax

                  do l=2,lmax

                     local_legen = cost*lm_factor1(l,0)*local_legen_old
     .                                - lm_factor2(l,0)*local_legen_old_old
                     local_legen_old_old = local_legen_old
                     local_legen_old = local_legen

                     aux1 = aux1 + dsqrt_factor2(l)*local_legen

                  end do
!                  
!                 Non-axisymmetric terms.
!
                  if (nonaxi) then

                     aux2 = sqrt((one-cost)*(one+cost))

!                    This general loop can only be used until l=lmax-2 since
!                    the recursion relations requires 2 starting values.

                     do m=1,lmax-2
! L=M
                        legen_base = -legen_base*aux2*lm_factor3(m)     

                        aux1 = aux1 + legen_base*cc_cosa(m,m)

                        if (.not.refy) then
                           aux1 = aux1 + legen_base*cs_sina(m,m)
                        end if

                        local_legen_old_old = legen_base
! L=M+1
                        local_legen_old = cost*lm_factor3(m+1)         
     .                                  *legen_base

                        aux1 = aux1 + local_legen_old*cc_cosa(m+1,m)

                        if (.not.refy) then
                           aux1 = aux1 + local_legen_old*cs_sina(m+1,m)
                        end if

! L=M+2,lmax
                        do l=m+2,lmax                             

                           local_legen = cost*lm_factor1(l,m)*local_legen_old
     .                                      - lm_factor2(l,m)*local_legen_old_old
                           local_legen_old_old = local_legen_old
                           local_legen_old = local_legen

                           aux1 = aux1 + local_legen*cc_cosa(l,m)

                           if (.not.refy) then
                              aux1 = aux1 + local_legen*cs_sina(l,m)
                           end if

                        end do
                     end do

!                    Then do M=lmax-1

                     if (lmax.gt.1) then

! L=lmax-1
                        m=lmax-1                                   
                        legen_base = -legen_base*aux2*lm_factor3(m)   
                        aux1 = aux1 + legen_base*cc_cosa(m,m)

                        if (.not.refy) then
                           aux1 = aux1 + legen_base*cs_sina(m,m)
                        end if

! L=lmax
                        local_legen = cost*lm_factor3(lmax)*   
     .                                            legen_base
                        aux1 = aux1 + local_legen*cc_cosa(lmax,m)

                        if (.not.refy) then
                           aux1 = aux1 + local_legen*cs_sina(lmax,m)
                        end if

                     end if

!                    And finally do M=lmax.

! L=lmax
                     m=lmax                                         

                     legen_base = -legen_base*aux2*lm_factor3(m)
                     aux1 = aux1 + legen_base*cc_cosa(m,m)

                     if (.not.refy) then
                        aux1 = aux1 + legen_base*cs_sina(m,m)
                     end if

                  end if

               end if

!              Find horizon function.

               ahfgrid(i,j,k) = rp - aux1

            end do
         end do
      end do


!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder_fun

