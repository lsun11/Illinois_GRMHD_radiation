/*@@
  @file      AHFinder_calcsigma.F
  @date      July 1999
  @author    Lars Nerger
  @desc 
             Calculate weight sigma for Nakamura flow
  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      subroutine AHFinder_calcsigma(CCTK_ARGUMENTS,xp,yp,zp,gupij,sigma)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer nmax
      integer l,m
      integer indx1,indx2
      integer indx,idir,jdir
      
      CCTK_REAL zero, sigma
      CCTK_REAL ni(3),si(3),supi(3),pij(3,3)
      CCTK_REAL ylmi((lmax+1)**2,3)
      CCTK_REAL alm((lmax+1)**2)
      CCTK_REAL xp,yp,zp,er,sum
      CCTK_REAL rhoradius,costheta,sintheta,sinthetal(-1:lmax),prod
      CCTK_REAL zlm0(0:lmax,0:lmax),zlm1(0:lmax,0:lmax)
      CCTK_REAL cosphi,sinphi,cosmphiangle(0:lmax),sinmphiangle(0:lmax)
      CCTK_REAL thetax,thetay,thetaz
      CCTK_REAL phianglex,phiangley
      CCTK_REAL c1,c2,c3,metnorm
      CCTK_REAL fi(3)
      CCTK_REAL gupij(3,3)

!     Description of variables
!
!     sigma   Weight of Nakamura Flow


!     **************************
!     ***   DEFINE NUMBERS   ***
!     **************************

      zero    = 0.0D0
      nmax = (lmax+1)**2

!     ******************************************************
!     ***   Calculate gradient of Y_lm(x,y,z)            ***
!     ***   Adopted from thorn_SpectralAHF by Gundlach  ***
!     ******************************************************
      

!     Initialize arrays
      ylmi = zero
      zlm0 = zero
      zlm1 = zero
      ni = zero
      si = zero
      supi = zero
      pij = zero
      alm = zero
      sinthetal = zero
      cosmphiangle = zero
      sinmphiangle = zero


!     Calculate quantities in spherical coordinates.
      er = sqrt(xp**2 + yp**2 + zp**2)
      if (er .eq. 0.d0) then
         call CCTK_WARN(0,'r=0 in AHFinder_calcsigma')
      end if
     
C     ********PATCH******************************************
C     Keep away from the z-axis, artificially.
      if (xp**2 + yp**2 .lt. 1.d-16) then
C         write(6,*) 'patch on axis in ylmder2'
         xp = 0.d0
         yp = 1.d-8
      end if

C     Polar angle theta.
      rhoradius = sqrt(xp**2 + yp**2)
      costheta = zp / er
      sintheta = rhoradius / er
C     Array sinthetal(l) contains [sin(theta)]^l
      sinthetal(-1) = 0.d0
      sinthetal(0) = 1.d0
      sinthetal(1) = sintheta
      do l=2,lmax
         sinthetal(l) = sinthetal(l-1) * sinthetal(1)
      end do

C     zlm0 is Y_lm stripped of exp(i m phiangle) and multiplied by sqrt(4pi).
C     zlm1 and zlm2 are the first and second derivative of zlm0 
C     with respect to theta. l=0 and m=l are treated separately.
C     We only need zlm0,1,2 for positive m.
C     Treat l=0 separately:
      zlm0(0,0) = 1.d0
      zlm1(0,0) = 0.d0

C     Other l>0.
      prod = 1.d0
      do l=1,lmax

C     Initialize m = l:
C     c1 = sqrt(2l+1) * sqrt(2l!) / 2^l / l! 
C     Y_ll = c1 * sin(theta)^l 
C     Y_llp = c1 * l * sin(theta)^(l-1) * costheta
         m = l
         prod = - prod * sqrt(dble(2*l-1) / dble(2*l)) 
         c1 = sqrt(dble(2*l+1)) * prod
         zlm0(l,m) = c1 * sinthetal(l)
         zlm1(l,m) = c1 * dble(l) * sinthetal(l-1) * costheta

C     Recursion relations for the other m:
C     c2 = sqrt[(2l+1)(2l-1)/(l^2-m^2)]
C     c3 = sqrt[((l-1)^2-m^2)(2l=1)/(l^2-m^2)/(2l-3)]
C     Y_lm = c2 Y_l-1,m cos(theta) + c3 Y_l-2,m
C     Y_lmp = c2 Y_l-1,mp cos(theta) + c2 Y_l-1,m + c3 Y_l-2,mp
         m = l-1
         c2 = sqrt(dble((2*l+1)*(2*l-1)) 
     $        / dble(l**2-m**2)) 
         zlm0(l,m) = c2 * zlm0(l-1,m) * costheta
         zlm1(l,m) = c2 * (zlm1(l-1,m) * costheta 
     $           - zlm0(l-1,m) * sintheta) 

         do m=0,l-2
            c2 = sqrt(dble((2*l+1)*(2*l-1)) 
     $           / dble(l**2-m**2)) 
            c3 =sqrt(dble((2*l+1)*((l-1)**2-m**2))
     $           / dble((l**2 - m**2)*(2*l-3))) 
            zlm0(l,m) = c2 * zlm0(l-1,m) * costheta 
     $           - c3 * zlm0(l-2,m)
            zlm1(l,m) = c2 * (zlm1(l-1,m) * costheta 
     $           - zlm0(l-1,m) * sintheta)
     $           - c3 * zlm1(l-2,m)
         end do
      end do


      if (rhoradius .ne. 0.d0) then
C     Azimuth angle phiangle is defined.
         cosphi = xp / rhoradius
         sinphi = yp / rhoradius
C     First and second derivatives of theta.
         thetax = costheta * cosphi / er
         thetay = costheta * sinphi / er
         thetaz = - sintheta / er
C     First and second derivatives of phiangle, the angle in the xy-plane.
         phianglex = - sinphi / rhoradius
         phiangley = cosphi / rhoradius
C     Put together arrays of basis functions in the xy plane.
C     The sqrt(2) factor is needed in going from the basis exp(+/-imphiangle)
C     to the basis cos/sin(mphiangle).
         cosmphiangle(0) = 1.d0
         sinmphiangle(0) = 1.d0
         cosmphiangle(1) = sqrt(2.d0) * cosphi
         sinmphiangle(1) = sqrt(2.d0) * sinphi
         do m=2,lmax
            cosmphiangle(m) = cosmphiangle(m-1) * cosphi - sinmphiangle(m-1) * sinphi
            sinmphiangle(m) = cosmphiangle(m-1) * sinphi + sinmphiangle(m-1) * cosphi
         end do
      else
C     See above for a patch that avoids ever getting here.
          call CCTK_WARN(0,'ylmder2 cannot handle axis')
C     For x = y = 0 the following are not defined.
!         cosphi = 7.d77
!         sinphi = 7.d77
!         thetax = 7.d77
!         thetay = 7.d77
!         phianglex = 7.d77
!         phiangley = 7.d77
!         cosmphiangle = 7.d77
!         sinmphiangle = 7.d77
C     Only these are defined.
         thetaz = 0.d0
         cosmphiangle(0) = 1.d0
         sinmphiangle(0) = 1.d0
      end if

C     Assemble the Y_lm from zlm0 and the basis of functions of phiangle.
C     Assemble the Y_lm,i and Y_lm,ij in the corresponding manner.
      do l=0,lmax
         do m=0,l
C     These are the values of indx corresponding to (l,m) and (l,-m).
            indx1 = l**2 + 1 + l + m
            indx2 = l**2 + 1 + l - m
C     Ylm,x
            ylmi(indx1,1) = thetax * zlm1(l,m) * cosmphiangle(m)
     $           - dble(m) * phianglex * zlm0(l,m) * sinmphiangle(m)
            ylmi(indx2,1) = thetax * zlm1(l,m) * sinmphiangle(m) 
     $           + dble(m) * phianglex * zlm0(l,m) * cosmphiangle(m)
C     Ylm,y
            ylmi(indx1,2) = thetay * zlm1(l,m) * cosmphiangle(m)
     $           - dble(m) * phiangley * zlm0(l,m) * sinmphiangle(m)
            ylmi(indx2,2) = thetay * zlm1(l,m) * sinmphiangle(m) 
     $           + dble(m) * phiangley * zlm0(l,m) * cosmphiangle(m)
C     Ylm,z
            ylmi(indx1,3) = thetaz * zlm1(l,m) * cosmphiangle(m)
            ylmi(indx2,3) = thetaz * zlm1(l,m) * sinmphiangle(m)
         end do
      end do

 
!     ********************************
!     ***  First derivatives of r  ***
!     ********************************
    
      ni(1) = xp / er
      ni(2) = yp / er
      ni(3) = zp / er


!     *************************************************************
!     ***   compress c0, cc and cs into alm with single index   ***
!     *************************************************************

      do l=0, lmax
         alm(l**2+l+1) = c0(l)
         do m=1,lmax
            alm(l**2+l+1+m) = cc(l,m)
            alm(l**2+l+1-m) = cs(l,m)
         end do
      end do


!     *******************************************************************
!     ***   Assemble derivatives of f. Recall f = r - h(theta,phiangle):   ***
!     *******************************************************************

C     Derivatives of r...
C     r,i = n_i and r,ij = (delta_ij - n_i n_j) / r = n_ij
      do idir=1,3
         fi(idir) = ni(idir)
      end do
C     ...and derivatives of -h.
      do indx=1,nmax
         do idir=1,3
            fi(idir) = fi(idir) - alm(indx) * ylmi(indx,idir)
         end do
      end do


!     *************************************************************
!     ***   Assemble normal vector s on surfaces of constant f  ***
!     ***   with indices down.  That is, normal with respect    ***
!     ***   to the physical metric g_ij.                        ***
!     *************************************************************

      metnorm = 0.d0
      do idir=1,3
         do jdir=1,3
            metnorm = metnorm + fi(idir) * fi(jdir) 
     $           * gupij(idir,jdir)
         end do
      end do
      metnorm = sqrt(metnorm)
      do idir=1,3
         si(idir) = fi(idir) / metnorm
      end do

C     Now with indices up, raised by g^ij.
      do idir=1,3
         supi(idir) = 0.d0
         do jdir=1,3
            supi(idir) = supi(idir) 
     $           + gupij(idir,jdir) * si(jdir) 
         end do
      end do

!     *******************************************************
!     ***   Assemble projection operator with indices up  ***
!     ***            p^ij = g^ij - s^i s^j                ***
!     *******************************************************

      do idir=1,3
         do jdir=1,3
            pij(idir,jdir) = gupij(idir,jdir)
     $           - supi(idir) * supi(jdir)
         end do
      end do


!     *************************************************
!     ***   Calculate weight sigma                  ***
!     ***   2 * r^2 / (p^ij (delta_ij - n_i n_j))   ***
!     *************************************************

      sum = 0.d0
      do idir=1,3
         sum = sum + pij(idir,idir)
         do jdir=1,3
            sum = sum 
     $           - pij(idir,jdir) * ni(idir) * ni(jdir)
         end do
      end do
      sigma = 2.d0 * er**2 / sum


!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder_calcsigma

