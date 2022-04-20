/*@@
  @file      AHFinder_dat.F
  @date      April 1998
  @author    Miguel Alcubierre
  @desc 
             F90 module containing data to be shared by
             the different subroutines.

             For the meaning of the spectral components of
             the different flows see:

             C. Gundlach, Phys. Rev. D 57, 863 (1998).

  @enddesc 
@@*/

#include "cctk.h"

      module AHFinder_dat

      implicit none

      logical firstfun,firstleg,firstint
      logical nonaxi,refx,refy,refz,cartoon
      logical offset,wander
      logical minarea,flow,find_trapped_surface
      logical logfile,verbose,veryver,guessverbose
      logical interror,guessold
      logical find3,sloppy,inner,guess_absmin,manual_guess
      logical status_old,status_old_0,status_old_1,status_old_2

      integer lmax,stepx,stepy,stepz
      integer ntheta,nphi
      integer mfind,ahf_ncall
      integer nx,ny,nz
      integer nprocs,myproc,ierr

      CCTK_INT nfile
      CCTK_INT interror1,interror2,interror3
      CCTK_INT inside_min_count,inside_min_neg_count

      CCTK_REAL xc,yc,zc
      CCTK_REAL xmn,ymn,zmn,xmx,ymx,zmx
      CCTK_REAL xc_old,yc_old,zc_old
      CCTK_REAL intexp,intexp2,intexpdel2,intarea
      CCTK_REAL rhmax,rhmin
      CCTK_REAL circ_eq,meri_p1,meri_p2

      CCTK_REAL hw,cw,nw
      CCTK_REAL dx,dy,dz

      CCTK_REAL avgx,avgy,avgz

      CCTK_REAL, allocatable, dimension (:)   :: c0,c0_0,c0_1,c0_2,c0_old
      CCTK_REAL, allocatable, dimension (:,:) :: cc,cc_0,cc_1,cc_2,cc_old
      CCTK_REAL, allocatable, dimension (:,:) :: cs,cs_0,cs_1,cs_2,cs_old

      CCTK_REAL, allocatable, dimension (:)   :: hflow0,cflow0,nflow0
      CCTK_REAL, allocatable, dimension (:,:) :: hflowc,cflowc,nflowc
      CCTK_REAL, allocatable, dimension (:,:) :: hflows,cflows,nflows

      character(len=200) :: filestr

      data ahf_ncall / 0 /
      data status_old   / .false. /
      data status_old_0 / .false. /
      data status_old_1 / .false. /
      data status_old_2 / .false. /

      save ahf_ncall
      save c0,c0_0,c0_1,c0_2
      save cc,cc_0,cc_1,cc_2
      save cs,cs_0,cs_1,cs_2
      save status_old,status_old_0,status_old_1,status_old_2

!     Description of variables:
!
!     firstfun   First call to function FUNC?
!     firstleg   First call to function LEGEN?
!     firstint   First call to subroutine AHFinder_int?
!
!     nonaxi     Is the data non-axisymmetric?
!
!     refx       Reflection symmetry x->-x?
!
!     refy       Reflection symmetry y->-y?
!
!     refz       Reflection symmetry z->-z?
!
!     cartoon    Are we using cartoon_2d?
!
!     inner      Look for inner horizon?
!
!     sloppy     Sloppy initial guess?
!
!     offset     Is the center offset from the origin?
!     wander     Do we allow the center to wander?
!
!     minarea    Do we want to minimize area?
!
!     flow       Use fast flow algorithm instead of minimization?
!
!     logfile    Write log file?
!
!     verbose    Write messages to screen?
!
!     veryver    Write lots of messages to screen?
!
!     guessverbose   Write info on initial guess points?
!
!     lmax       Maximum number of terms in theta expansion.
!
!     stepx      1 if (refx=1), 0 otherwise.
!     stepy      1 if (refy=1), 0 otherwise.
!     stepz      1 if (refz=1), 0 otherwise.
!
!     ntheta     Number of subdivisions in theta.
!     nphi       Number of subdivisions in phi.
!
!     xc         Coordinate x of centre of surface.
!     yc         Coordinate y of centre of surface.
!     zc         Coordinate z of centre of surface.
!
!     xmn        Minimum value of x over the grid.
!     xmx        Maximum value of x over the grid.
!
!     ymn        Minimum value of y over the grid.
!     ymx        Maximum value of y over the grid.
!
!     zmn        Minimum value of z over the grid.
!     zmx        Maximum value of z over the grid.
!
!     intexp     Integral of expansion over the surface.
!     intexp2    Integral of expansion squared over the surface.
!     intarea    Area of the surface.
!
!     rhmax      Maximum radius of horizon.
!     rhmin      Minimum radius of horizon.
!
!     circ_eq    Equatorial circumference.
!     meri_p1    Length of meridian at phi=0.
!     meri_p2    Length of meridian at phi=pi/2.
!
!     c0(l)      Coefficients of LEGEN(l,0,cos(theta))
!     cc(l,m)    Coefficients of LEGEN(l,m,cos(theta))*cos(m*phi)
!     cs(l,m)    Coefficients of LEGEN(l,m,cos(theta))*sin(m*phi)
!
!     hflow*     Spectral components for "H" flow.
!     cflow*     Spectral components for "C" flow.
!     nflow*     Spectral components for "N" flow.
!
!     interror1  Different from zero if radius is negative.
!     interror2  Different from zero if outside computational domain.
!     interror3  Different from zero if inside mask
!
!     interror   True if any of the 3 errors above are non-zero.
!
!     find3      Look for 3 horizons?
!     mfind      Counter to see for which horizon we are looking.
!
!     inside_min_count        Number of elements in integral (should always
!                             be ntheta*nphi)
!
!     inside_min_neg_count    Number of elements in integral that have
!                             negative expansion.
!
!     avgx        X-position of centroid of apparent horizon
!     avgy        Y-position of centroid of apparent horizon
!     avgz        Z-position of centroid of apparent horizon


!     ***************
!     ***   END   ***
!     ***************

      end module AHFinder_dat


