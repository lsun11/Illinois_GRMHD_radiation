/*@@
  @file      AHFinder_gau.F
  @date      October 1998
  @author    Miguel Alcubierre
  @desc
             Find gaussian curvature of surface.  The gaussian
             curvature is defined as R/2, where R is the Ricci
             scalar of the induced 2-geometry of the surface.

             As an extra "goody", this routine also integrates
             the equatorial circumference of the horizon and
             two polar circumferences at phiangle=0 and phiangle=pi/2.
  @enddesc
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine AHFinder_gau(CCTK_ARGUMENTS)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      logical firstgau
      logical firstcal(4)

      integer i,j,k,l,m,n,pp
      integer npoints,vindex
      integer param_table_handle,interp_handle,coord_system_handle,sum_handle
      integer ierror
      character(len=128) :: operator
      CCTK_INT nchars

      CCTK_INT rerror,error1,error2

      character(30) options_string

      CCTK_REAL LEGEN
      CCTK_REAL theta,phiangle,xp,yp,zp,rp
      CCTK_REAL cost,sint,cosp,sinp
      CCTK_REAL dtheta,dphiangle,dtp,idtheta,idphiangle,phianglestart
      CCTK_REAL trr,ttt,tpp,trt,trp,ttp,ft,fp
      CCTK_REAL deta,ideta
      CCTK_REAL nlx,nly,nlz,nux,nuy,nuz
      CCTK_REAL trxi,xi2
      CCTK_REAL zero,half,one,two,three,four,pi
      CCTK_REAL aux,sina,cosa
c      CCTK_REAL dt,idx,idy,idz
c      CCTK_REAL i2dx,i2dy,i2dz,idxx,idyy,idzz,idxy,idxz,idyz

      CCTK_REAL, dimension(3,3)     :: ug,xi
      CCTK_REAL, dimension(2,2)     :: ga,ua
      CCTK_REAL, dimension(2,2,2)   :: d1a,gammad,gamma
      CCTK_REAL, dimension(2,2,2,2) :: d2a,gamma2

      CCTK_POINTER, dimension(3) :: interp_coords
      CCTK_INT, dimension(7) :: in_array_indices
      CCTK_POINTER, dimension(7) :: out_arrays
      CCTK_INT, dimension(7) :: out_array_type_codes

      CCTK_REAL, allocatable, dimension(:,:) :: rr,xa,ya,za
      CCTK_REAL, allocatable, dimension(:,:) :: txx,tyy,tzz,txy,txz,tyz
      CCTK_REAL, allocatable, dimension(:,:) :: g11,g22,g12

      CCTK_REAL, dimension(3) :: ahf_centroid, rahf_centroid

      CCTK_REAL :: m1p, m2p, x1p, x2p,  y1p, y2p, z1p, z2p
      CCTK_REAL :: r12, r22, sigma
      CCTK_REAL, dimension(nx,ny,nz) :: fac
      CCTK_INT :: use_att, apower
      logical :: use_rot_att
      logical :: str_comp
      logical :: gaussf_exists
      external str_comp
!      CCTK_INT :: CCTK_IsFunctionAliased

      character(len=32) :: tmpstr

      character(len=200) :: gaussf

!     Declarations for macros.

c#include "CactusEinstein/ADMMacros/src/macro/ADM_Spacing_declare.h"
       CCTK_REAL :: dt
       CCTK_REAL :: idx, idy, idz
       CCTK_REAL :: i2dx, i2dy, i2dz
       CCTK_REAL :: i12dx, i12dy, i12dz
       CCTK_REAL :: idxx, idxy, idxz, idyy, idyz, idzz
       CCTK_REAL :: i12dxx, i12dyy, i12dzz
       CCTK_REAL :: i36dxy, i36dxz, i36dyz

#include "CactusEinstein/ADMMacros/src/macro/UPPERMET_declare.h"
#include "CactusEinstein/ADMMacros/src/macro/CHR2_declare.h"
#include "CactusEinstein/ADMMacros/src/macro/RICCI_declare.h"
#include "CactusEinstein/ADMMacros/src/macro/TRRICCI_declare.h"

!     Data.

      data firstgau    / .true. /
      data firstcal(1) / .true. /
      data firstcal(2) / .true. /
      data firstcal(3) / .true. /
      data firstcal(4) / .true. /

      save firstcal
      save firstgau

!     Description of variables:
!
!     i,j,k,
!     l,m,n,pp    Counters.
!
!     npoints    Number of points to interpolate.
!
!     error1     Different from zero if radius is negative.
!     error2     Different from zero if outside computational domain.
!
!     theta      Latitude.
!     phiangle        Longitude.
!
!     cost       cos(theta)
!     sint       sin(theta)
!     cosp       cos(phiangle)
!     sinp       sin(phiangle)
!
!     dtheta     Grid spacing in theta.
!     dphiangle       Grid spacing in phiangle.
!     dtp        dtheta*dphiangle.
!
!     phianglestart   Origin for phiangle (normally 0, but for some symmetries -pi/2).
!
!     idtheta    1/(2 dtheta)
!     idphiangle      1/(2 dphiangle)
!
!     xp         x coordinate from surface centre.
!     yp         y coordinate from surface centre.
!     zp         z coordinate from surface centre.
!
!     rp         Radius.
!     rr         Radius array.
!
!     xa         Array with x coordinates from grid centre.
!     ya         Array with y coordinates from grid centre.
!     za         Array with z coordinates from grid centre.
!
!     txx        Interpolated gxx metric component.
!     tyy        Interpolated gyy metric component.
!     tzz        Interpolated gzz metric component.
!     txy        Interpolated gxy metric component.
!     txz        Interpolated gxz metric component.
!     tyz        Interpolated gyz metric component.
!
!     trr        3-metric component {r,r}.
!     ttt        3-metric component {theta,theta}.
!     tpp        3-metric component {phiangle,phiangle}.
!     trt        3-metric component {r,theta}.
!     trp        3-metric component {r,phiangle}.
!     ttp        3-metric component {theta,phiangle}.
!
!     ft         dr/dtheta
!     fp         dr/dphiangle
!
!     g11        2-metric component {theta,theta}.
!     g22        2-metric component {phiangle,phiangle}.
!     g12        2-metric component {theta,phiangle}.
!
!     gaussian   Gaussian curvature.
!
!     ug         3x3 array with inverse 3-metric.
!
!     xi_ij      3x3 array for extrinsic curvature of level sets
!                of horizon function.
!     trxi       Trace of xi.
!     xi2        Square of xi.
!
!     ga_ij      2x2 array with angular metric,
!     ua_ij      2x2 array with inverse angular metric.
!
!     deta       Determinant of angular metric.
!     ideta      1/deta
!
!     d1a_ijk    d ga
!                 i  kl
!
!     d2a_ijkl   d d ga
!                 i j  kl
!
!     gammad     Christoffel symbol with 3 indices down.
!
!     gamma      Christoffel symbol with first index up.
!
!     gamma2     Product of Christoffel symbols.


!     **************************
!     ***   DEFINE NUMBERS   ***
!     **************************

      zero  = 0.0D0
      half  = 0.5D0
      one   = 1.0D0
      two   = 2.0D0
      three = 3.0D0
      four  = 4.0D0

      pi = acos(-one)

!     Number of points to interpolate.

      if (myproc.eq.0) then
         npoints = (ntheta+1)*(nphi+1)
      else
         npoints = 1
      end if

!     ******************************************************
!     ***   Get the reduction handle for sum operation   ***
!     ******************************************************

      call CCTK_ReductionArrayHandle(sum_handle,"sum")

      if (sum_handle.lt.0) then
        call CCTK_WARN(1,"Cannot get handle for sum reduction ! Forgot to activate an implementation providing reduction operators ??")
      end if

!     **************************************
!     ***   ALLOCATE MEMORY FOR ARRAYS   ***
!     **************************************

      if (myproc.eq.0) then

         allocate(rr(1:ntheta+1,1:nphi+1))

         allocate(xa(1:ntheta+1,1:nphi+1))
         allocate(ya(1:ntheta+1,1:nphi+1))
         allocate(za(1:ntheta+1,1:nphi+1))

         allocate(txx(1:ntheta+1,1:nphi+1))
         allocate(tyy(1:ntheta+1,1:nphi+1))
         allocate(tzz(1:ntheta+1,1:nphi+1))
         allocate(txy(1:ntheta+1,1:nphi+1))
         allocate(txz(1:ntheta+1,1:nphi+1))
         allocate(tyz(1:ntheta+1,1:nphi+1))

         allocate(g11(1:ntheta+1,1:nphi+1))
         allocate(g22(1:ntheta+1,1:nphi+1))
         allocate(g12(1:ntheta+1,1:nphi+1))

      else

         allocate(rr(1,1))

         allocate(xa(1,1))
         allocate(ya(1,1))
         allocate(za(1,1))

         allocate(txx(1,1))
         allocate(tyy(1,1))
         allocate(tzz(1,1))
         allocate(txy(1,1))
         allocate(txz(1,1))
         allocate(tyz(1,1))

         allocate(g11(1,1))
         allocate(g22(1,1))
         allocate(g12(1,1))

      end if

      rr = zero

      xa = zero
      ya = zero
      za = zero

      txx = zero
      tyy = zero
      tzz = zero
      txy = zero
      txz = zero
      tyz = zero

      g11 = zero
      g22 = zero
      g12 = zero

      gaussian = zero


!     *********************************
!     ***   INITIALIZE PARAMETERS   ***
!     *********************************

!     Find {phianglestart,dtheta,dphiangle,dtp}.

      phianglestart = zero

      dtheta = pi/dble(ntheta)

      if (cartoon) then
!         write(*,*) "ahfinder_gau cartoon ON: ",cartoon
         dphiangle    = zero
         dtp     = dtheta
         idphiangle   = one
         idtheta = half/dtheta

      else

         dphiangle = two*pi/dble(nphi)

         if (refz) dtheta = half*dtheta

         if (refx.and.refy) then
            dphiangle = 0.25D0*dphiangle
         else if (refx) then
            dphiangle = half*dphiangle
            phianglestart = - half*pi
         else if (refy) then
            dphiangle = half*dphiangle
         end if

         dtp = dtheta*dphiangle

         idtheta = half/dtheta
         idphiangle   = half/dphiangle

      end if

!     Initialize circumferences.

      circ_eq = zero
      meri_p1 = zero
      meri_p2 = zero

#include "CactusEinstein/ADMMacros/src/macro/ADM_Spacing.h"


!     *******************************************************
!     ***   FIND GAUSSIAN CURVATURE AS 3D GRID FUNCTION   ***
!     *******************************************************

!     Here I calculate the gaussian curvature of the level sets
!     of the horizon function as a 3D grid function.  Later I will
!     interpolate it onto the horizon.  This is cleaner than my old
!     method (see below).  It depends on less interpolations, and
!     the interpolations are done only after all derivatives have
!     been calculated.

      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1

!              Find covariant normal unit vector.

               aux = one/ahfgradn(i,j,k)

               nlx = ahfgradx(i,j,k)*aux
               nly = ahfgrady(i,j,k)*aux
               nlz = ahfgradz(i,j,k)*aux

!              Find upper spatial metric using standard macro.

#include "CactusEinstein/ADMMacros/src/macro/UPPERMET_guts.h"

               ug(1,1) = UPPERMET_UXX
               ug(2,2) = UPPERMET_UYY
               ug(3,3) = UPPERMET_UZZ

               ug(1,2) = UPPERMET_UXY
               ug(1,3) = UPPERMET_UXZ
               ug(2,3) = UPPERMET_UYZ

               ug(2,1) = ug(1,2)
               ug(3,1) = ug(1,3)
               ug(3,2) = ug(2,3)

!              Find contravariant normal unit vector.

               nux = ug(1,1)*nlx + ug(1,2)*nly + ug(1,3)*nlz
               nuy = ug(2,1)*nlx + ug(2,2)*nly + ug(2,3)*nlz
               nuz = ug(3,1)*nlx + ug(3,2)*nly + ug(3,3)*nlz

!              Find Christoffel symbols using standard macros.

#include "CactusEinstein/ADMMacros/src/macro/CHR2_guts.h"

!              Find extrinsic curvature of the 2-surfaces defined
!              as the level sets of the horizon function.  The
!              extrinsic curvature is defined as:
!                         __
!              xi   =  -  \/  n
!                ij         (i j)
!                    __
!              where \/ is the 3-dimensional covariant derivative and
!              n_i is the unit normal vector to the level sets of the
!              horizon function.

               xi(1,1) = - (ahfgradx(i+1,j,k)/ahfgradn(i+1,j,k)
     .                 - ahfgradx(i-1,j,k)/ahfgradn(i-1,j,k))/(two*dx)
     .                 + (CHR2_XXX*nlx + CHR2_YXX*nly + CHR2_ZXX*nlz)

               xi(2,2) = - (ahfgrady(i,j+1,k)/ahfgradn(i,j+1,k)
     .                 - ahfgrady(i,j-1,k)/ahfgradn(i,j-1,k))/(two*dy)
     .                 + (CHR2_XYY*nlx + CHR2_YYY*nly + CHR2_ZYY*nlz)

               xi(3,3) = - (ahfgradz(i,j,k+1)/ahfgradn(i,j,k+1)
     .                 - ahfgradz(i,j,k-1)/ahfgradn(i,j,k-1))/(two*dz)
     .                 + (CHR2_XZZ*nlx + CHR2_YZZ*nly + CHR2_ZZZ*nlz)

               xi(1,2) = - (ahfgradx(i,j+1,k)/ahfgradn(i,j+1,k)
     .                 - ahfgradx(i,j-1,k)/ahfgradn(i,j-1,k))/(four*dy)
     .                 -   (ahfgrady(i+1,j,k)/ahfgradn(i+1,j,k)
     .                 - ahfgrady(i-1,j,k)/ahfgradn(i-1,j,k))/(four*dx)
     .                 + (CHR2_XXY*nlx + CHR2_YXY*nly + CHR2_ZXY*nlz)

               xi(1,3) = - (ahfgradx(i,j,k+1)/ahfgradn(i,j,k+1)
     .                 - ahfgradx(i,j,k-1)/ahfgradn(i,j,k-1))/(four*dz)
     .                 -   (ahfgradz(i+1,j,k)/ahfgradn(i+1,j,k)
     .                 - ahfgradz(i-1,j,k)/ahfgradn(i-1,j,k))/(four*dx)
     .                 + (CHR2_XXZ*nlx + CHR2_YXZ*nly + CHR2_ZXZ*nlz)

               xi(2,3) = - (ahfgrady(i,j,k+1)/ahfgradn(i,j,k+1)
     .                 - ahfgrady(i,j,k-1)/ahfgradn(i,j,k-1))/(four*dz)
     .                 -   (ahfgradz(i,j+1,k)/ahfgradn(i,j+1,k)
     .                 - ahfgradz(i,j-1,k)/ahfgradn(i,j-1,k))/(four*dy)
     .                 + (CHR2_XYZ*nlx + CHR2_YYZ*nly + CHR2_ZYZ*nlz)

               xi(2,1) = xi(1,2)
               xi(3,1) = xi(1,3)
               xi(3,2) = xi(2,3)

!              Find trace of xi.

               aux = zero

               do l=1,3
                  do m=1,3
                     aux = aux + ug(l,m)*xi(l,m)
                  end do
               end do

               trxi = aux

!              Find square of xi.

               aux = zero

               do l=1,3
                  do m=1,3
                     do n=1,3
                        do pp=1,3
                           aux = aux + ug(l,m)*ug(n,pp)*xi(l,n)*xi(m,pp)
                        end do
                     end do
                  end do
               end do

               xi2 = aux

!              Find 3-Ricci and 3-Ricci scalar using standard macros.

#include "CactusEinstein/ADMMacros/src/macro/RICCI_guts.h"
#include "CactusEinstein/ADMMacros/src/macro/TRRICCI_guts.h"

!              Find 2-Ricci scalar using the contracted Gauss-Codazzi
!              relations:
!
!              (2)     (3)         a  b (3)              2           ab
!                 R  =    R  -  2 n  n     R   +  (tr xi)   -  xi  xi
!                                           ab                   ab

               aux = TRRICCI_TRRICCI - two*(nux**2*RICCI_RXX
     .             + nuy**2*RICCI_RYY + nuz**2*RICCI_RZZ
     .             + two*(nux*nuy*RICCI_RXY + nux*nuz*RICCI_RXZ
     .             + nuy*nuz*RICCI_RYZ))

               aux = aux + trxi**2 - xi2

!              Find gaussian curvature.  The gaussian curvature is
!              defined as R/2, with R the Ricci scalar of the surfaces.

               ahfgauss(i,j,k) = half*aux

!              Undefine macros!

#include "CactusEinstein/ADMMacros/src/macro/UPPERMET_undefine.h"
#include "CactusEinstein/ADMMacros/src/macro/CHR2_undefine.h"
#include "CactusEinstein/ADMMacros/src/macro/RICCI_undefine.h"
#include "CactusEinstein/ADMMacros/src/macro/TRRICCI_undefine.h"

            end do
         end do
      end do

!     Boundaries.

      ahfgauss(1,:,:) = ahfgauss(2,:,:)
      ahfgauss(:,1,:) = ahfgauss(:,2,:)
      ahfgauss(:,:,1) = ahfgauss(:,:,2)

      ahfgauss(nx,:,:) = ahfgauss(nx-1,:,:)
      ahfgauss(:,ny,:) = ahfgauss(:,ny-1,:)
      ahfgauss(:,:,nz) = ahfgauss(:,:,nz-1)

!     Synchronize.

      call CCTK_SyncGroup(ierror,cctkGH,"ahfinder::ahfinder_gauss")


!     *************************************
!     ***   FIND INTERPOLATING POINTS   ***
!     *************************************

!     Initialize {error1,error2}.

      error1 = 0
      error2 = 0

!     Notice that everything here is done on processor zero!

      if (myproc.eq.0) then

         avgx = 0.0D0
         avgy = 0.0D0
         avgz = 0.0D0

         do j=1,nphi+1
            do i=1,ntheta+1

!              Find {theta,phiangle}.

               theta = dtheta*dble(i-1)
               phiangle   = dphiangle*dble(j-1) + phianglestart

!              Find sines and cosines.

               cost = cos(theta)
               sint = sin(theta)

               cosp = cos(phiangle)
               sinp = sin(phiangle)

!              Find radius rp.

               rp = c0(0)

               do l=1+stepz,lmax,1+stepz
                  rp = rp + c0(l)*LEGEN(l,0,cost)
               end do

!              Notice how the sum over m is first.  This will allow
!              me to use the recursion relations to avoid having to
!              start from scratch every time.  Also, I sum over all
!              l s even if I do not want some terms.  This is
!              because in order to use the recursion relations I
!              need all polynomials.

               if (nonaxi) then
                  do m=1,lmax
                     aux = dble(m)*phiangle
                     sina = sin(aux)
                     cosa = cos(aux)
                     do l=m,lmax
                        aux = LEGEN(l,m,cost)
                        rp = rp + aux*cc(l,m)*cosa
                        if (.not.refy) then
                           rp = rp + aux*cs(l,m)*sina
                        end if
                     end do
                  end do
               end if

!              Check for negative radius.

               if (rp.le.zero) then
                  error1 = 1
               end if

!              Save radius array.

               rr(i,j) = rp

!              Find {xa,ya,za}

               xa(i,j) = rp*sint*cosp + xc
               ya(i,j) = rp*sint*sinp + yc
               za(i,j) = rp*cost + zc

!              Check if we are within bounds.

               xp = xa(i,j)
               yp = ya(i,j)
               zp = za(i,j)

               if ((xp.gt.xmx).or.(xp.lt.xmn).or.
     .             (yp.gt.ymx).or.(yp.lt.ymn).or.
     .             (zp.gt.zmx).or.(zp.lt.zmn)) then
                  error2 = 1
               end if

!              Add to sums for average position.

               avgx = avgx + xp
               avgy = avgy + yp
               avgz = avgz + zp

            end do
         end do

!        Find position of horizon centroid.

         avgx = avgx/dble(npoints)
         avgy = avgy/dble(npoints)
         avgz = avgz/dble(npoints)

         ahf_centroid(1) = avgx
         ahf_centroid(2) = avgy
         ahf_centroid(3) = avgz

!        Take symmetries into account.

         if (refx) ahf_centroid(1) = xc
         if (refy) ahf_centroid(2) = yc
         if (refz) ahf_centroid(3) = zc

!     Other processors.

      else

         xa = half*(xmx+xmn)
         ya = half*(ymx+ymn)
         za = half*(zmx+zmn)

         ahf_centroid(1) = zero
         ahf_centroid(2) = zero
         ahf_centroid(3) = zero
      end if


!     Reduce the errors across processors (all processors must
!     know about this since all will participate on the interpolation
!     below).

      call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,sum_handle,
     .   error1,rerror,CCTK_VARIABLE_INT)

      if (ierror.ne.0) then
         call CCTK_WARN(1,"Reduction failed!")
      end if
      error1 = rerror

      call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,sum_handle,
     .   error2,rerror,CCTK_VARIABLE_INT)

      if (ierror.ne.0) then
         call CCTK_WARN(1,"Reduction failed!")
      end if
      error2 = rerror

      call CCTK_ReduceLocArrayToArray1D(ierror,cctkGH,-1,sum_handle,
     .   ahf_centroid,rahf_centroid,3,CCTK_VARIABLE_REAL)

      if (ierror.ne.0) then
         call CCTK_WARN(1,"Reduction failed!")
      end if

      if ( ( horizon_to_announce_centroid .gt. 0 ) .and. 
     .     ( horizon_to_announce_centroid .eq. mfind ) ) then

         call CCTK_IsFunctionAliased(ierror, 
     .                                "SetDriftCorrectPosition")
         if ( ierror .eq. 1 ) then
            call CCTK_INFO("Announcing centroid to DriftCorrect")
            call SetDriftCorrectPosition ( cctkGH, ahf_centroid(1),
     .                                             ahf_centroid(2), 
     .                                             ahf_centroid(3) )
         end if
      end if

      if ( ( horizon_to_output_centroid .gt. 0 ) .and.
     .     ( horizon_to_output_centroid .eq. mfind ) ) then

         ahf_centroid_x = rahf_centroid(1)
         ahf_centroid_y = rahf_centroid(2)
         ahf_centroid_z = rahf_centroid(3)

         if (verbose) then
            write(*,*)
            write(6,*) 'AHFinder_gau: ahf_centroid_x = ',
     .                                             ahf_centroid_x
            write(6,*) 'AHFinder_gau: ahf_centroid_y = ',
     .                                             ahf_centroid_y
            write(6,*) 'AHFinder_gau: ahf_centroid_z = ',
     .                                             ahf_centroid_z
         end if
      end if


!     If there was an error then return from the subroutine
!     (but remember to deallocate the arrays first!).

      if ((error1.ne.0).or.(error2.ne.0)) then
         deallocate(rr)
         deallocate(xa,ya,za)
         deallocate(txx,tyy,tzz,txy,txz,tyz)
         deallocate(g11,g22,g12)
         return
      end if


!     ***********************
!     ***   INTERPOLATE   ***
!     ***********************

!     Here I should really only interpolate the gaussian curvature.
!     The interpolation of the metric is only needed for the old
!     gaussian curvature calculation (which I dont use any more),
!     and for the calculation of the circumferences (which really
!     should be done in the routine AHFinder_int.F).  I need to fix
!     this sometime soon! (October 2000).

!     parameter, local interpolator, and coordinate system handle.
      param_table_handle = -1
      interp_handle = -1
      coord_system_handle = -1

      options_string = "order = " // char(ichar('0') + interpolation_order)
      call Util_TableCreateFromString (param_table_handle, options_string)
      if (param_table_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot create parameter table for interpolator")
      endif

      call CCTK_FortranString (nchars, interpolation_operator, operator)
      call CCTK_InterpHandle (interp_handle, operator)
      if (interp_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
      endif

      call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
      if (coord_system_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot get handle for cart3d coordinate system ! Forgot to activate an implementation providing coordinates ??")
      endif


!     fill in the input/output arrays for the interpolator
      interp_coords(1) = CCTK_PointerTo(xa)
      interp_coords(2) = CCTK_PointerTo(ya)
      interp_coords(3) = CCTK_PointerTo(za)

      call CCTK_VarIndex (vindex, "bssn::gxx")
      in_array_indices(1) = vindex
      call CCTK_VarIndex (vindex, "bssn::gyy")
      in_array_indices(2) = vindex
      call CCTK_VarIndex (vindex, "bssn::gzz")
      in_array_indices(3) = vindex
      call CCTK_VarIndex (vindex, "bssn::gxy")
      in_array_indices(4) = vindex
      call CCTK_VarIndex (vindex, "bssn::gxz")
      in_array_indices(5) = vindex
      call CCTK_VarIndex (vindex, "bssn::gyz")
      in_array_indices(6) = vindex
      call CCTK_VarIndex (vindex, "ahfinder::ahfgauss")
      in_array_indices(7) = vindex

      out_arrays(1) = CCTK_PointerTo(txx)
      out_arrays(2) = CCTK_PointerTo(tyy)
      out_arrays(3) = CCTK_PointerTo(tzz)
      out_arrays(4) = CCTK_PointerTo(txy)
      out_arrays(5) = CCTK_PointerTo(txz)
      out_arrays(6) = CCTK_PointerTo(tyz)
      out_arrays(7) = CCTK_PointerTo(gaussian)

      out_array_type_codes = CCTK_VARIABLE_REAL


!     Interpolation.
      call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     .                            param_table_handle, coord_system_handle,
     .                            npoints, CCTK_VARIABLE_REAL, interp_coords,
     .                            7, in_array_indices,
     .                            7, out_array_type_codes, out_arrays)
      if (ierror < 0) then
        call CCTK_WARN (1, "AHFinder_gau.F: interpolator call returned an error code");

        stop
!        write(*,*) "out_arraystxx = ",txx
!        write(*,*) "out_arraystyy = ",tyy
!        write(*,*) "out_arraystzz = ",tzz
!        write(*,*) "out_arraystxy = ",txy
!        write(*,*) "out_arraystxz = ",txz
!        write(*,*) "out_arraystyz = ",tyz
!        write(*,*) "out_arraysgaussian = ",gaussian
      endif

!     release parameter table
      call Util_TableDestroy (ierror, param_table_handle)


!     *********************************
!     ***   FIND SPHERICAL METRIC   ***
!     *********************************

      if (myproc.eq.0) then

         do j=1,nphi+1
            do i=1,ntheta+1

!              Find {theta,phiangle}.

               theta = dtheta*dble(i-1)
               phiangle   = dphiangle*dble(j-1) + phianglestart

!              Find {rp}.

               rp = rr(i,j)

!              Find sines and cosines.

               cost = cos(theta)
               sint = sin(theta)

               cosp = cos(phiangle)
               sinp = sin(phiangle)

!              Find metric in spherical coordinates.

               trr = sint**2*(txx(i,j)*cosp**2
     .             + tyy(i,j)*sinp**2) + tzz(i,j)*cost**2
     .             + two*sint*(txy(i,j)*sint*cosp*sinp
     .             + cost*(txz(i,j)*cosp + tyz(i,j)*sinp))

               ttt = rp**2*(cost**2*(txx(i,j)*cosp**2
     .             + tyy(i,j)*sinp**2) + tzz(i,j)*sint**2
     .             + two*cost*(txy(i,j)*cost*cosp*sinp
     .             - sint*(txz(i,j)*cosp + tyz(i,j)*sinp)))

               tpp = (rp*sint)**2*(txx(i,j)*sinp**2 + tyy(i,j)*cosp**2
     .             - two*txy(i,j)*cosp*sinp)

               trt = rp*(sint*cost*(txx(i,j)*cosp**2 + tyy(i,j)*sinp**2
     .             - tzz(i,j) + two*txy(i,j)*sinp*cosp)
     .             + (cost**2-sint**2)*(txz(i,j)*cosp + tyz(i,j)*sinp))

               trp = rp*sint*(sint*(cosp*sinp*(tyy(i,j) - txx(i,j))
     .             + txy(i,j)*(cosp**2 - sinp**2))
     .             + cost*(cosp*tyz(i,j) - sinp*txz(i,j)))

               ttp = rp**2*sint*(cost*(sinp*cosp*(tyy(i,j) - txx(i,j))
     .             + (cosp**2 - sinp**2)*txy(i,j))
     .             + sint*(sinp*txz(i,j) - cosp*tyz(i,j)))

!              Find derivatives  {ft,fp}.  For interior points I
!              use centered differences, and for the boundary points
!              I use second order one-sided differences.

               if ((i.ne.1).and.(i.ne.ntheta+1)) then
                  ft = idtheta*(rr(i+1,j) - rr(i-1,j))
               else if (i.eq.1) then
                  ft = - idtheta*(three*rr(1,j)
     .               - four*rr(2,j) + rr(3,j))
               else
                  ft = + idtheta*(three*rr(ntheta,j)
     .               - four*rr(ntheta,j) + rr(ntheta-1,j))
               end if

               if (nonaxi) then
                  if ((j.ne.1).and.(j.ne.nphi+1)) then
                     fp = idphiangle*(rr(i,j+1) - rr(i,j-1))
                  else if (j.eq.1) then
                     fp = - idphiangle*(three*rr(i,1)
     .                  - four*rr(i,2) + rr(i,3))
                  else
                     fp = + idphiangle*(three*rr(i,nphi)
     .                  - four*rr(i,nphi) + rr(i,nphi-1))
                  end if
               else
                  fp = zero
               end if

!              Find induced metric on the surface.

               g11(i,j) = ttt + trr*ft**2 + two*trt*ft
               g22(i,j) = tpp + trr*fp**2 + two*trp*fp
               g12(i,j) = ttp + trr*ft*fp + trp*ft + trt*fp

            end do
         end do

      end if


!     ************************************
!     ***   CALCULATE CIRCUMFERENCES   ***
!     ************************************

      if (myproc.eq.0) then

!        Equatorial circumference.

         if (refz) then

            do j=1,nphi
               circ_eq = circ_eq + dphiangle
     .            *sqrt(abs(half*(g22(ntheta+1,j) + g22(ntheta+1,j+1))))
            end do

         else

            aux = half*pi/dtheta
            i = int(aux)
            aux = aux - int(aux)

            if (cartoon) then
               circ_eq = 6.283185307D0*sqrt(g22(i,0))
            else
               do j=1,nphi
                  circ_eq = circ_eq + dphiangle
     .               *((one-aux)*sqrt(abs(half*(g22(i,j) + g22(i,j+1))))
     .               + aux*sqrt(abs(half*(g22(i+1,j) + g22(i+1,j+1)))))
               end do
            end if
         end if

!        Meridians.

         do i=1,ntheta
            meri_p1 = meri_p1 + dtheta
     .         *sqrt(abs(half*(g11(i,1) + g11(i+1,1))))
         end do

         if (cartoon) then
            meri_p2 = meri_p1
         else if (refx.and.refy) then
            do i=1,ntheta
               meri_p2 = meri_p2 + dtheta
     .            *sqrt(abs(half*(g11(i,nphi+1) + g11(i+1,nphi+1))))
            end do
         else
            aux = half*pi/dphiangle
            j = int(aux)
            aux = aux - int(aux)
            do i=1,ntheta
               meri_p2 = meri_p2 + dtheta
     .            *((one-aux)*sqrt(abs(half*(g11(i,j) + g11(i,j))))
     .            + aux*sqrt(abs(half*(g11(i,j+1) + g11(i,j+1)))))
            end do
         end if

!        Rescale the integrals according to symmetries.

         if (.not.cartoon) then
            if (refx) circ_eq = two*circ_eq
            if (refy) circ_eq = two*circ_eq
         end if

         if (refz) then
            meri_p1 = two*meri_p1
            meri_p2 = two*meri_p2
         end if

      end if


!     *************************************************
!     ***   OLD CALCULATION OF GAUSSIAN CURVATURE   ***
!     *************************************************

!     This is the old calculation of the gaussian curvature.
!     It is correct as far as I know, but it depends on taking
!     second derivatives of the interpolated metric, and this
!     seems not to behave well numerically (we would need higher
!     order interpolation).  I did not want to delete this code
!     since it might be useful in the future, so here I just jump
!     over it.

      if (.false.) then

!     The Gaussian curvature of the surface is defined as R/2
!     where R is the Ricci scalar of the two-geometry.  I can
!     now calculate the Ricci scalar using the spherical metric
!     components {ttt,tpp,ttp}.

!     Initialize arrays.

      gaussian = zero

      ga(1,1) = one
      ga(2,2) = one
      ga(1,2) = zero
      ga(2,1) = zero

      ua(1,1) = one
      ua(2,2) = one
      ua(1,2) = zero
      ua(2,1) = zero

      d1a = zero
      d2a = zero

      gamma  = zero
      gamma2 = zero
      gammad = zero

      if (myproc.eq.0) then

!        Interior points.

         do j=2,nphi
            do i=2,ntheta

!              Find angular metric.

               ga(1,1) = g11(i,j)
               ga(2,2) = g22(i,j)
               ga(1,2) = g12(i,j)
               ga(2,1) = g12(i,j)

!              Find determinant of angular metric.

               deta  = ga(1,1)*ga(2,2) - ga(1,2)**2
               ideta = one/deta

!              Find inverse angular metric.

               ua(1,1) = + ga(2,2)*ideta
               ua(2,2) = + ga(1,1)*ideta
               ua(1,2) = - ga(1,2)*ideta
               ua(2,1) = - ga(2,1)*ideta

!              First derivatives of angular metric.

               d1a(1,1,1) = idtheta*(g11(i+1,j) - g11(i-1,j))
               d1a(1,2,2) = idtheta*(g22(i+1,j) - g22(i-1,j))
               d1a(1,1,2) = idtheta*(g12(i+1,j) - g12(i-1,j))
               d1a(1,2,1) = d1a(1,1,2)

               if (nonaxi) then
                  d1a(2,1,1) = idphiangle*(g11(i,j+1) - g11(i,j-1))
                  d1a(2,2,2) = idphiangle*(g22(i,j+1) - g22(i,j-1))
                  d1a(2,1,2) = idphiangle*(g12(i,j+1) - g12(i,j-1))
                  d1a(2,2,1) = d1a(2,1,2)
               else
                  d1a(2,1,1) = zero
                  d1a(2,2,2) = zero
                  d1a(2,1,2) = zero
                  d1a(2,2,1) = zero
               end if

!              Second derivatives of angular metric.

               d2a(1,1,1,1) = four*idtheta**2*(g11(i+1,j)
     .                      - two*g11(i,j) + g11(i-1,j))
               d2a(1,1,2,2) = four*idtheta**2*(g22(i+1,j)
     .                      - two*g22(i,j) + g22(i-1,j))
               d2a(1,1,1,2) = four*idtheta**2*(g12(i+1,j)
     .                      - two*g12(i,j) + g12(i-1,j))

               if (nonaxi) then
                  d2a(2,2,1,1) = four*idphiangle**2*(g11(i,j+1)
     .                         - two*g11(i,j) + g11(i,j-1))
                  d2a(2,2,2,2) = four*idphiangle**2*(g22(i,j+1)
     .                         - two*g22(i,j) + g22(i,j-1))
                  d2a(2,2,1,2) = four*idphiangle**2*(g12(i,j+1)
     .                         - two*g12(i,j) + g12(i,j-1))

                  d2a(1,2,1,1) = (g11(i+1,j+1) - g11(i+1,j-1)
     .                         - g11(i-1,j+1) + g11(i-1,j-1))
     .                         *idtheta*idphiangle
                  d2a(1,2,2,2) = (g22(i+1,j+1) - g22(i+1,j-1)
     .                         - g22(i-1,j+1) + g22(i-1,j-1))
     .                         *idtheta*idphiangle
                  d2a(1,2,1,2) = (g12(i+1,j+1) - g12(i+1,j-1)
     .                         - g12(i-1,j+1) + g12(i-1,j-1))
     .                         *idtheta*idphiangle
               else
                  d2a(2,2,1,1) = zero
                  d2a(2,2,2,2) = zero
                  d2a(2,2,1,2) = zero
                  d2a(1,2,1,1) = zero
                  d2a(1,2,2,2) = zero
                  d2a(1,2,1,2) = zero
               end if

               d2a(1,1,2,1) = d2a(1,1,1,2)
               d2a(2,2,2,1) = d2a(2,2,1,2)
               d2a(1,2,2,1) = d2a(1,2,1,2)

               d2a(2,1,1,1) = d2a(1,2,1,1)
               d2a(2,1,2,2) = d2a(1,2,2,2)
               d2a(2,1,1,2) = d2a(1,2,1,2)
               d2a(2,1,2,1) = d2a(1,2,2,1)

!              Add correction terms to derivatives of g_{phiangle,phiangle}.
!              Near the axis, it turns out that the lower order
!              contributions to the derivatives of g_{phiangle,phiangle}
!              cancel out, and the gaussian curvature is made only
!              of higher order terms that the centered finite difference
!              approximations get wrong!  So here I add higher order
!              corrections to make sure that I get correctly the
!              first two terms of the Taylor series (assuming that
!              for small theta  g_{phiangle,phiangle} = f(r,phiangle) sin(theta)^2).

               theta = dtheta*dble(i-1)

               d1a(1,2,2) = d1a(1,2,2) + (two/three)*sin(theta)
     .                    *(g22(i+1,j) - two*g22(i,j) + g22(i-1,j))

               d2a(1,1,2,2) = d2a(1,1,2,2) + one/three
     .                    *(g22(i+1,j) - two*g22(i,j) + g22(i-1,j))

               if (nonaxi) then
                  d2a(1,2,2,2) = d2a(1,2,2,2) + (two/three)*sin(theta)
     .               *((g22(i+1,j+1) - two*g22(i,j+1) + g22(i-1,j+1))
     .               - (g22(i+1,j-1) - two*g22(i,j-1) + g22(i-1,j-1)))
     .               *idphiangle
                  d2a(2,1,2,2) = d2a(1,2,2,2)
               end if

!              Find Christoffel symbols.

               do l=1,2
                  do m=1,2
                     do n=1,2
                        gammad(l,m,n) = half*(d1a(m,n,l) + d1a(n,m,l)
     .                                - d1a(l,m,n))
                     end do
                  end do
               end do

               do l=1,2
                  do m=1,2
                     do n=1,2
                        aux = zero
                        do k=1,2
                           aux = aux + ua(l,k)*gammad(k,m,n)
                        end do
                        gamma(l,m,n) = aux
                     end do
                  end do
               end do

               do k=1,2
                  do l=1,2
                     do m=1,2
                        do n=1,2
                           aux = zero
                           do pp=1,2
                              aux = aux + gamma(pp,k,l)*gammad(pp,m,n)
                           end do
                           gamma2(k,l,m,n) = aux
                        end do
                     end do
                  end do
               end do

!              Find gaussian curvature.

               aux = zero

               do k=1,2
                  do l=1,2
                     do m=1,2
                        do n=1,2
                           aux = aux + ua(k,l)*ua(m,n)
     .                         *(d2a(k,n,l,m) - d2a(k,l,m,n)
     .                         + gamma2(k,m,l,n) - gamma2(k,l,m,n))
                        end do
                     end do
                  end do
               end do

               gaussian(i,j) = half*aux

            end do
         end do

!        Boundaries.

         gaussian(1,:) = two*gaussian(2,:) - gaussian(3,:)
         gaussian(ntheta+1,:) = two*gaussian(ntheta,:)
     .                      - gaussian(ntheta-1,:)

         gaussian(:,1) = two*gaussian(:,2) - gaussian(:,3)
         gaussian(:,nphi+1) = two*gaussian(:,nphi)
     .                    - gaussian(:,nphi-1)

      end if

      end if


!     ************************************
!     ***   INITIALIZE FIRSTGAU FLAG   ***
!     ************************************

      if (find3) then
         firstgau = firstcal(mfind)
         firstcal(mfind) = .false.
      else
         firstgau = firstcal(4)
         firstcal(4) = .false.
      end if


!     ***********************************
!     ***   SAVE GAUSSIAN CURVATURE   ***
!     ***********************************

      if (myproc.eq.0) then

         if (find3) then
            if (mfind.eq.1) then
               gaussf = filestr(1:nfile)//"/ahf_0.gauss"
            else if (mfind.eq.2) then
               gaussf = filestr(1:nfile)//"/ahf_1.gauss"
            else
               gaussf = filestr(1:nfile)//"/ahf_2.gauss"
            end if
         else
            gaussf = filestr(1:nfile)//"/ahf.gauss"
         end if

         inquire(file=gaussf, exist=gaussf_exists)

         if (firstgau) then
            if (find3) then
               if (mfind.eq.3) then
                  firstgau = .false.
               end if
            else
               firstgau = .false.
            end if
         end if

!        On first call replace file.

         if (.not.gaussf_exists) then

           call AHFinder_WriteGaussHeader(CCTK_PASS_FTOF, gaussf)

!        On subsequent calls, append to file.

         else

            open(11,file=gaussf,form='formatted', status='old',
     .          position='append')

            write(11,*)
            write(11,"(A11,I4)")    '# Time step',cctk_iteration
            write(11,"(A6,ES11.3)") '# Time',cctk_time
            write(11,"(A6,I4)")     '# Call',ahf_ncall

         end if

!        Save data points.

         do i=1,ntheta+1
            do j=1,nphi+1
               write(11,"(ES11.3)") gaussian(i,j)
            end do
         end do

!        Close file.

         close(11)

      end if


!     *****************************
!     ***   DEALLOCATE MEMORY   ***
!     *****************************

      deallocate(rr)
      deallocate(xa,ya,za)
      deallocate(txx,tyy,tzz,txy,txz,tyz)
      deallocate(g11,g22,g12)


!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder_gau


!     ******************************************
!     *** String comparison utility function ***
!     ******************************************

      function str_comp(a1,a2)

      logical :: str_comp
      character(len=*),intent(in) :: a1, a2

      if ( ( verify(trim(a1),trim(a2)).eq.0 )  .and.
     .                     ( len_trim(a1) .eq. len_trim(a2) ) ) then
        str_comp = .true.
      else
        str_comp = .false.
      end if

      end function str_comp
