/*@@
  @file      AHFinder_exp.F
  @date      April 1998
  @author    Miguel Alcubierre
  @desc 
            Here I calculate the (normalized) expansion of the
            horizon function. The expansion is defined by the
            following expression:

                    __2                     2 /  ab  __a__b      \
            exp  =  \/ f / u  +  d f d f / u  | K  - \/ \/ f / u |  -  trK
                                  a   b       \                  /

            where:

                  /  mn         \ 1/2
            u  =  | g   d f d f |
                  \      m   n  /

  @enddesc
  @version   $Header: /cactus/CactusEinstein/AHFinder/src/AHFinder_exp.F,v 1.18 2002/10/08 20:43:33 miguel Exp $
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      subroutine AHFinder_exp(CCTK_ARGUMENTS)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j,k
      integer CCTK_Equals,index

      CCTK_REAL det,idet
      CCTK_REAL gdxx,gdyy,gdzz,gdxy,gdxz,gdyz
      CCTK_REAL guxx,guyy,guzz,guxy,guxz,guyz
      CCTK_REAL kdxx,kdyy,kdzz,kdxy,kdxz,kdyz
      CCTK_REAL ddxxx,ddxyy,ddxzz,ddxxy,ddxxz,ddxyz
      CCTK_REAL ddyxx,ddyyy,ddyzz,ddyxy,ddyxz,ddyyz
      CCTK_REAL ddzxx,ddzyy,ddzzz,ddzxy,ddzxz,ddzyz
      CCTK_REAL d2fxx,d2fyy,d2fzz,d2fxy,d2fxz,d2fyz
      CCTK_REAL c2fxx,c2fyy,c2fzz,c2fxy,c2fxz,c2fyz
      CCTK_REAL dfx,dfy,dfz,dfux,dfuy,dfuz
      CCTK_REAL aux,T0,T1,T2,T3,T4
      CCTK_REAL idx,idy,idz
      CCTK_REAL idx2,idy2,idz2,idxy,idxz,idyz
      CCTK_REAL zero,half,one,two

!     Description of variables:
!
!     i,j,k,          counters
!
!     gdij            g
!                      ij
!                      ij
!     guij            g
!
!
!     det             det(g)
!
!     idet            1/det
!
!
!     kdij            K
!                      ij
!
!     ddijk           d g
!                      i jk
!
!     dfi             d f
!                      i
!                      ij
!     dfui            g  d f
!                         j
!
!     d2fij           d d f
!                      i j
!
!     c2fij           f
!                      ;ij
!
!     aux,T*          auxiliary variables


!     *******************
!     ***   NUMBERS   ***
!     *******************

      zero = 0.0d0
      half = 0.5d0
      one  = 1.0d0
      two  = 2.0d0


!     **************************
!     ***   FIND EXPANSION   ***
!     **************************

      idx = half/dx
      idy = half/dy
      idz = half/dz

      idx2 = one/dx**2
      idy2 = one/dy**2
      idz2 = one/dz**2

      idxy = idx*idy
      idxz = idx*idz
      idyz = idy*idz

!      stop

      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1

!              Find spatial metric.

               gdxx = gxx(i,j,k)
               gdyy = gyy(i,j,k)
               gdzz = gzz(i,j,k)
               gdxy = gxy(i,j,k)
               gdxz = gxz(i,j,k)
               gdyz = gyz(i,j,k)

!              Find extrinsic curvature.

               kdxx = kxx(i,j,k)
               kdyy = kyy(i,j,k)
               kdzz = kzz(i,j,k)
               kdxy = kxy(i,j,k)
               kdxz = kxz(i,j,k)
               kdyz = kyz(i,j,k)

!              Find derivatives of metric using finite differences.

               ddxxx = idx*(gxx(i+1,j,k) - gxx(i-1,j,k))
               ddxyy = idx*(gyy(i+1,j,k) - gyy(i-1,j,k))
               ddxzz = idx*(gzz(i+1,j,k) - gzz(i-1,j,k))
               ddxxy = idx*(gxy(i+1,j,k) - gxy(i-1,j,k))
               ddxxz = idx*(gxz(i+1,j,k) - gxz(i-1,j,k))
               ddxyz = idx*(gyz(i+1,j,k) - gyz(i-1,j,k))

               ddyxx = idy*(gxx(i,j+1,k) - gxx(i,j-1,k))
               ddyyy = idy*(gyy(i,j+1,k) - gyy(i,j-1,k))
               ddyzz = idy*(gzz(i,j+1,k) - gzz(i,j-1,k))
               ddyxy = idy*(gxy(i,j+1,k) - gxy(i,j-1,k))
               ddyxz = idy*(gxz(i,j+1,k) - gxz(i,j-1,k))
               ddyyz = idy*(gyz(i,j+1,k) - gyz(i,j-1,k))

               ddzxx = idz*(gxx(i,j,k+1) - gxx(i,j,k-1))
               ddzyy = idz*(gyy(i,j,k+1) - gyy(i,j,k-1))
               ddzzz = idz*(gzz(i,j,k+1) - gzz(i,j,k-1))
               ddzxy = idz*(gxy(i,j,k+1) - gxy(i,j,k-1))
               ddzxz = idz*(gxz(i,j,k+1) - gxz(i,j,k-1))
               ddzyz = idz*(gyz(i,j,k+1) - gyz(i,j,k-1))

!              Find determinant of spatial metric.

               det = gdxx*gdyy*gdzz + two*gdxy*gdxz*gdyz &
                    - gdxx*gdyz**2 - gdyy*gdxz**2 - gdzz*gdxy**2

!               write(*,*) det,gdxx,gdyy,gdzz,gdxy,gdxz,gdyz, &
!                  gdxx,gdyz,gdyy,gdxz,gdzz,gdxy

!              If determinant is not zero proceed.

               if (det.gt.zero) then

                  idet = one/det

!                 Find inverse spatial metric.

                  guxx = idet*(gdyy*gdzz - gdyz**2)
                  guyy = idet*(gdxx*gdzz - gdxz**2)
                  guzz = idet*(gdxx*gdyy - gdxy**2)

                  guxy = idet*(gdxz*gdyz - gdzz*gdxy)
                  guxz = idet*(gdxy*gdyz - gdyy*gdxz)
                  guyz = idet*(gdxy*gdxz - gdxx*gdyz)

!                 Find spatial derivatives of f.
                  
                  T0 = two*ahfgrid(i,j,k)
              
                  T1 = ahfgrid(i+1,j,k)
                  T2 = ahfgrid(i-1,j,k)
           
                  dfx   = (T1 - T2)*idx
                  d2fxx = (T1 - T0 + T2)*idx2

                  T1 = ahfgrid(i,j+1,k)
                  T2 = ahfgrid(i,j-1,k)

                  dfy   = (T1 - T2)*idy
                  d2fyy = (T1 - T0 + T2)*idy2

                  T1 = ahfgrid(i,j,k+1)
                  T2 = ahfgrid(i,j,k-1)

                  dfz   = (T1 - T2)*idz
                  d2fzz = (T1 - T0 + T2)*idz2

!                 Save gradient of horizon function and its norm
!                 (they will be needed later in the flow algorithm
!                 and in the gaussian curvature).  Here I also try
!                 to avoid possible division by zero later by not
!                 allowing   ahfgradn(i,j,k) = 0.  This should only
!                 ever happen far from the horizon, so resetting this
!                 to 1 should have no important effects.

                  ahfgradx(i,j,k) = dfx
                  ahfgrady(i,j,k) = dfy
                  ahfgradz(i,j,k) = dfz

                  aux = guxx*dfx**2 + guyy*dfy**2 + guzz*dfz**2 &
                       + two*(guxy*dfx*dfy + guxz*dfx*dfz + guyz*dfy*dfz)

                  ahfgradn(i,j,k) = sqrt(aux)
                  if (ahfgradn(i,j,k).eq.zero) ahfgradn(i,j,k) = one

!                 Find crossed derivatives.

                  d2fxy = (ahfgrid(i+1,j+1,k) + ahfgrid(i-1,j-1,k) &
                       -  ahfgrid(i+1,j-1,k) - ahfgrid(i-1,j+1,k))*idxy 
                  d2fxz = (ahfgrid(i+1,j,k+1) + ahfgrid(i-1,j,k-1) &
                       -  ahfgrid(i+1,j,k-1) - ahfgrid(i-1,j,k+1))*idxz
                  d2fyz = (ahfgrid(i,j+1,k+1) + ahfgrid(i,j-1,k-1) &
                       -  ahfgrid(i,j+1,k-1) - ahfgrid(i,j-1,k+1))*idyz

!                 Raise indices in first derivatives.

                  dfux = guxx*dfx + guxy*dfy + guxz*dfz
                  dfuy = guxy*dfx + guyy*dfy + guyz*dfz
                  dfuz = guxz*dfx + guyz*dfy + guzz*dfz
                  
!                 Find second covariant derivatives of f.

                  c2fxx = d2fxx - half*(dfux*ddxxx &
                       + dfuy*(two*ddxxy - ddyxx) &
                       + dfuz*(two*ddxxz - ddzxx))
                  
                  c2fyy = d2fyy - half*(dfuy*ddyyy &
                       + dfux*(two*ddyxy - ddxyy) &
                       + dfuz*(two*ddyyz - ddzyy))

                  c2fzz = d2fzz - half*(dfuz*ddzzz &
                       + dfux*(two*ddzxz - ddxzz) &
                       + dfuy*(two*ddzyz - ddyzz))

                  c2fxy = d2fxy - half*(dfux*ddyxx + dfuy*ddxyy &
                       + dfuz*(ddxyz + ddyxz - ddzxy))

                  c2fxz = d2fxz - half*(dfux*ddzxx + dfuz*ddxzz &
                       + dfuy*(ddxyz + ddzxy - ddyxz))

                  c2fyz = d2fyz - half*(dfuy*ddzyy + dfuz*ddyzz &
                       + dfux*(ddyxz + ddzxy - ddxyz))

!                              /  m      \ 1/2
!                 Find:  u  =  | d f d f |
!                              \      m  /

                  T0 = dfx*dfux + dfy*dfuy + dfz*dfuz

!                 if T0 is positive proceed.

                  if (T0.gt.zero) then

                     T0 = one/sqrt(T0)

!                           __2
!                    Find:  \/ f / u

                     T1 = guxx*c2fxx + guyy*c2fyy + guzz*c2fzz &
                          + two*(guxy*c2fxy + guxz*c2fxz + guyz*c2fyz)

                     T1 = T1*T0

!                                a   b       2
!                    Find:  K   d f d f  /  u
!                            ab
                     
                     T2 = kdxx*dfux**2 + kdyy*dfuy**2 + kdzz*dfuz**2 &
                          + two*(dfux*(kdxy*dfuy + kdxz*dfuz) + kdyz*dfuy*dfuz)

                     T2 = T2*T0**2

!                            __a__b                3
!                    Find:   \/ \/ f  d f d f  /  u
!                                      a   b

                     T3 = c2fxx*dfux**2 + c2fyy*dfuy**2 + c2fzz*dfuz**2 &
                          + two*(dfux*(c2fxy*dfuy + c2fxz*dfuz) &
                          + c2fyz*dfuy*dfuz)

                     T3 = T3*T0**3

!                    Find:  trK

                     T4 = guxx*kdxx + guyy*kdyy + guzz*kdzz &
                          + two*(guxy*kdxy + guxz*kdxz + guyz*kdyz)

!                    Find the expansion.

                     ahf_exp(i,j,k) = T1 + T2 - T3 - T4

!                 T0 was not positive.

                  else
                     ahf_exp(i,j,k) = zero
                  end if

!              Determinant was zero.

               else
                  ahf_exp(i,j,k) = zero
               end if

            end do
         end do
      end do

!     Zero out the edges before doing the boundaries.

      ahf_exp( 1, 1, :) = 0.0d0
      ahf_exp( 1,ny, :) = 0.0d0
      ahf_exp( 1, :, 1) = 0.0d0
      ahf_exp( 1, :,nz) = 0.0d0
      ahf_exp(nx, 1, :) = 0.0d0
      ahf_exp(nx,ny, :) = 0.0d0
      ahf_exp(nx, :, 1) = 0.0d0
      ahf_exp(nx, :,nz) = 0.0d0
      ahf_exp( :, 1, 1) = 0.0d0
      ahf_exp( :,ny, 1) = 0.0d0
      ahf_exp( :, 1,nz) = 0.0d0
      ahf_exp( :,ny,nz) = 0.0d0

      ahfgradx( 1, 1, :) = 0.0d0
      ahfgradx( 1,ny, :) = 0.0d0
      ahfgradx( 1, :, 1) = 0.0d0
      ahfgradx( 1, :,nz) = 0.0d0
      ahfgradx(nx, 1, :) = 0.0d0
      ahfgradx(nx,ny, :) = 0.0d0
      ahfgradx(nx, :, 1) = 0.0d0
      ahfgradx(nx, :,nz) = 0.0d0
      ahfgradx( :, 1, 1) = 0.0d0
      ahfgradx( :,ny, 1) = 0.0d0
      ahfgradx( :, 1,nz) = 0.0d0
      ahfgradx( :,ny,nz) = 0.0d0

      ahfgrady( 1, 1, :) = 0.0d0
      ahfgrady( 1,ny, :) = 0.0d0
      ahfgrady( 1, :, 1) = 0.0d0
      ahfgrady( 1, :,nz) = 0.0d0
      ahfgrady(nx, 1, :) = 0.0d0
      ahfgrady(nx,ny, :) = 0.0d0
      ahfgrady(nx, :, 1) = 0.0d0
      ahfgrady(nx, :,nz) = 0.0d0
      ahfgrady( :, 1, 1) = 0.0d0
      ahfgrady( :,ny, 1) = 0.0d0
      ahfgrady( :, 1,nz) = 0.0d0
      ahfgrady( :,ny,nz) = 0.0d0

      ahfgradz( 1, 1, :) = 0.0d0
      ahfgradz( 1,ny, :) = 0.0d0
      ahfgradz( 1, :, 1) = 0.0d0
      ahfgradz( 1, :,nz) = 0.0d0
      ahfgradz(nx, 1, :) = 0.0d0
      ahfgradz(nx,ny, :) = 0.0d0
      ahfgradz(nx, :, 1) = 0.0d0
      ahfgradz(nx, :,nz) = 0.0d0
      ahfgradz( :, 1, 1) = 0.0d0
      ahfgradz( :,ny, 1) = 0.0d0
      ahfgradz( :, 1,nz) = 0.0d0
      ahfgradz( :,ny,nz) = 0.0d0

      ahfgradn( 1, 1, :) = 0.0d0
      ahfgradn( 1,ny, :) = 0.0d0
      ahfgradn( 1, :, 1) = 0.0d0
      ahfgradn( 1, :,nz) = 0.0d0
      ahfgradn(nx, 1, :) = 0.0d0
      ahfgradn(nx,ny, :) = 0.0d0
      ahfgradn(nx, :, 1) = 0.0d0
      ahfgradn(nx, :,nz) = 0.0d0
      ahfgradn( :, 1, 1) = 0.0d0
      ahfgradn( :,ny, 1) = 0.0d0
      ahfgradn( :, 1,nz) = 0.0d0
      ahfgradn( :,ny,nz) = 0.0d0

!     Boundaries on x direction.

      ahf_exp(1,2:ny-1,2:nz-1)  = 2.0D0*ahf_exp(2,2:ny-1,2:nz-1) &
           - ahf_exp(3,2:ny-1,2:nz-1)
      ahfgradx(1,2:ny-1,2:nz-1) = 2.0D0*ahfgradx(2,2:ny-1,2:nz-1) &
           - ahfgradx(3,2:ny-1,2:nz-1)
      ahfgrady(1,2:ny-1,2:nz-1) = 2.0D0*ahfgrady(2,2:ny-1,2:nz-1) &
           - ahfgrady(3,2:ny-1,2:nz-1)
      ahfgradz(1,2:ny-1,2:nz-1) = 2.0D0*ahfgradz(2,2:ny-1,2:nz-1) &
           - ahfgradz(3,2:ny-1,2:nz-1)
      ahfgradn(1,2:ny-1,2:nz-1) = 2.0D0*ahfgradn(2,2:ny-1,2:nz-1) &
           - ahfgradn(3,2:ny-1,2:nz-1)

      ahf_exp(nx,2:ny-1,2:nz-1)  = 2.0D0*ahf_exp(nx-1,2:ny-1,2:nz-1) &
           - ahf_exp(nx-2,2:ny-1,2:nz-1)
      ahfgradx(nx,2:ny-1,2:nz-1) = 2.0D0*ahfgradx(nx-1,2:ny-1,2:nz-1) &
           - ahfgradx(nx-2,2:ny-1,2:nz-1)
      ahfgrady(nx,2:ny-1,2:nz-1) = 2.0D0*ahfgrady(nx-1,2:ny-1,2:nz-1) &
           - ahfgrady(nx-2,2:ny-1,2:nz-1)
      ahfgradz(nx,2:ny-1,2:nz-1) = 2.0D0*ahfgradz(nx-1,2:ny-1,2:nz-1) &
           - ahfgradz(nx-2,2:ny-1,2:nz-1)
      ahfgradn(nx,2:ny-1,2:nz-1) = 2.0D0*ahfgradn(nx-1,2:ny-1,2:nz-1) &
           - ahfgradn(nx-2,2:ny-1,2:nz-1)


!     Boundaries on y direction handled below next section. -Zach


!     Boundaries on z direction.

      ahf_exp(2:nx-1,2:ny-1,1)  = 2.0D0*ahf_exp(2:nx-1,2:ny-1,2) &
           - ahf_exp(2:nx-1,2:ny-1,3)
      ahfgradx(2:nx-1,2:ny-1,1) = 2.0D0*ahfgradx(2:nx-1,2:ny-1,2) &
           - ahfgradx(2:nx-1,2:ny-1,3)
      ahfgrady(2:nx-1,2:ny-1,1) = 2.0D0*ahfgrady(2:nx-1,2:ny-1,2) &
           - ahfgrady(2:nx-1,2:ny-1,3)
      ahfgradz(2:nx-1,2:ny-1,1) = 2.0D0*ahfgradz(2:nx-1,2:ny-1,2) &
           - ahfgradz(2:nx-1,2:ny-1,3)
      ahfgradn(2:nx-1,2:ny-1,1) = 2.0D0*ahfgradn(2:nx-1,2:ny-1,2) &
           - ahfgradn(2:nx-1,2:ny-1,3)

      ahf_exp(2:nx-1,2:ny-1,nz)  = 2.0D0*ahf_exp(2:nx-1,2:ny-1,nz-1) &
           - ahf_exp(2:nx-1,2:ny-1,nz-2)
      ahfgradx(2:nx-1,2:ny-1,nz) = 2.0D0*ahfgradx(2:nx-1,2:ny-1,nz-1) &
           - ahfgradx(2:nx-1,2:ny-1,nz-2)
      ahfgrady(2:nx-1,2:ny-1,nz) = 2.0D0*ahfgrady(2:nx-1,2:ny-1,nz-1) &
           - ahfgrady(2:nx-1,2:ny-1,nz-2)
      ahfgradz(2:nx-1,2:ny-1,nz) = 2.0D0*ahfgradz(2:nx-1,2:ny-1,nz-1) &
           - ahfgradz(2:nx-1,2:ny-1,nz-2)
      ahfgradn(2:nx-1,2:ny-1,nz) = 2.0D0*ahfgradn(2:nx-1,2:ny-1,nz-1) &
           - ahfgradn(2:nx-1,2:ny-1,nz-2)

!     Boundaries on y direction.
      if(nonaxi .and. ny.ne.3) then
         ahf_exp(2:nx-1,ny,2:nz-1)  = 2.0D0*ahf_exp(2:nx-1,ny-1,2:nz-1) &
              - ahf_exp(2:nx-1,ny-2,2:nz-1)
         ahfgradx(2:nx-1,ny,2:nz-1) = 2.0D0*ahfgradx(2:nx-1,ny-1,2:nz-1) &
              - ahfgradx(2:nx-1,ny-2,2:nz-1)
         ahfgrady(2:nx-1,ny,2:nz-1) = 2.0D0*ahfgrady(2:nx-1,ny-1,2:nz-1) &
              - ahfgrady(2:nx-1,ny-2,2:nz-1)
         ahfgradz(2:nx-1,ny,2:nz-1) = 2.0D0*ahfgradz(2:nx-1,ny-1,2:nz-1) &
              - ahfgradz(2:nx-1,ny-2,2:nz-1)
         ahfgradn(2:nx-1,ny,2:nz-1) = 2.0D0*ahfgradn(2:nx-1,ny-1,2:nz-1) &
              - ahfgradn(2:nx-1,ny-2,2:nz-1)
         
         ahf_exp(2:nx-1,1,2:nz-1)  = 2.0D0*ahf_exp(2:nx-1,2,2:nz-1) &
              - ahf_exp(2:nx-1,3,2:nz-1)
         ahfgradx(2:nx-1,1,2:nz-1) = 2.0D0*ahfgradx(2:nx-1,2,2:nz-1) &
              - ahfgradx(2:nx-1,3,2:nz-1)
         ahfgrady(2:nx-1,1,2:nz-1) = 2.0D0*ahfgrady(2:nx-1,2,2:nz-1) &
              - ahfgrady(2:nx-1,3,2:nz-1)
         ahfgradz(2:nx-1,1,2:nz-1) = 2.0D0*ahfgradz(2:nx-1,2,2:nz-1) &
              - ahfgradz(2:nx-1,3,2:nz-1)
         ahfgradn(2:nx-1,1,2:nz-1) = 2.0D0*ahfgradn(2:nx-1,2,2:nz-1) &
              - ahfgradn(2:nx-1,3,2:nz-1)
!     Zach: Sorry for the double negative here... It's needed.
      else 
!     call ezaxibc_scalar(cctk_lsh,X,Y,Z,ahf_exp,excision_zone_gf,Symmetry)
         call BndCartoon2DVN(ierr, cctkGH, 0, "ahfinder::ahf_exp")
         call CCTK_VarIndex(index,'ahfinder::ahfgradx')
         call BndCartoon2DVI(ierr, cctkGH, 1, index)
      end if

    !  write(*,*) "hiiiiiii", 2.0D0*ahf_exp(2,2,2)- ahf_exp(2,3,2)
      do k=1,nz
         do j=1,ny
            do i=1,nx
              ! write(*,*) "AHFEXP",i,j,k,ahf_exp(i,j,k)
            end do
         end do
      end do

!     Synchronize.


      call CCTK_SyncGroup(ierr,cctkGH,"ahfinder::ahfgradient")
      call CCTK_SyncGroup(ierr,cctkGH,"ahfinder::ahfinderexp")

      call CartSymGN(ierr,cctkGH,"ahfinder::ahfgradient")
      call CartSymGN(ierr,cctkGH,"ahfinder::ahfinderexp")


!     ***************
!     ***   END   ***
!     ***************

      return
      end
