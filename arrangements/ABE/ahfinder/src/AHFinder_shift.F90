/*@@
  @file      AHFinder_shift.F
  @date      May 2000
  @author    Miguel Alcubierre
  @desc 
             Set up horizon (area) locking shift.  The area
             locking shift is obtained from the following
             expression:

                 a        __a
             beta   =  b  \/ f


             where f is the "horizon function", and where b is given by:


                         /             ab               2 \
             b  =  alpha | trK   -    K      d f d f / u  |
                         \                    a   b       /

                         / __2      __a__b              2 \
                     /   | \/ f  -  \/ \/ f  d f d f / u  |
                         \                    a   b       /


             with "alpha" the lapse function and:


              2     mn
             u  =  g   d f d f
                        m   n

  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_DefineThorn.h"

      subroutine AHFinder_shift(CCTK_ARGUMENTS)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k
      integer order_save

      CCTK_REAL det,idet
      CCTK_REAL gdxx,gdyy,gdzz,gdxy,gdxz,gdyz
      CCTK_REAL guxx,guyy,guzz,guxy,guxz,guyz
      CCTK_REAL kdxx,kdyy,kdzz,kdxy,kdxz,kdyz
      CCTK_REAL ddxxx,ddxyy,ddxzz,ddxxy,ddxxz,ddxyz
      CCTK_REAL ddyxx,ddyyy,ddyzz,ddyxy,ddyxz,ddyyz
      CCTK_REAL ddzxx,ddzyy,ddzzz,ddzxy,ddzxz,ddzyz
      CCTK_REAL dfx,dfy,dfz,dfux,dfuy,dfuz
      CCTK_REAL d2fxx,d2fyy,d2fzz,d2fxy,d2fxz,d2fyz
      CCTK_REAL c2fxx,c2fyy,c2fzz,c2fxy,c2fxz,c2fyz
      CCTK_REAL idx,idy,idz,idx2,idy2,idz2,idxy,idxz,idyz
      CCTK_REAL zero,half,one,two
      CCTK_REAL aux,T0,T1,T2,T3,T4,T5

      CCTK_REAL, DIMENSION(nx,ny,nz) :: mask,betax_old,betay_old,betaz_old


!     *******************
!     ***   NUMBERS   ***
!     *******************

      zero = 0.0D0
      half = 0.5D0
      one  = 1.0D0
      two  = 2.0D0


!     *********************************
!     ***   SAY WHAT WE ARE DOING   ***
!     *********************************

      write(*,*)
      write(*,*) 'AHFinder: Setting up horizon locking shift.'
      write(*,*)

      betax_old = betax
      betay_old = betay
      betaz_old = betaz


!     *******************************
!     ***   SET UP SHIFT VECTOR   ***
!     *******************************

      idx = half/dx
      idy = half/dy
      idz = half/dz

      idxy = idx*idy
      idxz = idx*idz
      idyz = idy*idz

      idx2 = one/dx**2
      idy2 = one/dy**2
      idz2 = one/dz**2

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

!              Find determinant of spatial metric.

               det = gdxx*gdyy*gdzz + two*gdxy*gdxz*gdyz
     .             - gdxx*gdyz**2 - gdyy*gdxz**2 - gdzz*gdxy**2

               idet = one/det

!              Find inverse spatial metric.

               guxx = idet*(gdyy*gdzz - gdyz**2)
               guyy = idet*(gdxx*gdzz - gdxz**2)
               guzz = idet*(gdxx*gdyy - gdxy**2)

               guxy = idet*(gdxz*gdyz - gdzz*gdxy)
               guxz = idet*(gdxy*gdyz - gdyy*gdxz)
               guyz = idet*(gdxy*gdxz - gdxx*gdyz)

!              Find first spatial derivatives of f.

               dfx = idx*(ahfgrid(i+1,j,k) - ahfgrid(i-1,j,k))
               dfy = idy*(ahfgrid(i,j+1,k) - ahfgrid(i,j-1,k))
               dfz = idz*(ahfgrid(i,j,k+1) - ahfgrid(i,j,k-1))

!              Raise indices in first derivatives.

               dfux = guxx*dfx + guxy*dfy + guxz*dfz
               dfuy = guxy*dfx + guyy*dfy + guyz*dfz
               dfuz = guxz*dfx + guyz*dfy + guzz*dfz

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

!              Find second spatial derivatives of f.

               d2fxx = (ahfgrid(i+1,j,k) - two*ahfgrid(i,j,k)
     .               +  ahfgrid(i-1,j,k))*idx2
               d2fyy = (ahfgrid(i,j+1,k) - two*ahfgrid(i,j,k)
     .               +  ahfgrid(i,j-1,k))*idy2
               d2fzz = (ahfgrid(i,j,k+1) - two*ahfgrid(i,j,k)
     .               +  ahfgrid(i,j,k-1))*idz2

               d2fxy = (ahfgrid(i+1,j+1,k) + ahfgrid(i-1,j-1,k)
     .               -  ahfgrid(i+1,j-1,k) - ahfgrid(i-1,j+1,k))*idxy
               d2fxz = (ahfgrid(i+1,j,k+1) + ahfgrid(i-1,j,k-1)
     .               -  ahfgrid(i+1,j,k-1) - ahfgrid(i-1,j,k+1))*idxz
               d2fyz = (ahfgrid(i,j+1,k+1) + ahfgrid(i,j-1,k-1)
     .               -  ahfgrid(i,j+1,k-1) - ahfgrid(i,j-1,k+1))*idyz

!              Find second covariant derivatives of f.

               c2fxx = d2fxx - half*(dfux*ddxxx
     .               + dfuy*(two*ddxxy - ddyxx)
     .               + dfuz*(two*ddxxz - ddzxx))

               c2fyy = d2fyy - half*(dfuy*ddyyy
     .               + dfux*(two*ddyxy - ddxyy)
     .               + dfuz*(two*ddyyz - ddzyy))

               c2fzz = d2fzz - half*(dfuz*ddzzz
     .               + dfux*(two*ddzxz - ddxzz)
     .               + dfuy*(two*ddzyz - ddyzz))

               c2fxy = d2fxy - half*(dfux*ddyxx + dfuy*ddxyy
     .               + dfuz*(ddxyz + ddyxz - ddzxy))

               c2fxz = d2fxz - half*(dfux*ddzxx + dfuz*ddxzz
     .               + dfuy*(ddxyz + ddzxy - ddyxz))

               c2fyz = d2fyz - half*(dfuy*ddzyy + dfuz*ddyzz
     .               + dfux*(ddyxz + ddzxy - ddxyz))

!                                 2    /  m      \ - 1
!              Find:   T0 =  1 / u  =  | d f d f |
!                                      \      m  /

               T0 = dfx*dfux + dfy*dfuy + dfz*dfuz

               if (T0.gt.zero) then
                  T0 = one/T0
               else
                  betax(i,j,k) = zero
                  betay(i,j,k) = zero
                  betaz(i,j,k) = zero
                  goto 10
               end if

!                            __2
!              Find:  T1  =  \/ f
 
               T1 = guxx*c2fxx + guyy*c2fyy + guzz*c2fzz
     .            + two*(guxy*c2fxy + guxz*c2fxz + guyz*c2fyz)

!                            __ __     a   b     2
!              Find:  T2  =  \/ \/ f  d f d f / u
!                              a  b 

               T2 = c2fxx*dfux**2 + c2fyy*dfuy**2 + c2fzz*dfuz**2
     .            + two*(c2fxy*dfux*dfuy + c2fxz*dfux*dfuz
     .            + c2fyz*dfuy*dfuz)

               T2 = T2*T0

!              Find:  T3  =  trK

               T3 = guxx*kdxx + guyy*kdyy + guzz*kdzz
     .            + two*(guxy*kdxy + guxz*kdxz + guyz*kdyz)

!                                 a   b     2
!              Find:  T4  =  K   d f d f / u
!                             ab

               T4 = kdxx*dfux**2 + kdyy*dfuy**2 + kdzz*dfuz**2
     .            + two*(kdxy*dfux*dfuy + kdxz*dfux*dfuz
     .            + kdyz*dfuy*dfuz)

               T4 = T4*T0

!              Find shift.

               aux = alp(i,j,k)*(T3 - T4)/(T1 - T2)

               betax(i,j,k) = aux*dfux
               betay(i,j,k) = aux*dfuy
               betaz(i,j,k) = aux*dfuz

 10            continue

            end do
         end do
      end do

!     Boundaries.

      if (CCTK_EQUALS(ahf_shiftbound,'extrap')) then

!        X direction.

         betax(1,:,:) = two*betax(2,:,:) - betax(3,:,:)
         betay(1,:,:) = two*betay(2,:,:) - betay(3,:,:)
         betaz(1,:,:) = two*betaz(2,:,:) - betaz(3,:,:)

         betax(nx,:,:) = two*betax(nx-1,:,:) - betax(nx-2,:,:)
         betay(nx,:,:) = two*betay(nx-1,:,:) - betay(nx-2,:,:)
         betaz(nx,:,:) = two*betaz(nx-1,:,:) - betaz(nx-2,:,:)

!        Y direction.

         betax(:,1,:) = two*betax(:,2,:) - betax(:,3,:)
         betay(:,1,:) = two*betay(:,2,:) - betay(:,3,:)
         betaz(:,1,:) = two*betaz(:,2,:) - betaz(:,3,:)

         betax(:,ny,:) = two*betax(:,ny-1,:) - betax(:,ny-2,:)
         betay(:,ny,:) = two*betay(:,ny-1,:) - betay(:,ny-2,:)
         betaz(:,ny,:) = two*betaz(:,ny-1,:) - betaz(:,ny-2,:)

!        Z direction.

         betax(:,:,1) = two*betax(:,:,2) - betax(:,:,3)
         betay(:,:,1) = two*betay(:,:,2) - betay(:,:,3)
         betaz(:,:,1) = two*betaz(:,:,2) - betaz(:,:,3)

         betax(:,:,nz) = two*betax(:,:,nz-1) - betax(:,:,nz-2)
         betay(:,:,nz) = two*betay(:,:,nz-1) - betay(:,:,nz-2)
         betaz(:,:,nz) = two*betaz(:,:,nz-1) - betaz(:,:,nz-2)

      end if


!     ********************
!     ***   EXCISION   ***
!     ********************

      if (ahf_shiftexcise.eq.1) then

#ifdef EXCISION_SIMPLEEXCISION

         order_save = order
         order = 1

         call SimpleExcision(_CCTK_ARGUMENTS,x,y,z,r,mask,
     .      betax,betax_old,zero,one,'extrapolation')
         call SimpleExcision(_CCTK_ARGUMENTS,x,y,z,r,mask,
     .      betay,betay_old,zero,one,'extrapolation')
         call SimpleExcision(_CCTK_ARGUMENTS,x,y,z,r,mask,
     .      betaz,betaz_old,zero,one,'extrapolation')

         order = order_save

#else
         call CCTK_WARN(0,"You have not compiled with SimpleExcision")
#endif

      end if


!     ********************************************
!     ***   APPLY SYMMETRIES AND SYNCHRONIZE   ***
!     ********************************************

      call CartSymGN(ierr,cctkGH,"shift::shift_evolve")
      call CCTK_SyncGroup(ierr,cctkGH,"shift::shift_evolve")


!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder_shift

