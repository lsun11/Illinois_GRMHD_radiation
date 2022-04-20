#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "GenericFD.h"

subroutine PsiKadelia_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i,j,k,ierr
  integer,dimension(3) :: ones
  real*8  :: dX,dY,dZ,psi4th
! 1st of 2 needed #includes for GenericFD.h:
#include "../../../GenFD_decl_varF90.h"

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !     Initialise finite differencing variables
  ! 2nd of 2 needed #includes for GenericFD.h:
#include "../../../GenFD_set_varF90.h"

  if(MOD(cctk_iteration,out_every).eq.0) then
     !  if(MOD(cctk_iteration,4).eq.0) then

     write(*,*) "howdy...",dX,cctk_lsh

     !if(CCTK_DELTA_SPACE(1)==cctk_delta_space(1)) then
     !     if(cctk_iteration==0) then
     !        !$omp parallel do
     !        do k=1,cctk_lsh(3)
     !           do j=1,cctk_lsh(2)
     !              do i=1,cctk_lsh(1)
     !                 psi0re(i,j,k) = 0.D0
     !                 psi0im(i,j,k) = 0.D0
     !              end do
     !           end do
     !        end do
     !        !$omp end parallel do
     !     end if

     ! Note: we use PsiNorm,PsiRes,PsiTau, MRsx, MRsy, and MRsz for temporary storage!
     !$omp parallel do private(psi4th)
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              Psi(i,j,k) = exp(phi(i,j,k))
              psi4th = Psi(i,j,k)*Psi(i,j,k)*Psi(i,j,k)*Psi(i,j,k)
              Kxx(i,j,k) = psi4th * (Axx(i,j,k) + (1.D0/3.D0) * gxx(i,j,k) * trK(i,j,k))
              Kxy(i,j,k) = psi4th * (Axy(i,j,k) + (1.D0/3.D0) * gxy(i,j,k) * trK(i,j,k))
              Kxz(i,j,k) = psi4th * (Axz(i,j,k) + (1.D0/3.D0) * gxz(i,j,k) * trK(i,j,k))
              Kyy(i,j,k) = psi4th * (Ayy(i,j,k) + (1.D0/3.D0) * gyy(i,j,k) * trK(i,j,k))
              Kyz(i,j,k) = psi4th * (Ayz(i,j,k) + (1.D0/3.D0) * gyz(i,j,k) * trK(i,j,k))
              Kzz(i,j,k) = psi4th * (Azz(i,j,k) + (1.D0/3.D0) * gzz(i,j,k) * trK(i,j,k))

              !We use  gxxx, gxxy, gxxz, gxyx, gxyy, gxyz as temporary storage...
              gxxx(i,j,k) = gxx(i,j,k)*psi4th
              gxxy(i,j,k) = gxy(i,j,k)*psi4th
              gxxz(i,j,k) = gxz(i,j,k)*psi4th
              gxyx(i,j,k) = gyy(i,j,k)*psi4th
              gxyy(i,j,k) = gyz(i,j,k)*psi4th
              gxyz(i,j,k) = gzz(i,j,k)*psi4th
           end do
        end do
     end do
     !$omp end parallel do

     call PsiKadelia_lowlevel(cctkGH,cctk_lsh,cctk_nghostzones, &
          dX,dY,dZ,r,x,y,z, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gxxx,gxxy,gxxz,gxyx,gxyy,gxyz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          phi,psi, &
          Kxx,Kxy,Kxz,Kyy,Kyz,Kzz, & 
          Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, & 
          Gammaxxx,Gammayxx,Gammazxx,Gammaxyy,Gammayyy,Gammazyy, &
          Gammaxzz,Gammayzz,Gammazzz,Gammaxxy,Gammayxy,Gammazxy, &
          Gammaxxz,Gammayxz,Gammazxz,Gammaxyz,Gammayyz,Gammazyz, &
          psi0re,psi0im)

     ones = 1

     !     call BndFlatGN(ierr,cctkGH,ones,"gw_extraction::Psi4_gw")
     !     call CCTK_SyncGroup(ierr,cctkGH,"gw_extraction::Psi4_gw")
     !     write(*,*) "ierrb",ierr
     call CartSymGN(ierr,cctkGH,"gw_extraction::Psi4_gw")

     ! Finally, we restore the gij,k's used as temporary storage above:
     !$omp parallel do
     do k=cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j=cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i=cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              psi0re_p(i,j,k) = psi0re(i,j,k)
              psi0re_p_p(i,j,k) = psi0re(i,j,k)
              psi0im_p(i,j,k) = psi0im(i,j,k)
              psi0im_p_p(i,j,k) = psi0im(i,j,k)

              gxxx(i,j,k) = D1gf(gxx,i,j,k)
              gxxy(i,j,k) = D2gf(gxx,i,j,k)
              gxxz(i,j,k) = D3gf(gxx,i,j,k)
              gxyx(i,j,k) = D1gf(gxy,i,j,k)
              gxyy(i,j,k) = D2gf(gxy,i,j,k)
              gxyz(i,j,k) = D3gf(gxy,i,j,k)
           end do
        end do
     end do
     !$omp end parallel do

!!$     ! Next, we fill in the other two timelevels.
!!$     !$omp parallel do
!!$     do k=1,cctk_lsh(3)
!!$        do j=1,cctk_lsh(2)
!!$           do i=1,cctk_lsh(1)
!!$              ! BAD IDEA:
!!$              !psi0re_p_p(i,j,k) = psi0re(i,j,k)
!!$              !psi0im_p_p(i,j,k) = psi0im(i,j,k)
!!$              !              psi0re_p(i,j,k) = psi0re(i,j,k)
!!$              !              psi0im_p(i,j,k) = psi0im(i,j,k)
!!$
!!$              if(k.lt.4.and.j.gt.4.and.i.gt.4.and.isnan(psi0re(i,j,k)).and.CCTK_DELTA_SPACE(1).ne.cctk_delta_space(1)) then
!!$!                 write(*,*) "found nan at",i,j,k,CCTK_DELTA_SPACE(1)
!!$!                 stop
!!$              end if
!!$           end do
!!$        end do
!!$     end do
!!$     !$omp end parallel do
  end if

end subroutine PsiKadelia_driver
