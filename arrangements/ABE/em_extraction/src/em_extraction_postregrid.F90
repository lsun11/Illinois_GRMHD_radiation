#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine em_extraction_postregrid(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i,j,k,ierr
  integer,dimension(3) :: ones
  real*8  :: dX,dY,dZ

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

!!$  if(MOD(cctk_iteration,out_every).eq.0) then
     
     !write(*,*) "howdy...in em_extraction postregrid",dX,cctk_lsh

!!$
!!$     NPphi0re_p = 1.d100
!!$     NPphi0im_p = 1.d100
!!$     NPphi2re_p = 1.d100
!!$     NPphi2im_p = 1.d100
!!$
!!$     NPphi0re_p_p = 1.d100
!!$     NPphi0im_p_p = 1.d100
!!$     NPphi2re_p_p = 1.d100
!!$     NPphi2im_p_p = 1.d100

     
     call em_extraction_lowlevel(cctkGH,cctk_lsh,cctk_nghostzones, &
          dX,dY,dZ,r,x,y,z, &
          gxx,gxy,gxz,gyy,gyz,gzz,lapm1,shiftx,shifty,shiftz, &
          phi,psi, &
          Bx,By,Bz,Ex,Ey,Ez, &
          NPphi0re,NPphi0im,NPphi2re,NPphi2im)
     
     ones = 1
     
     ! Vasilis says: I 'm not so sure about the following, but it's probably OK
     
     call CartSymGN(ierr,cctkGH,"em_extraction::NPphi2_em")
     
     ! Finally, we set the previous level stuff
!!$     do k=cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
!!$        do j=cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
!!$           do i=cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
!!$     !$omp parallel do
!!$     do k=1,cctk_lsh(3)
!!$        do j=1,cctk_lsh(2)
!!$           do i=1,cctk_lsh(1)
!!$              phi0re_p(i,j,k) = phi0re(i,j,k)
!!$              phi0re_p_p(i,j,k) = phi0re(i,j,k)
!!$              phi0im_p(i,j,k) = phi0im(i,j,k)
!!$              phi0im_p_p(i,j,k) = phi0im(i,j,k)
!!$              
!!$              phi2re_p(i,j,k) = phi2re(i,j,k)
!!$              phi2re_p_p(i,j,k) = phi2re(i,j,k)
!!$              phi2im_p(i,j,k) = phi2im(i,j,k)
!!$              phi2im_p_p(i,j,k) = phi2im(i,j,k)
!!$           end do
!!$        end do
!!$     end do
!!$     !$omp end parallel do         
!!$
!!$  end if

end subroutine em_extraction_postregrid
