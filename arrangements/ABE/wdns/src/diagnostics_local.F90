#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for wdns thorn
!-----------------------------------------------------------------------------
subroutine wdns_diagnostics_local(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                  :: i,j,k
  
  ! If you want to debug, don't run through this every out_every timesteps
!!$  if(MOD(cctk_iteration,out_every)==0) then
!!$     !write(*,*) "COMPUTING B2 INSIDE BHNS THORN, diagnostics_local.F90, dx=",CCTK_DELTA_SPACE(1)
!!$
!!$     !$omp parallel do
!!$     do k = 1,cctk_lsh(3)
!!$        do j = 1,cctk_lsh(2)
!!$           do i = 1,cctk_lsh(1)
!!$              temp7(i,j,k) = 0.D0
!!$              temp8(i,j,k)    = 0.D0
!!$              MONOPOLE(i,j,k) = 0.D0
!!$           end do
!!$        end do
!!$     end do
!!$     !$omp end parallel do
!!$
!!$     !  call bhns_compute_B_from_A(CCTK_PASS_FTOF)
!!$
!!$     call bhns_compute_b2_cpp(cctkGH,cctk_lsh, phi, lapm1, &
!!$          shiftx,shifty,shiftz,vx,vy,vz,Bx,By,Bz, & 
!!$          gxx, gxy, gxz, gyy, gyz, gzz, temp1)
!!$
!!$!!!$omp parallel do
!!$     do k = 1,cctk_lsh(3)
!!$        do j = 1,cctk_lsh(2)
!!$           do i = 1,cctk_lsh(1)
!!$
!!$              if(refbd(i,j,k).gt.0.5D0) then
!!$                 temp8(i,j,k) = -2.D0
!!$              else 
!!$                 temp8(i,j,k) = temp1(i,j,k)
!!$              end if
!!$
!!$              !  do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
!!$              !     do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
!!$              !        do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
!!$              !           if(abs(emask(i,j,k)-1.D0) .lt. 1.D-8) then 
!!$              !temp8(i,j,k)    = temp1(i,j,k)
!!$              !              MONOPOLE(i,j,k) = temp1(i,j,k)/P(i,j,k)
!!$              temp7(i,j,k) = sqrt(Ax(i,j,k)**2+Ay(i,j,k)**2+Az(i,j,k)**2)
!!$
!!$              !if(temp8(i,j,k).gt.0.1 .or. temp8(i,j,k).lt.0.D0) then
!!$              !   write(*,*) "BAD B2... Horizon centroid:",bh_posn_x(1),bh_posn_y(1),bh_posn_z(1)
!!$              !   write(*,*) "BAD B2... Dist from horizon centroid:",sqrt((x(i,j,k)-bh_posn_x(1))**2 + (y(i,j,k)-bh_posn_y(1))**2 + (z(i,j,k)-bh_posn_z(1))**2)
!!$              !   write(*,*) "BAD B2... INSIDE vol_integrand-mag_energies.F90: B2,x,y,z:",temp8(i,j,k),X(i,j,k),Y(i,j,k),Z(i,j,k)
!!$              !end if
!!$
!!$              !           end if
!!$           end do ! i-loop
!!$        end do ! j-loop
!!$     end do ! k-loop
!!$!!!$omp end parallel do
!!$
!!$  end if
end subroutine wdns_diagnostics_local
