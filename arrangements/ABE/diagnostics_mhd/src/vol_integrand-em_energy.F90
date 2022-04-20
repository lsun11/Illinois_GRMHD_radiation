!-----------------------------------------------------------
! Set up integrand for the EM field energy: (E^2 + B^2)/8pi
!-----------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine em_energy_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: AXISYM

  real*8                             :: xcenter,ycenter,zcenter,horizdirn_x,horizdirn_y,horizdirn_z,distfrombhcenter,horiz_radius
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)

  PI = 3.14159265358979323846D0

  WhichIntegral = 103

!  write(*,*) "Start vol_integrand-em_energy.F90"
  
  if(cctk_iteration==19456 .and. 1==0) then
     write(*,*) "RESETTING A FIELDS BY FACTOR OF 1/6.43!"
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              Ax(i,j,k) = 1.D0/6.43*Ax(i,j,k)
              Ay(i,j,k) = 1.D0/6.43*Ay(i,j,k)
              Az(i,j,k) = 1.D0/6.43*Az(i,j,k)

              Ax_p(i,j,k) = 1.D0/6.43*Ax_p(i,j,k)
              Ay_p(i,j,k) = 1.D0/6.43*Ay_p(i,j,k)
              Az_p(i,j,k) = 1.D0/6.43*Az_p(i,j,k)

              Ax_p_p(i,j,k) = 1.D0/6.43*Ax_p_p(i,j,k)
              Ay_p_p(i,j,k) = 1.D0/6.43*Ay_p_p(i,j,k)
              Az_p_p(i,j,k) = 1.D0/6.43*Az_p_p(i,j,k)
           end do
        end do
     end do
     !$omp end parallel do
  end if

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     !-----------------------------------------------------------------------------
     ! Compute integrand
     !-----------------------------------------------------------------------------
     !Reset VolIntegrand and VolIntegrand2 quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              VolIntegrand(i,j,k) = 0.D0
              VolIntegrand2(i,j,k) = 0.D0
              VolIntegrand3(i,j,k) = 0.D0
              VolIntegrand4(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do

     if (num_BHs .ne. 0) then
        xcenter = bh_posn_x(1)
        ycenter = bh_posn_y(1)
        zcenter = bh_posn_z(1)
        
        !     Get horizon radius in direction where it is likely to be minimized:
        horizdirn_x = 0.D0
        horizdirn_y = 0.D0
        horizdirn_z = 100000.D0
        call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,horiz_radius);
     end if


     !$omp parallel do
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              if (num_BHs .ne. 0) then
                 VolIntegrand(i,j,k) = exp(10.d0*phi(i,j,k)) * &
                      ( gxx(i,j,k)*(Ex(i,j,k)**2 + Bx(i,j,k)**2) + &
                      2.d0*gxy(i,j,k)*(Ex(i,j,k)*Ey(i,j,k) + Bx(i,j,k)*By(i,j,k)) + &
                      2.d0*gxz(i,j,k)*(Ex(i,j,k)*Ez(i,j,k) + Bx(i,j,k)*Bz(i,j,k)) + &
                      gyy(i,j,k)*(Ey(i,j,k)**2 + By(i,j,k)**2) + &
                      2.d0*gyz(i,j,k)*(Ey(i,j,k)*Ez(i,j,k) + By(i,j,k)*Bz(i,j,k)) + &
                      gzz(i,j,k)*(Ez(i,j,k)**2 + Bz(i,j,k)**2) ) / (8.D0*PI)
                 if(abs(emask(i,j,k)-1.D0) .lt. 1.D-8) then
                    VolIntegrand2(i,j,k) = VolIntegrand(i,j,k)
                 end if

                 distfrombhcenter= sqrt((x(i,j,k)-xcenter)**2 + (y(i,j,k)-ycenter)**2 + (z(i,j,k)-zcenter)**2)
                 if(distfrombhcenter.gt.1.7D0) then
                    VolIntegrand3(i,j,k) = VolIntegrand(i,j,k)
                 end if
                 if(distfrombhcenter.gt.3.4D0) then
                    VolIntegrand4(i,j,k) = VolIntegrand(i,j,k)
                 end if
              else
                 VolIntegrand(i,j,k) = exp(10.d0*phi(i,j,k)) * &
                      ( gxx(i,j,k)*(Ex(i,j,k)**2 + Bx(i,j,k)**2) + &
                      2.d0*gxy(i,j,k)*(Ex(i,j,k)*Ey(i,j,k) + Bx(i,j,k)*By(i,j,k)) + &
                      2.d0*gxz(i,j,k)*(Ex(i,j,k)*Ez(i,j,k) + Bx(i,j,k)*Bz(i,j,k)) + &
                      gyy(i,j,k)*(Ey(i,j,k)**2 + By(i,j,k)**2) + &
                      2.d0*gyz(i,j,k)*(Ey(i,j,k)*Ez(i,j,k) + By(i,j,k)*Bz(i,j,k)) + &
                      gzz(i,j,k)*(Ez(i,j,k)**2 + Bz(i,j,k)**2) ) / (8.D0*PI)
              end if
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
     !$omp end parallel do
     if (Symmetry == AXISYM) then
        VolIntegrand = X*VolIntegrand
     end if
  end if

  return
end subroutine em_energy_integrand
