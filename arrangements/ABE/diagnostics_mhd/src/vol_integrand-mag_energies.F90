!-----------------------------------------------------------
! Set up integrand for the EM field energy: (b^2/2) u^0 \sqrt{-g}
!-----------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mag_energies_integrand(CCTK_ARGUMENTS)
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
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)

  PI = 3.14159265358979323846D0

  WhichIntegral = 109

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     !-----------------------------------------------------------------------------
     ! Compute integrand
     !-----------------------------------------------------------------------------
     !Reset VolIntegrand, VolIntegrand2, VolIntegrand3, VolIntegrand4 quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              VolIntegrand(i,j,k) = 0.D0
              VolIntegrand2(i,j,k) = 0.D0
              VolIntegrand3(i,j,k) = 0.D0
              VolIntegrand4(i,j,k) = 0.D0
!              temp8(i,j,k) = 0.D0
!              MONOPOLE(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do

     call compute_mag_energy_dens_cpp(cctkGH,cctk_lsh,phi,lapm1,shiftx,shifty, &
                shiftz,vx,vy,vz,Bx,By,Bz,gxx,gxy,gxz,gyy,gyz,gzz,VolIntegrand,VolIntegrand3,temp1)

     write(*,*) "INSIDE mag_energies_integrand.  SHOULD BE BEFORE IOASCII JUNK!"

     if (num_BHs .ne. 0) then
        !$omp parallel do
        do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
           do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
              do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
                 if(abs(emask(i,j,k)-1.D0) .lt. 1.D-8) then 
!                    temp8(i,j,k) = temp1(i,j,k)
!                    MONOPOLE(i,j,k) = temp1(i,j,k)/P(i,j,k)
!
!                    if(temp8(i,j,k).gt.0.1 .or. temp8(i,j,k).lt.0.D0) then
!                       write(*,*) "BAD B2... INSIDE vol_integrand-mag_energies.F90",temp8(i,j,k)
!                    end if

                    VolIntegrand2(i,j,k) = VolIntegrand(i,j,k)
                    VolIntegrand4(i,j,k) = VolIntegrand3(i,j,k)
                 end if
              end do ! i-loop
           end do ! j-loop
        end do ! k-loop
        !$omp end parallel do
     end if
  end if

  return
end subroutine mag_energies_integrand
