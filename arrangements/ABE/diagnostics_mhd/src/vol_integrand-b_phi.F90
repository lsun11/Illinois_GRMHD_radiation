!---------------------------------------------------------
! Integrand for integral of |B_{\phi}| over proper volume
!---------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine b_phi_integrand(CCTK_ARGUMENTS)
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
  real*8                             :: Bsubphi

  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)

  WhichIntegral = 105

!  write(*,*) "Start vol_integrand-b_phi"
  
  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     !-----------------------------------------------------------------------------
     ! Compute integrand
     !-----------------------------------------------------------------------------
     !Reset VolIntegrand quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              VolIntegrand(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do

     !!!!$omp parallel do
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              Bsubphi = X(i,1,1)*exp(4.d0*phi(i,j,k))*(gxy(i,j,k)*Bx(i,j,k) + &
                   gyy(i,j,k)*By(i,j,k) + gyz(i,j,k)*Bz(i,j,k))
              VolIntegrand(i,j,k) = rho_star(i,j,k)*abs(Bsubphi)
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
     !!!!$omp end parallel do
     if (Symmetry == AXISYM) then
        VolIntegrand = X*VolIntegrand
     end if
  end if
end subroutine b_phi_integrand
