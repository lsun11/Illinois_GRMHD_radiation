!-----------------------------------------------------------------------------
! Compute \rho_star \epsilon to get the total internal energy integrand 
!  This version is hybrid EOS - compatible
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine minternal_hybrid_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> ABE Functions used:
  interface
     subroutine compute_pcold_epscold(rhob, P_cold, eps_cold, &
          neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse)
       implicit none
       integer :: neos,ergo_star, enable_OS_collapse
       real*8  :: rhob, P_cold, eps_cold, ergo_sigma
       real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
       real*8, dimension(neos+1) :: k_tab, gamma_tab
     end subroutine compute_pcold_epscold
  end interface

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: dV,cs,sn,pi,FOUR,ONE, rho_tiny
  real*8                             :: rhob, P_cold, eps_cold, eps
  integer                            :: i,j,k
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)

  WhichIntegral = 102
  
  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then

     rho_tiny = rho_b_atm

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
!$omp parallel do private (rhob,eps,P_cold,eps_cold)
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)

              rhob = rho_b(i,j,k)
              if (rhob .gt. rho_tiny) then
                 call compute_pcold_epscold(rhob, P_cold, eps_cold, &
                      neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse)
                 if (enable_OS_collapse.eq.1) then
                    eps = eps_cold;
                 else
                 eps = eps_cold + (P(i,j,k)-P_cold)/(gamma_th-1.d0)/rhob
                 end if	
              else
                 eps = 0.d0
              end if

              if (Symmetry == AXISYM) then
                 VolIntegrand(i,j,k) = rho_star(i,j,k)*eps* abs(X(i,1,1))
              else
                 VolIntegrand(i,j,k) = rho_star(i,j,k)*eps
              end if
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
!$omp end parallel do
  end if
end subroutine minternal_hybrid_integrand
