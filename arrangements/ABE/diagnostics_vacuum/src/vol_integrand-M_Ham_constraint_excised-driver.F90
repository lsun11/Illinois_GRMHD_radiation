!---------------------------------------------------------------------
! Sets the Hamiltonian constraint integrand to zero on the ghostzones
!                 *** Excised BHs version ***
!---------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine M_Ham_constraint_excised_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  integer :: i,j,k
  integer :: AXISYM
  real*8, parameter :: pi = 3.14159265358979323846d0
  real*8 :: f1o2pi


  AXISYM=4

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0 .and. enable_M_constraint==1) then
     WhichIntegral = 8
     f1o2pi = 0.5d0/pi

     !Reset VolIntegrand's quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              VolIntegrand(i,j,k) = 0.D0
              VolIntegrand2(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do

     !$omp parallel do
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
	      ! The integrand is rho_con psi^5, which is the contribution 
	      ! of ADM mass due to matter source term. 
   	      ! Here rho_con = |H|/(16pi) = |PsiRes|/(2pi), 
              ! H = R - K_ij K^ij +K^2 - 16 pi rho
              VolIntegrand(i,j,k)  = abs(PsiRes(i,j,k))*f1o2pi
           end do
        end do
     end do
     !$omp end parallel do

     if (Symmetry==AXISYM) then
        VolIntegrand  = VolIntegrand *X
     end if

     if(num_BHs.gt.0) then
        call excise_bhs_VolInt(CCTK_PASS_FTOF)
     end if

     !$omp parallel do
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              !!if(r(i,j,k).lt.inner_volInt_radius) then
              !!   VolIntegrand2(i,j,k) = VolIntegrand(i,j,k)
              !!end if
	      if ( (x(i,j,k)-bh_posn_x(1))**2+(y(i,j,k)-bh_posn_y(1))**2+(z(i,j,k)-bh_posn_z(1))**2 & 
		   .lt. inner_volInt_radius**2) VolIntegrand2(i,j,k) = VolIntegrand(i,j,k)
	      if (num_BHs.gt.1 .and. ( (x(i,j,k)-bh_posn_x(2))**2+(y(i,j,k)-bh_posn_y(2))**2+ & 
                     (z(i,j,k)-bh_posn_z(2))**2 .lt. inner_volInt_radius**2) ) & 
                  VolIntegrand2(i,j,k) = VolIntegrand(i,j,k)
           end do
        end do
     end do
     !$omp end parallel do

     if (Symmetry==AXISYM) then
        VolIntegrand2 = VolIntegrand2*X
     end if
  end if

end subroutine M_Ham_constraint_excised_integrand
