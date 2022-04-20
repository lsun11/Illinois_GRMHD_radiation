!---------------------------------------------------------------------
! Sets the Hamiltonian constraint integrand to zero on the ghostzones
!                 *** Excised BHs version ***
!---------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine Ham_constraint_excised_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  integer :: i,j,k
  integer :: AXISYM
  AXISYM=4

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     WhichIntegral = 4

     !Reset VolIntegrand's quickly, using OpenMP
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

     !$omp parallel do
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              VolIntegrand(i,j,k)  = PsiRes(i,j,k) *PsiRes(i,j,k)
              VolIntegrand2(i,j,k) = PsiNorm(i,j,k)*PsiNorm(i,j,k)
           end do
        end do
     end do
     !$omp end parallel do

     if (Symmetry==AXISYM) then
        VolIntegrand  = VolIntegrand *X
        VolIntegrand2 = VolIntegrand2*X
     end if

     if(num_BHs.gt.0) then
        call excise_bhs_VolInt(CCTK_PASS_FTOF)
     end if

     !$omp parallel do
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              if(r(i,j,k).lt.inner_volInt_radius) then
                 VolIntegrand3(i,j,k) = VolIntegrand(i,j,k)
                 VolIntegrand4(i,j,k) = VolIntegrand2(i,j,k)
              end if
           end do
        end do
     end do
     !$omp end parallel do

     if (Symmetry==AXISYM) then
        VolIntegrand3 = VolIntegrand3*X
        VolIntegrand4 = VolIntegrand4*X
     end if
  end if

end subroutine Ham_constraint_excised_integrand
