!----------------------------
! Set up rest mass integrand
!----------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine rest_mass_inside_AH_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: dX, dY, dZ
  integer                            :: i,j,k
  integer                            :: AXISYM
  parameter(AXISYM = 4)

  WhichIntegral = 107

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

     if (num_BHs .ne. 0) then
        !$omp parallel do
        do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
           do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
              do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
                 if(abs(emask(i,j,k)-1.D0) .ge. 1.D-8) then
                    VolIntegrand(i,j,k) = rho_star(i,j,k) 
                 end if
              end do ! i-loop
           end do ! j-loop
        end do ! k-loop
        !$omp end parallel do
        if (Symmetry == AXISYM) then
           VolIntegrand = X*VolIntegrand
        end if
     end if

!     if(b2_gridfunction(i,j,k).gt.0.1) then
!        write(*,*) "BAD B2... INSIDE vol_integrand-rest_mass_inside_AH.F90",b2_gridfunction(i,j,k)
!     end if

  end if

  !Need to reset mask.  Otherwise AH finder will basically turn off.
  !Reset emask quickly, using OpenMP
  if (num_BHs .ne. 0) then
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              emask(i,j,k) = 1.D0
           end do
        end do
     end do
     !$omp end parallel do
  end if
end subroutine rest_mass_inside_AH_integrand
