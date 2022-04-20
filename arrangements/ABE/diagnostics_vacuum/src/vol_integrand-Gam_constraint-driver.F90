!--------------------------------------------
! Driver for the gamma constraint integrands
!--------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine gam_constraint_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8 :: dX,dY,dZ
  integer :: i,j,k
  integer :: AXISYM
  AXISYM=4

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     WhichIntegral = 7

     !Reset VolIntegrand's quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              VolIntegrand(i,j,k) = 0.D0
              VolIntegrand2(i,j,k) = 0.D0
              VolIntegrand3(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do

     call gamcheck(cctkGH,cctk_lsh, cctk_nghostzones, Symmetry, &
          dx,  dy,  dz,  &
          gconx,gcony,gconz, &
          Gammax, Gammay, Gammaz,  &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)

     !$omp parallel do
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              if (Symmetry==AXISYM) then
                 VolIntegrand(i,j,k)  = X(i,j,k)*gconx(i,j,k)*gconx(i,j,k)
                 VolIntegrand2(i,j,k) = X(i,j,k)*gcony(i,j,k)*gcony(i,j,k)
                 VolIntegrand3(i,j,k) = X(i,j,k)*gconz(i,j,k)*gconz(i,j,k)
              else
                 VolIntegrand(i,j,k)  = gconx(i,j,k)*gconx(i,j,k)
                 VolIntegrand2(i,j,k) = gcony(i,j,k)*gcony(i,j,k)
                 VolIntegrand3(i,j,k) = gconz(i,j,k)*gconz(i,j,k)
              end if
           end do
        end do
     end do
     !$omp end parallel do

     if(num_BHs.gt.0) then
        call excise_bhs_VolInt(CCTK_PASS_FTOF)
     end if

  end if
end subroutine gam_constraint_integrand
