!-----------------------------------------------------------------------------
! Integrand of T (total kinetic energy) = Omega S_phi
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine kinetic_energy_T_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,R_axi,cs,sn,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - one )

  WhichIntegral = 100
  
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
     ! Remove OpenMP support until we remove the temporary variables
     !!!$omp parallel do
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              if (Symmetry == AXISYM) then
                 VolIntegrand(i,j,k) = 0.5D0 * vy(i,j,k) * st_y(i,j,k) * abs(X(i,1,1))
              else
                 ! ADDED 1.D-10 in following line to avoid NaN's in cs, sn when R_axi==0
                 R_axi = sqrt(X(i,1,1)*X(i,1,1) + Y(1,j,1)*Y(1,j,1)) + 1.D-10
                 cs = X(i,1,1)/R_axi
                 sn = Y(1,j,1)/R_axi
                 VolIntegrand(i,j,k) = 0.5D0 * (cs*vy(i,j,k) - sn*vx(i,j,k)) &
                      * (cs*st_y(i,j,k) - sn*st_x(i,j,k)) 
              end if
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
     !!!$omp end parallel do
  end if

end subroutine kinetic_energy_T_integrand
