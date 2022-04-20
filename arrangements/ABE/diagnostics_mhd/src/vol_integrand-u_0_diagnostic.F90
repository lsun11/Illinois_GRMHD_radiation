!----------------------------
! Set up rest mass integrand
!----------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine u_0_diagnostic_inside_AH(CCTK_ARGUMENTS)
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

  WhichIntegral = 1070

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

              ! u_0 = u^0 * (-alpha^2 + gamma_ij beta^i (beta^j + v^j))
              ! -1 - u_0 > 0 --> unbound, otherwise bound.
              VolIntegrand2(i,j,k) = -1.D0 - u0(i,j,k)*(-(lapm1(i,j,k)+1.D0)**2 + &
                   exp(4.D0*phi(i,j,k)) * ( &
                   gxx(i,j,k)*shiftx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) + &
                   2.D0*gxy(i,j,k)*shiftx(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                   2.D0*gxz(i,j,k)*shiftx(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) + &
                   gyy(i,j,k)*shifty(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                   2.D0*gyz(i,j,k)*shifty(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) + &
                   gzz(i,j,k)*shiftz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k))))

              VolIntegrand3(i,j,k) = 0.D0
              VolIntegrand4(i,j,k) = 0.D0

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
                    if(VolIntegrand2(i,j,k).gt.0.D0) then
                       VolIntegrand3(i,j,k) = rho_star(i,j,k)
                    else
                       VolIntegrand4(i,j,k) = rho_star(i,j,k)
                    end if
                 end if
              end do ! i-loop
           end do ! j-loop
        end do ! k-loop
        !$omp end parallel do
        if (Symmetry == AXISYM) then
           VolIntegrand = X*VolIntegrand
        end if
     end if

  end if

end subroutine u_0_diagnostic_inside_AH
