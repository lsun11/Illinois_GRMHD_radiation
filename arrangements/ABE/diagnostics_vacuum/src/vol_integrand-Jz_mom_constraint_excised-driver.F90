!-------------------------------------------------------------------
! Sets the momentum constraint integrands to zero on the ghostzones
!                 *** Excised BHs version ***
!-------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine Jz_mom_constraint_excised_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  integer :: i,j,k
  integer, parameter :: AXISYM = 4
  real*8, parameter :: pi = 3.14159265358979323846d0
  real*8 :: f1o8pi

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     WhichIntegral = 9
     f1o8pi = 0.125d0/pi

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
              ! Note that MRsi is exp(-6 phi) M^i, where M^i is the
              ! quantity defined in Eq. (17) of gr-qc/0209102
              VolIntegrand(i,j,k)  = f1o8pi*exp(10.d0*phi(i,j,k)) * abs( X(i,j,k) * ( & 
		gxy(i,j,k)*MRsx(i,j,k) + gyy(i,j,k)*MRsy(i,j,k) + & 
		gyz(i,j,k)*MRsz(i,j,k) ) - Y(i,j,k) * ( & 
                gxx(i,j,k)*MRsx(i,j,k) + gxy(i,j,k)*MRsy(i,j,k) + &
                gxz(i,j,k)*MRsz(i,j,k) ) )
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
              if(r(i,j,k).lt.inner_volInt_radius) then
                 VolIntegrand2(i,j,k) = VolIntegrand(i,j,k)
              end if
           end do
        end do
     end do
     !$omp end parallel do

     if (Symmetry==AXISYM) then
        VolIntegrand2 = VolIntegrand2*X
     end if

  end if

end subroutine Jz_mom_constraint_excised_integrand
