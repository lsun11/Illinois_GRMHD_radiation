
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_ridX(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8 		:: dX,dY,dZ, f1odX, f1odY, f1odZ, b2bt
  integer               :: i,j,k,dummy

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)
  f1odX = 1.D0/dX
  f1odY = 1.D0/dY
  f1odZ = 1.D0/dZ

  ext = cctk_lsh

  if(isnan(Bx(1,1,1)).eqv..false.) then

     write(*,*) "============================================="
     write(*,*) "INSIDE RIDX"
     write(*,*) "============================================="

     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservatives')

     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              !Note that the primitives solver requires Bx,By,Bz as input.
              Bx(i,j,k) = Bxtilde(i,j,k)
              By(i,j,k) = Bytilde(i,j,k)
              Bz(i,j,k) = Bztilde(i,j,k)
           end do
        end do
     end do
     !$omp end parallel do

     ! Notice that we set B^i = tildeB^i above.
     b2bt = -1.D0
     call convert_b(ext,Bx,By,Bz,phi,b2bt)

     !$omp parallel do
     do i=1,cctk_lsh(1)
        do j=1,cctk_lsh(2)
           do k=1,cctk_lsh(3)
              temp11(i,j,k) = exp(6.D0*phi(i,j,k))
              temp7(i,j,k) = temp11(i,j,k)*Bx(i,j,k)
              temp8(i,j,k) = temp11(i,j,k)*By(i,j,k)
              temp9(i,j,k) = temp11(i,j,k)*Bz(i,j,k)
           end do
        end do
     end do
     !$omp end parallel do

     !First, compute the divergence of B everywhere, except stencil ghostzones:
     !$omp parallel do
     do k = 2,cctk_lsh(3)
        do j = 2,cctk_lsh(2)
           do i = 2,cctk_lsh(1)
              temp10(i,j,k) = &
                   0.5D0*(0.5D0*(f1odX*( &
                   (temp7(i,j,k) - temp7(i-1,j,k)) + &
                   (temp7(i,j-1,k) - temp7(i-1,j-1,k)) + &
                   (temp7(i,j,k-1) - temp7(i-1,j,k-1)) + &
                   (temp7(i,j-1,k-1) - temp7(i-1,j-1,k-1)) &
                   ))) + &
                   0.5D0*(0.5D0*(f1odY*( &
                   (temp8(i,j,k) - temp8(i,j-1,k)) + &
                   (temp8(i-1,j,k) - temp8(i-1,j-1,k)) + &
                   (temp8(i,j,k-1) - temp8(i,j-1,k-1)) + &
                   (temp8(i-1,j,k-1) - temp8(i-1,j-1,k-1)) &
                   ))) + &
                   0.5D0*(0.5D0*(f1odZ*( &
                   (temp9(i,j,k) - temp9(i,j,k-1)) + &
                   (temp9(i,j-1,k) - temp9(i,j-1,k-1)) + &
                   (temp9(i-1,j,k) - temp9(i-1,j,k-1)) + &
                   (temp9(i-1,j-1,k) - temp9(i-1,j-1,k-1)) &
                   )))
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
     !$omp end parallel do

     call ridX_lowlevel(cctkGH,ext,dX,dY,dZ, &
          temp7,temp8,temp9,temp10)

     !$omp parallel do
     do i=1,cctk_lsh(1)
        do j=1,cctk_lsh(2)
           do k=1,cctk_lsh(3)
              Bx(i,j,k) = temp7(i,j,k)/temp11(i,j,k)
              By(i,j,k) = temp8(i,j,k)/temp11(i,j,k)
              Bz(i,j,k) = temp9(i,j,k)/temp11(i,j,k)
              Bxtilde(i,j,k) = Bx(i,j,k)
              Bytilde(i,j,k) = By(i,j,k)
              Bztilde(i,j,k) = Bz(i,j,k)
           end do
        end do
     end do
     !$omp end parallel do

     ! Here, we convert B^i to tilde B^i
     b2bt = 1.D0
     call convert_b(ext,Bxtilde,Bytilde,Bztilde,phi,b2bt)

     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservatives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')

     write(*,*) "============================================="
     write(*,*) "LEAVING RIDX"
     write(*,*) "============================================="
  end if

end subroutine driver_ridX
