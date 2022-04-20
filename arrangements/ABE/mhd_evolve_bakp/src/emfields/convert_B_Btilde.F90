!--------------------------------------------------------------------! 
! Conversion between normal B^i and \tilde{B}^i                      !
!--------------------------------------------------------------------!
!
subroutine convert_b(ext,Bx,By,Bz,phi,b2bt)
  implicit none
  integer, dimension(3)                          :: ext
  real*8, dimension(ext(1),ext(2),ext(3))           :: Bx,By,Bz,phi
  real*8                                         :: sqrtg
  integer                                        :: i,j,k
  real*8					 :: b2bt
  !
  !$omp parallel do
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           sqrtg = exp(6.D0*phi(i,j,k)*b2bt)
           Bx(i,j,k) = Bx(i,j,k) * sqrtg
           By(i,j,k) = By(i,j,k) * sqrtg
           Bz(i,j,k) = Bz(i,j,k) * sqrtg
        end do
     end do
  end do
  !$omp end parallel do

end subroutine convert_b
