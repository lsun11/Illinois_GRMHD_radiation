subroutine komar_mass_integrand(ex,integrand,X,phi,alpha,rho,S, &
        betax,betay,betaz, Sx,Sy,Sz, Symmetry)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X
  real*8, dimension(ex(1),ex(2),ex(3))        :: phi,alpha,rho,S
  real*8, dimension(ex(1),ex(2),ex(3))        :: betax,betay,betaz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Sx,Sy,Sz
  integer                                     :: Symmetry
! Output oarameter
!
  real*8, dimension(ex(1),ex(2),ex(3))        :: integrand
! 
! Other parameters
!
  integer                            :: AXISYM
  integer                            :: i,j,k
  parameter (AXISYM = 4)
!
! Set up integrand
!
  do k=1,ex(3)
     do j=1,ex(2)
        do i=1,ex(1)
           !write(*,*) hi.,i,j,k,phi(i,j,k),alpha(i,j,k),rho(i,j,k),S(i,j,k),betax(i,j,k),Sx(i,j,k)
           integrand(i,j,k) = exp(6.d0*phi(i,j,k)) * ( (1.d0+alpha(i,j,k))* &
                (rho(i,j,k)+S(i,j,k)) - 2.d0 * (betax(i,j,k)*Sx(i,j,k) &
                + betay(i,j,k)*Sy(i,j,k) + betaz(i,j,k)*Sz(i,j,k) ) )
           if (Symmetry==AXISYM) integrand(i,j,k) = X(i,1,1)*integrand(i,j,k)
           if(integrand(i,j,k).ne.0.D0) write(*,*) "hi.",i,j,k,integrand(i,j,k),phi(i,j,k),alpha(i,j,k),rho(i,j,k),S(i,j,k),betax(i&
  &,j,k),Sx(i,j,k)
        end do
     end do
  end do
end subroutine komar_mass_integrand
