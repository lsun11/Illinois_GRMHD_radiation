#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine OS_initialdata_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8 		:: dX,dY,dZ
  integer               :: i,j,k
  character             :: varname*30
  real*8                :: pgxx,pgxy,pgxz,pgyy,pgyz,pgzz
  real*8                :: alpha,psi2
  real*8                :: v,gamma,Radius,rho0
   REAL*8, PARAMETER :: PI_D=3.141592653589793238462643383279502884197

 
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh
  write(*,*) "cctk_lsh(1): ",cctk_lsh(1)
  write(*,*) "cctk_lsh(2): ",cctk_lsh(2)
  write(*,*) "cctk_lsh(3): ",cctk_lsh(3)

  rho0=1.d0/(4.d0/3.0*PI_D*R_edge**3)
 
    do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
	   psi2=psi(i,j,k)*psi(i,j,k)
	   alpha = lapm1(i,j,k)+1.d0
	   Radius=sqrt(X(i,j,k)*X(i,j,k) &
			+ Y(i,j,k)*Y(i,j,k) &
			+ Z(i,j,k)*Z(i,j,k))
	   !set rho_b
	   if (Radius*psi2 .lt. R_edge) then
             rho_b(i,j,k)=rho0
	   else
             rho_b(i,j,k)=0.d0
	   endif
	   vx(i,j,k) = 0.d0
	   vy(i,j,k) = 0.d0
	   vz(i,j,k) = 0.d0
	   u0(i,j,k) = 1.d0/alpha
        end do
     end do
  end do
  
 ! !set K_poly
 ! K_poly=1.d-6
  write(*,*) "K_poly: ",K_poly
  !eos stuff
  n_poly=1.0d0/(gamma_th-1.0d0)
  neos=1
  rho_tab(1)=0.03
  P_tab(1) = k_poly * rho_tab(1)**gamma_th
  eps_tab(1)=n_poly*P_tab(1)/rho_tab(1)
  do i=1,2
     k_tab(i)=K_poly
     gamma_tab(i)=gamma_th
  enddo

end subroutine OS_initialdata_driver
