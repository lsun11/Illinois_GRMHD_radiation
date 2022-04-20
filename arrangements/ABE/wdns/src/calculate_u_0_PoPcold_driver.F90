
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine calc_u_0_PoPcold_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8 		:: dX,dY,dZ
  real*8                :: gamma1, gamma2, gamma3
  real*8                :: kappa1, kappa2, kappa3
  real*8                :: rhoo1, rhoo2
!!  real*8                :: an1, an2, an3, c2
!!  integer               :: i,ierr



  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  rhoo1=rho_tab(1)
  rhoo2=rho_tab(2)
  gamma1 = gamma_tab(1)
  gamma2 = gamma_tab(2)
  gamma3 = gamma_tab(3)
  kappa1 = k_tab(1)
  kappa2 = k_tab(2)
  kappa3 = k_tab(3)
  
  write(*,*) "rhoo1=", rhoo1,   "rhoo2=", rhoo2   
  write(*,*) "gamma1=", gamma1, "gamma2=", gamma2, "gamma3=", gamma3
  write(*,*) "kappa1=", kappa1, "kappa2=", kappa2, "kappa3=", kappa3


!!  write(*,*) "n_ijk:",cctk_lsh
!!  write(*,*) "ZZZ:",Z(1,1,:)
!!  write(*,*) "X:",X(:,1,1)


     call calc_u_0_PoPcold(cctkGH, X(1,1,1),Y(1,1,1),Z(1,1,1), &
          dX,dY,dZ,cctk_lsh, &
          lapm1,shiftx,shifty,shiftz, &
          psi,gxx,gxy,gxz,gyy,gyz, gzz, &
          rho_b, u0, P, vx,vy,vz, temp8, temp9, gamma1, gamma2, gamma3, kappa1, kappa2, kappa3, rhoo1, rhoo2)



end subroutine calc_u_0_PoPcold_driver
