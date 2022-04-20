!--------------------------------------------------------
! Okay, we've read in the initial data from files.
! Now we set up all other required variables, including:
!  emfields, BSSN variables, primitives, etc.
!--------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine magnetar_setup_remaining_matter_quantities(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: delta_bar

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  integer :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)

  write(*,*) "part3: set up remaining matter quantitites"

  ext = cctk_lsh

  n_poly = Pr(1,1,1)
  RADEQUAT = Pr(1,1,2)

  !delta_bar sets up a bar-mode, I believe
  delta_bar = 0.D0 

  !
  ! Compute the initial matter profile
  !
  if(enable_alt_atmosphere==0) then
     call compute_magnetar_hybrid(ext,P_deplete, RADEQUAT, delta_bar, &
          X,Y,Z, &
          neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
          lapm1,phi,shiftx,shifty, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          rho,S,Sx,Sy,Sz, &
          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
          rho_star,tau,st_x,st_y,st_z, &
          P,w,vx,vy,vz,rho_b, &
          u0,h,rho_b_atm,Omega_Frame,PhysicalRadius,eps_flag,K_poly,n_poly)
  else
     write(*,*) "BEFORE COMP ALT ATM"
     call compute_magnetar_hybrid_alt_atmosphere(ext,P_deplete, RADEQUAT, delta_bar, &
          X,Y,Z, &
          neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
          lapm1,phi,shiftx,shifty, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          rho,S,Sx,Sy,Sz, &
          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
          rho_star,tau,st_x,st_y,st_z, &
          P,w,vx,vy,vz,rho_b, &
          u0,h,rho_b_atm_gf,pfloor_gf,Omega_Frame,PhysicalRadius,eps_flag,K_poly,n_poly)
     write(*,*) "AFTER COMP ALT ATM"
  end if

  write(*,*) "Pressua",P(1,1,1)

end subroutine magnetar_setup_remaining_matter_quantities
