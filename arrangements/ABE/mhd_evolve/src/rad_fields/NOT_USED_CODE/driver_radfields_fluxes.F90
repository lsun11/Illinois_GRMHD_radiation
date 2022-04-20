!--------------------------------------------------
! Compute B^i fluxes from induction equation: v1.5
!--------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine radfields_fluxes(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS: 
  real*8                :: dX,dY,dZ
  integer               :: pow_axi,i,j,k
  integer, dimension(3) :: ext

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)


!!---------------------construction---------------

if(1==0) then
  call flux_rad_cpp (m, cctkGH, ext, X, Y, Z,&
       tau_rad_flux, tau_rad_flux_diag,& 
       S_radx_flux, S_rady_flux, S_radz_flux, &
       E_radr, E_radl, &
       F_radxr,F_radxl, &
       F_radyr, F_radyl, &
       F_radzr, F_radzl, &
       Pr, Pl,&
       rho_br, rho_bl,&
       Bxr, Bxl, Byr, Byl, Bzr, Bzl,&
       v02r, v02l,&
       vxr,vxl,vyr,vyl,vzr,vzl, &
       gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f, &
       gupxx_f, gupyy_f, gupzz_f,gupxy_f, gupxz_f, gupyz_f, &
       cmax,cmin, &
       shiftx_f, shifty_f, shiftz_f, &
       lapm1_f, phi_f, &
       pow_axi,Symmetry, rad_closure_scheme,&
       enable_OS_collapse,&
       neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th,&
       Erad_atm_cut)
end if


if(excision_enable == 1) then
	call scalar_excision_bc(ext,X,Y,Z,tau_rad_flux,Symmetry,excision_zone_gf);
	call scalar_excision_bc(ext,X,Y,Z,S_radx_flux,Symmetry,excision_zone_gf);
	call scalar_excision_bc(ext,X,Y,Z,S_rady_flux,Symmetry,excision_zone_gf);
	call scalar_excision_bc(ext,X,Y,Z,S_radz_flux,Symmetry,excision_zone_gf);
end if 

end subroutine radfields_fluxes
