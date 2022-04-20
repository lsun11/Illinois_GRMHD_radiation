!|----------------------------------------------------+
!| Function to compute fluxes of fluid and EM fields 
!!   across a spherical surface of radius surf_radius.
!!   and centered at the origin
!|----------------------------------------------------+
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mhd_surface_fluxes_2(cctkGH, F_Mesc)
  implicit none
  
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  ! Variables in the function call
  CCTK_POINTER                             :: cctkGH
  real*8                                   :: F_Mesc
  ! Variables needed for interpolation
  real*8, dimension(N_theta*N_phi,3)       :: pointcoords
  integer                                  :: vindex
  ! Output arrays, dummy indices, parameters
  real*8, dimension(N_theta*N_phi)         :: Psi6int,vxint,vyint,vzint,Pint,u0int
  real*8, dimension(N_theta*N_phi)         :: sbxint,sbyint,sbzint,hint,rho_starint,sbtint
  real*8, dimension(N_theta*N_phi)         :: Psi4int,gxxint,gxyint,gxzint,gyyint,gyzint,gzzint
  real*8, dimension(N_theta*N_phi)         :: lapm1int,lapseconst,shiftxint,shiftyint,shiftzint
  real*8                                   :: phiangle,costheta,sintheta,PI
  real*8                                   :: bx_i,by_i,bz_i,bt_i,vx_i,vy_i,vz_i,ut_i,lap_i
  real*8                                   :: beta_x_i,beta_y_i,beta_z_i,betax_i,betay_i,betaz_i
  real*8                                   :: betabeta,betav,Bbeta,u_t_i,b_t_i,bb,u_phi_i,b_phi_i
  real*8                                   :: gxx_i,gxy_i,gxz_i,gyy_i,gyz_i,gzz_i
  real*8                                   :: rhos_i,h_i,P_i
  real*8                                   :: F_E_emx,F_E_emy,F_E_emz
  real*8                                   :: F_J_emx,F_J_emy,F_J_emz
  real*8                                   :: F_J_fx,F_J_fy,F_J_fz
  real*8                                   :: u_x_i,u_y_i,u_z_i,b_x_i,b_y_i,b_z_i
  real*8                                   :: nx_i,ny_i,nz_i,nx_g,ny_g,nz_g,F_M0_i,U_i,vrad_i
  integer                                  :: i,j,n
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  PI = 3.14159265358979323846D0

  !| Psi6int <~~~ exp(6.D0*phi)
  !| Psi4int   <~~~ exp(4.D0*phi),
  !| lapseconst   <~~~ 1.D0/(sqrt(4.0*PI)*(lapm1+1.D0))
  !| lapm1int   <~~~ lapm1,
  !| shiftxint  <~~~ shiftx, shiftyint <~~~ shifty, shiftzint <~~~ shiftz 
  !| gxxint  <~~~ gxx, gxyint <~~~ gxy, gxzint <~~~ gxz, 
  !| gyyint  <~~~ gyy, gyzint <~~~ gyz, gzzint <~~~ gzz, 
  !| sbxint   <~~~ bx,  sbyint  <~~~ by,  sbzint  <~~~ bz,  sbtint  <~~~ bt,
  !| vxint  <~~~ vx, vyint <~~~ vy, vzint <~~~ vz, u0int <~~~ u0
  !| rho_starint  <~~~ rho_star, hint <~~~ h, Pint <~~~ P
 
  n = 1
  do i=1,N_theta
     costheta = 1.0 - (i - 0.5)*dcostheta
     sintheta = sqrt(1.0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5)*dphi
        if(N_phi==1) phiangle = 0.D0
        pointcoords(n,1) = surf_radius*sintheta*cos(phiangle)
        pointcoords(n,2) = surf_radius*sintheta*sin(phiangle)
        pointcoords(n,3) = surf_radius*costheta
        n = n + 1
     end do
  end do

  call CCTK_VarIndex(vindex,"bssn::phi")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,Psi6int)
  call CCTK_VarIndex(vindex,"lapse::lapm1")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,lapm1int)
  do i=1,N_theta*N_phi
     Psi6int(i) = exp(6.D0*Psi6int(i))
     Psi4int(i) = Psi6int(i)**(4.D0/6.D0)
     lapm1int(i) = lapm1int(i) + 1.D0
     lapseconst(i) = 1.D0/(sqrt(4.D0*PI)*(lapm1int(i)))
  end do
  call CCTK_VarIndex(vindex,"shift::shiftx")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,shiftxint)
  call CCTK_VarIndex(vindex,"shift::shifty")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,shiftyint)
  call CCTK_VarIndex(vindex,"shift::shiftz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,shiftzint)

  call CCTK_VarIndex(vindex,"bssn::gxx")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gxxint)
  call CCTK_VarIndex(vindex,"bssn::gxy")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gxyint)
  call CCTK_VarIndex(vindex,"bssn::gxz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gxzint)
  call CCTK_VarIndex(vindex,"bssn::gyy")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gyyint)
  call CCTK_VarIndex(vindex,"bssn::gyz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gyzint)
  call CCTK_VarIndex(vindex,"bssn::gzz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gzzint)

  call CCTK_VarIndex(vindex,"mhd_evolve::sbx")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,sbxint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sby")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,sbyint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,sbzint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbt")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,sbtint)

  do i=1,N_theta*N_phi
     gxxint(i) = gxxint(i) * Psi4int(i)
     gxyint(i) = gxyint(i) * Psi4int(i)
     gxzint(i) = gxzint(i) * Psi4int(i)
     gyyint(i) = gyyint(i) * Psi4int(i)
     gyzint(i) = gyzint(i) * Psi4int(i)
     gzzint(i) = gzzint(i) * Psi4int(i)
     sbxint(i) = sbxint(i) * lapseconst(i)
     sbyint(i) = sbyint(i) * lapseconst(i)
     sbzint(i) = sbzint(i) * lapseconst(i)
     sbtint(i) = sbtint(i) * lapseconst(i)
  end do
  call CCTK_VarIndex(vindex,"mhd_evolve::vx")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,vxint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vy")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,vyint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,vzint)
  call CCTK_VarIndex(vindex,"mhd_evolve::u0")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,u0int)

  call CCTK_VarIndex(vindex,"mhd_evolve::rho_star")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,rho_starint)
  call CCTK_VarIndex(vindex,"mhd_evolve::h")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,hint)
  call CCTK_VarIndex(vindex,"mhd_evolve::P")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,Pint)

  !Zero out our flux sums first.
  F_Mesc   = 0.D0
  !|~~~~~> Integrate the dot product of Dphi with the surface normal
  do i=1,N_theta*N_phi
     bx_i = sbxint(i)
     by_i = sbyint(i)
     bz_i = sbzint(i)
     bt_i = sbtint(i)

     vx_i = vxint(i)
     vy_i = vyint(i)
     vz_i = vzint(i)
     ut_i = u0int(i)

     lap_i = lapm1int(i)

     betax_i = shiftxint(i)
     betay_i = shiftyint(i)
     betaz_i = shiftzint(i)

     gxx_i = gxxint(i)
     gxy_i = gxyint(i)
     gxz_i = gxzint(i)
     gyy_i = gyyint(i)
     gyz_i = gyzint(i)
     gzz_i = gzzint(i)

     rhos_i = rho_starint(i)
     h_i = hint(i)
     P_i = Pint(i)

     nx_g = pointcoords(i,1)/surf_radius
     ny_g = pointcoords(i,2)/surf_radius
     nz_g = pointcoords(i,3)/surf_radius

     nx_i = pointcoords(i,1)
     ny_i = pointcoords(i,2)
     nz_i = pointcoords(i,3)

     vrad_i = nx_g*vx_i + ny_g*vy_i + nz_g*vz_i

     ! \beta_i
     beta_x_i = gxxint(i)*betax_i + gxyint(i)*betay_i + gxzint(i)*betaz_i
     beta_y_i = gxyint(i)*betax_i + gyyint(i)*betay_i + gyzint(i)*betaz_i
     beta_z_i = gxzint(i)*betax_i + gyzint(i)*betay_i + gzzint(i)*betaz_i

     ! \beta.\beta
     betabeta = betax_i*beta_x_i + betay_i*beta_y_i + betaz_i*beta_z_i

     ! \beta.v
     betav = vx_i*beta_x_i + vy_i*beta_y_i + vz_i*beta_z_i

     ! u_t
     u_t_i = u0int(i)*(-lap_i*lap_i + betabeta + betav)
     
     ! U = -1-u_0
     U_i = 1.D0 - u_t_i

     
     F_M0_i    = (nx_g*vx_i + ny_g*vy_i + nz_g*vz_i)*rhos_i
   
     !| Escaping mass
     if ((vrad_i > 0.0d0) .and. (U_i .gt. 0.0d0)) then
        F_Mesc    = F_Mesc + F_M0_i
     end if
  end do

  !|  Multiply by surface element.
  F_Mesc    = F_Mesc* surf_radius * surf_radius * dphi * dcostheta*sym_factor



end subroutine mhd_surface_fluxes_2
