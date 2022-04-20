
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_regridder_reader_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: numpoints,checksum
  real*8  :: dx,dy,dz,xmin,ymin,zmin

  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)

  xmin = x(1,1,1)
  ymin = y(1,1,1)
  zmin = z(1,1,1)

  numpoints = cctk_lsh(1)*cctk_lsh(2)*cctk_lsh(3)
  checksum = cctk_lsh(2)*301+cctk_lsh(3)*1+cctk_lsh(1)*1000000 + 1241*CCTK_MyProc(cctkGH)+int(999999*dx)+int(abs(xmin)*10.)+int(abs(zmin)*100.)

  ! Magnetic fields
  call bhns_regridder_in_1gf(numpoints,Ax,checksum)
  call bhns_regridder_in_1gf(numpoints,Ay,checksum)
  call bhns_regridder_in_1gf(numpoints,Az,checksum)
  call bhns_regridder_in_1gf(numpoints,psi6phi,checksum)

  ! Hydro
  call bhns_regridder_in_1gf(numpoints,rho_star,checksum)
  call bhns_regridder_in_1gf(numpoints,tau,checksum)
  call bhns_regridder_in_1gf(numpoints,mhd_st_x,checksum)
  call bhns_regridder_in_1gf(numpoints,mhd_st_y,checksum)
  call bhns_regridder_in_1gf(numpoints,mhd_st_z,checksum)

  ! Metric
  call bhns_regridder_in_1gf(numpoints,phi,checksum)
  call bhns_regridder_in_1gf(numpoints,chi,checksum)
  call bhns_regridder_in_1gf(numpoints,trK,checksum)

  call bhns_regridder_in_1gf(numpoints,gxx,checksum)
  call bhns_regridder_in_1gf(numpoints,gxy,checksum)
  call bhns_regridder_in_1gf(numpoints,gxz,checksum)
  call bhns_regridder_in_1gf(numpoints,gyy,checksum)
  call bhns_regridder_in_1gf(numpoints,gyz,checksum)
  call bhns_regridder_in_1gf(numpoints,gzz,checksum)

  call bhns_regridder_in_1gf(numpoints,Axx,checksum)
  call bhns_regridder_in_1gf(numpoints,Axy,checksum)
  call bhns_regridder_in_1gf(numpoints,Axz,checksum)
  call bhns_regridder_in_1gf(numpoints,Ayy,checksum)
  call bhns_regridder_in_1gf(numpoints,Ayz,checksum)
  call bhns_regridder_in_1gf(numpoints,Azz,checksum)

  call bhns_regridder_in_1gf(numpoints,Gammax,checksum)
  call bhns_regridder_in_1gf(numpoints,Gammay,checksum)
  call bhns_regridder_in_1gf(numpoints,Gammaz,checksum)

  ! Set previous timelevels.
  Ax_p = Ax; Ax_p_p = Ax
  Ay_p = Ay; Ay_p_p = Ay
  Az_p = Az; Az_p_p = Az
  psi6phi_p = psi6phi; psi6phi_p_p = psi6phi; 

  rho_star_p = rho_star; rho_star_p_p = rho_star
  tau_p = tau; tau_p_p = tau
  mhd_st_x_p = mhd_st_x; mhd_st_x_p_p = mhd_st_x
  mhd_st_y_p = mhd_st_y; mhd_st_y_p_p = mhd_st_y
  mhd_st_z_p = mhd_st_z; mhd_st_z_p_p = mhd_st_z

  phi_p = phi; phi_p_p = phi
  chi_p = chi; chi_p_p = chi
  trK_p = trK; trK_p_p = trK

  gxx_p = gxx; gxx_p_p = gxx
  gxy_p = gxy; gxy_p_p = gxy
  gxz_p = gxz; gxz_p_p = gxz
  gyy_p = gyy; gyy_p_p = gyy
  gyz_p = gyz; gyz_p_p = gyz
  gzz_p = gzz; gzz_p_p = gzz

  Axx_p = Axx; Axx_p_p = Axx
  Axy_p = Axy; Axy_p_p = Axy
  Axz_p = Axz; Axz_p_p = Axz
  Ayy_p = Ayy; Ayy_p_p = Ayy
  Ayz_p = Ayz; Ayz_p_p = Ayz
  Azz_p = Azz; Azz_p_p = Azz

  call bhns_compute_B_from_A(CCTK_PASS_FTOF)

end subroutine bhns_regridder_reader_driver
