
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_regridder_reader_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: numpoints,checksum,i,j,k
  real*8  :: dx,dy,dz,xmin,ymin,zmin

  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)

  xmin = x(1,1,1)
  ymin = y(1,1,1)
  zmin = z(1,1,1)

  numpoints = cctk_lsh(1)*cctk_lsh(2)*cctk_lsh(3)
  write(*,*) "PARAMS:",cctk_lsh,CCTK_MyProc(cctkGH),dx,xmin,ymin,zmin
  checksum = cctk_lsh(2)*301+cctk_lsh(3)*1+cctk_lsh(1)*1000000 + 1241*CCTK_MyProc(cctkGH)+int(999999*dx)+int(abs(xmin)*10.)+int(abs(zmin)*100.)

  ! Read in Magnetic fields, Conservatives, & Metric
  call bhns_regridder_in_allgfs(numpoints,checksum, &
       Ax,Ay,Az,psi6phi, &
       rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z, &
       shiftx,shifty,shiftz,lapm1, &
       phi,chi,trK, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       Axx,Axy,Axz,Ayy,Ayz,Azz, &
       Gammax,Gammay,Gammaz, &
       rho_b_atm,tau_atm)
  
!!$  do k=1,cctk_lsh(3)
!!$     do j=1,cctk_lsh(2)
!!$        do i=1,cctk_lsh(1)
!!$           if(abs(y(i,j,k)).lt.dy*0.01 .and. abs(z(i,j,k)).lt.dz*0.01) then
!!$              if(abs(Ax(i,j,k) - 5*x(i,j,k)**2).gt.dx*0.01) then
!!$                 write(*,*) "BAD",x(i,j,k),Ax(i,j,k),5*x(i,j,k)**2
!!$                 stop
!!$              end if
!!$           end if
!!$        end do
!!$     end do
!!$  end do
  
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

  st_x = mhd_st_x
  st_y = mhd_st_y
  st_z = mhd_st_z
  
  shiftx_p = shiftx; shiftx_p_p = shiftx
  shifty_p = shifty; shifty_p_p = shifty
  shiftz_p = shiftz; shiftz_p_p = shiftz
  lapm1_p = lapm1; lapm1_p_p = lapm1

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

  Gammax_p = Gammax; Gammax_p_p = Gammax
  Gammay_p = Gammay; Gammay_p_p = Gammay
  Gammaz_p = Gammaz; Gammaz_p_p = Gammaz

  gamma_th = 2.D0
  ! K_poly=1 because Keisuke's datafiles are in polytropic units already!!!
  k_poly=1.0d0
  n_poly=1.0d0/(gamma_th-1.0d0)

  neos=1
  rho_tab(1)=1.0
  P_tab(1)=K_poly
  eps_tab(1)=K_poly/(gamma_th-1.0d0)
  do i=1,2
     k_tab(i)=K_poly
     gamma_tab(i)=gamma_th
  enddo

  ! First compute gupij, based on gij:
  temp1 = 1.D0/(gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz &
       - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz)
  
  gupxx =   ( gyy * gzz - gyz * gyz )* temp1
  gupxy = - ( gxy * gzz - gyz * gxz )* temp1
  gupxz =   ( gxy * gyz - gyy * gxz )* temp1
  gupyy =   ( gxx * gzz - gxz * gxz )* temp1
  gupyz = - ( gxx * gyz - gxy * gxz )* temp1
  gupzz =   ( gxx * gyy - gxy * gxy )* temp1

  !Finally, compute K_ij's for apparent horizon finder:
  Psi = exp(phi)
  Kxx = Psi*Psi*Psi*Psi * (Axx + (1.D0/3.D0) * gxx * trK)
  Kxy = Psi*Psi*Psi*Psi * (Axy + (1.D0/3.D0) * gxy * trK)
  Kxz = Psi*Psi*Psi*Psi * (Axz + (1.D0/3.D0) * gxz * trK)
  Kyy = Psi*Psi*Psi*Psi * (Ayy + (1.D0/3.D0) * gyy * trK)
  Kyz = Psi*Psi*Psi*Psi * (Ayz + (1.D0/3.D0) * gyz * trK)
  Kzz = Psi*Psi*Psi*Psi * (Azz + (1.D0/3.D0) * gzz * trK)

  !====================================
end subroutine bhns_regridder_reader_driver
