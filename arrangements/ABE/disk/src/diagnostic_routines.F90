#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! This version of fluxes uses an analytic determination of the event
! horizon radius.  For use with Cowling runs--such as disks around black
! holes a la McKinney and Gammie. 

! The strange thing about this is that the output is placed in
! steerable parameters instead of grid variables.
subroutine fluxes_ana(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: dphi, dcostheta, sym_factor, ONE
  real*8 :: F_M0_ijk, F_E_fluid_ijk, F_E_em_ijk, F_J_fluid_ijk, F_J_em_ijk
  real*8 :: psi6_i, psi4_i, lap_i, fact, btx_i, bty_i, btz_i
  real*8 :: gxx_i, gxy_i, gxz_i, gyy_i, gyz_i, gzz_i
  real*8 :: bt_i, bx_i, by_i, bz_i, ut_i, vx_i, vy_i, vz_i
  real*8 :: rhos_i, h_i, P_i, bt_x_i, bt_y_i, bt_z_i
  real*8 :: bb, bv, Bbeta, u_t_i, b_t_i, u_x_i, u_y_i, u_phi_i
  real*8 :: b_x_i, b_y_i, b_phi_i
  real*8 :: F_E_emx, F_E_emy, F_E_emz, F_J_fx, F_J_fy, F_J_fz
  real*8 :: F_J_emx, F_J_emy, F_J_emz, rr_i, dS, PI, dX, dY, dZ, f1ospi
  real*8 :: F_E_emxt, F_E_emyt, F_E_emzt
  integer :: i, interp_order
  integer, dimension(3)                     :: ext,global_ext
  real*8, dimension(n_tot,3)                :: pointcoords
  real*8, dimension(n_tot)                  :: phiint, utint, vxint, vyint, vzint
  real*8, dimension(n_tot)                  :: sbtint, sbxint, sbyint, sbzint
  real*8, dimension(n_tot)                  :: gintxx,gintxy,gintxz,gintyy,gintyz,gintzz
  real*8, dimension(n_tot)                  :: lapseint,shiftintx,shiftinty,shiftintz 
  real*8, dimension(n_tot)                  :: rhosint, hint, Pint
  ! testing
  real*8, dimension(n_tot)                  :: Bxint, Byint, Bzint
  real*8                                    :: lbx_i, lby_i, lbz_i


  PI = 3.14159265358979323844D0
  dphi = 2.0 * PI / N_phi
  dcostheta = 1.0 / N_theta
  interp_order = 1
  ONE = 1.0
  ext = cctk_lsh
  global_ext = cctk_gsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  do i=1,n_tot
     pointcoords(i,1) = xn_d(i)
     pointcoords(i,2) = yn_d(i)
     pointcoords(i,3) = zn_d(i)
  end do

  ! call interpolation routine.
  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,phi,phiint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,lapm1,lapseint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,u0,utint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,vx,vxint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,vy,vyint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,vz,vzint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,sbt,sbtint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,sbx,sbxint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,sby,sbyint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,sbz,sbzint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gxx,gintxx) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gxy,gintxy) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gxz,gintxz) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gyy,gintyy) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gyz,gintyz) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gzz,gintzz) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,shiftx,shiftintx) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,shifty,shiftinty) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,shiftz,shiftintz) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,rho_star,rhosint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,h,hint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,P,Pint) 

  ! testing
  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Bx,Bxint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,By,Byint) 

  call interpolate_pointset(cctkGH,cctk_nghostzones,n_tot,interp_order,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Bz,Bzint) 

  !
  ! integrate on the horizon
  !

  F_M0 = 0.d0
  F_E_fluid = 0.d0
  F_E_em = 0.d0
  F_J_fluid = 0.d0
  F_J_em = 0.d0

  F_E_emxt = 0.d0
  F_E_emyt = 0.d0
  F_E_emzt = 0.d0

  f1ospi = 1.d0/sqrt(4.d0*PI)
  do i=1,n_tot 
     !  do i=3,n_tot 
     ! read in data
     psi6_i = exp(6.d0*phiint(i))
     psi4_i = exp(4.d0*phiint(i))
     lap_i = lapseint(i) + 1.0
     fact = f1ospi/lap_i
     btx_i = shiftintx(i)
     bty_i = shiftinty(i)
     btz_i = shiftintz(i)
     gxx_i = gintxx(i)*psi4_i
     gxy_i = gintxy(i)*psi4_i
     gxz_i = gintxz(i)*psi4_i
     gyy_i = gintyy(i)*psi4_i
     gyz_i = gintyz(i)*psi4_i
     gzz_i = gintzz(i)*psi4_i
     bt_i  = sbtint(i)*fact
     bx_i  = sbxint(i)*fact
     by_i  = sbyint(i)*fact
     bz_i  = sbzint(i)*fact
     ut_i  = utint(i)
     vx_i  = vxint(i)
     vy_i  = vyint(i)
     vz_i  = vzint(i)
     rhos_i = rhosint(i)
     h_i  = hint(i)
     P_i  = Pint(i)

     ! testing
     lbx_i = Bxint(i)
     lby_i = Byint(i)
     lbz_i = Bzint(i)

     ! \beta_i
     bt_x_i = gxx_i*btx_i + gxy_i*bty_i + gxz_i*btz_i
     bt_y_i = gxy_i*btx_i + gyy_i*bty_i + gyz_i*btz_i
     bt_z_i = gxz_i*btx_i + gyz_i*bty_i + gzz_i*btz_i
     ! \beta.\beta
     bb = btx_i*bt_x_i + bty_i*bt_y_i + btz_i*bt_z_i
     ! \beta.v
     bv = vx_i*bt_x_i + vy_i*bt_y_i + vz_i*bt_z_i
     ! \beta.b
     Bbeta = bx_i*bt_x_i + by_i*bt_y_i + bz_i*bt_z_i

     ! u_t
     u_t_i = ut_i*(-lap_i*lap_i + bb + bv)

     ! b_t
     b_t_i = bt_i*(-lap_i*lap_i + bb) + Bbeta

     ! u_x, u_y, and u_phi
     u_x_i = ut_i*( gxx_i*(vx_i + btx_i) + gxy_i*(vy_i + bty_i) &
          + gxz_i*(vz_i + btz_i) )
     u_y_i = ut_i*( gxy_i*(vx_i + btx_i) + gyy_i*(vy_i + bty_i) &
          + gyz_i*(vz_i + btz_i) )
     u_phi_i = xn_d(i)*u_y_i - yn_d(i)*u_x_i

     ! b_x, b_y, and b_phi 
     b_x_i = bt_x_i*bt_i + gxx_i*bx_i + gxy_i*by_i + gxz_i*bz_i
     b_y_i = bt_y_i*bt_i + gxy_i*bx_i + gyy_i*by_i + gyz_i*bz_i
     b_phi_i = xn_d(i)*b_y_i - yn_d(i)*b_x_i

     ! b.b
     bb = -lap_i*lap_i*bt_i*bt_i &
          + gxx_i*(bx_i+btx_i*bt_i)*(bx_i+btx_i*bt_i) &
          + gyy_i*(by_i+bty_i*bt_i)*(by_i+bty_i*bt_i) &
          + gzz_i*(bz_i+btz_i*bt_i)*(bz_i+btz_i*bt_i) &
          + 2.0*( gxy_i*(bx_i+btx_i*bt_i)*(by_i+bty_i*bt_i) &
          + gxz_i*(bx_i+btx_i*bt_i)*(bz_i+btz_i*bt_i) &
          + gyz_i*(by_i+bty_i*bt_i)*(bz_i+btz_i*bt_i) )

     F_E_emx = bb*ut_i*vx_i*u_t_i - bx_i*b_t_i
     F_E_emy = bb*ut_i*vy_i*u_t_i - by_i*b_t_i
     F_E_emz = bb*ut_i*vz_i*u_t_i - bz_i*b_t_i

     F_J_fx = rhos_i*h_i*vx_i*u_phi_i - lap_i*psi6_i*P_i*yn_d(i)
     F_J_fy = rhos_i*h_i*vy_i*u_phi_i + lap_i*psi6_i*P_i*xn_d(i)
     F_J_fz = rhos_i*h_i*vz_i*u_phi_i

     F_J_emx = (ut_i*vx_i*u_phi_i - 0.5*yn_d(i))*bb - bx_i*b_phi_i
     F_J_emy = (ut_i*vy_i*u_phi_i + 0.5*xn_d(i))*bb - by_i*b_phi_i
     F_J_emz = ut_i*vz_i*u_phi_i*bb - bz_i*b_phi_i

     F_M0_ijk = (nx_d(i)*vx_i + ny_d(i)*vy_i + nz_d(i)*vz_i)*rhos_i
     F_E_fluid_ijk = -F_M0_ijk*h_i*u_t_i
     F_E_em_ijk = -(nx_d(i)*F_E_emx + ny_d(i)*F_E_emy + nz_d(i)*F_E_emz)*lap_i*psi6_i
     F_J_fluid_ijk = (nx_d(i)*F_J_fx + ny_d(i)*F_J_fy + nz_d(i)*F_J_fz)
     F_J_em_ijk = (nx_d(i)*F_J_emx + ny_d(i)*F_J_emy + nz_d(i)*F_J_emz)*lap_i*psi6_i

     rr_i = xn_d(i)*xn_d(i) + yn_d(i)*yn_d(i) + zn_d(i)*zn_d(i)

     if (Symmetry==AXISYM) then
        dS = rr_i * (4.0*PI) * dcostheta
     else
        dS = rr_i * dphi * dcostheta
     end if

     F_M0      = F_M0      + F_M0_ijk      * dS
     F_E_fluid = F_E_fluid + F_E_fluid_ijk * dS
     F_E_em    = F_E_em    + F_E_em_ijk    * dS
     F_J_fluid = F_J_fluid + F_J_fluid_ijk * dS
     F_J_em    = F_J_em    + F_J_em_ijk    * dS

     F_E_emxt = F_E_emxt -nx_d(i)*F_E_emx*lap_i*psi6_i*dS
     F_E_emyt = F_E_emyt -ny_d(i)*F_E_emy*lap_i*psi6_i*dS
     F_E_emzt = F_E_emzt -nz_d(i)*F_E_emz*lap_i*psi6_i*dS

     !      write(*,100) bb, ut_i, vz_i, u_t_i, bz_i, b_t_i
     ! Uncomment below line to enable velocity output
     !      write(*,100) ut_i, vx_i, vy_i, vz_i

  end do

  write(*,*) 'leaving fluxes_ana'
!!$   write(*,*) 'F_E_emxt = ', F_E_emxt
!!$   write(*,*) 'F_E_emyt = ', F_E_emyt
  write(*,*) 'F_E_emzt = ', F_E_emzt
100 FORMAT (6e12.3)
end subroutine fluxes_ana

subroutine bi_phys(ext,X,Y,Z,PhysR,RDerv,Bx,By,Bz,Bx_phys,By_phys,Bz_phys,Symmetry)
  implicit none
  ! Input parameters:
  integer, dimension(3)                       :: ext
  real*8, dimension(ext(1),ext(2),ext(3))     :: X,Y,Z,PhysR,RDerv
  integer                                     :: Symmetry
  real*8                                      :: rhob_cut
  real*8, dimension(ext(1),ext(2),ext(3))     :: Bx,By,Bz,Bx_phys,By_phys,Bz_phys

  ! Other variables 
  real*8                                      :: Jxx,Jxy,Jxz,Jyy,Jyz,Jzz
  real*8                                      :: Rp,dR,r,xi,yi,zi,fac
  real*8                                      :: Bxi,Byi,Bzi
  integer                                     :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  integer, parameter                          :: AXISYM = 4
  !
  imin = lbound(Bx,1)
  imax = ubound(Bx,1)
  kmin = lbound(Bx,3)
  kmax = ubound(Bx,3) 
  if (Symmetry==AXISYM) then
     jmin = 2
     jmax = 2
  else 
     jmin = lbound(Bx,2)
     jmax = ubound(Bx,2)
  end if
  do k=kmin,kmax
     do j=jmin,jmax
        do i=imin,imax
           xi = X(i,j,k)
           zi = Z(i,j,k)
           yi = Y(i,j,k)
           r = sqrt(xi*xi + yi*yi + zi*zi)
           Bxi = Bx(i,j,k)
           Byi = By(i,j,k)
           Bzi = Bz(i,j,k)
           Rp = PhysR(i,j,k)
           dR = RDerv(i,j,k)
           fac = dR - Rp/r
           Jxx = Rp/r + xi/r * xi/r * fac
           Jxy = xi/r * yi/r * fac
           Jxz = xi/r * zi/r * fac
           Jyy = Rp/r + yi/r * yi/r * fac
           Jyz = yi/r * zi/r * fac
           Jzz = Rp/r + zi/r * zi/r * fac
           Bx_phys(i,j,k) = Jxx*Bxi + Jxy*Byi + Jxz*Bzi
           By_phys(i,j,k) = Jxy*Bxi + Jyy*Byi + Jyz*Bzi
           Bz_phys(i,j,k) = Jxz*Bxi + Jyz*Byi + Jzz*Bzi
        end do
     end do
  end do

  if (Symmetry==AXISYM) then
     ! rotate those guys?
     call axibc_vector(ext,X,Y,Z,Bx_phys,By_phys,Bz_phys)
  end if

end subroutine bi_phys
!-----------------------------------------------------------------------------
!
! Integrand for the EM field energy: (E^2 + B^2)/8pi
!
!-----------------------------------------------------------------------------
!
subroutine em_energy_disk(ext, demenergy, &
     X, Y, Z, phi, gxx,gxy,gxz,gyy,gyz,gzz, Ex,Ey,Ez, &
     Bx,By,Bz, Symmetry, mskf, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
  ! Input parameters:
  integer, dimension(3)                        :: ext
  real*8, dimension(ext(1),ext(2),ext(3))      :: X,Y,Z, mskf
  real*8, dimension(ext(1),ext(2),ext(3))      :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))      :: Ex,Ey,Ez, Bx,By,Bz
  integer                                      :: adjimin, adjjmin, adjkmin
  integer                                      :: adjimax, adjjmax, adjkmax
  integer                                      :: Symmetry,disk_flag
  ! output:
  real*8                                       :: demenergy
  ! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE, temp
  integer                            :: i,j,k, ilow, klow
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - ONE )
  !
  ! Coordinate grid size, putting adjustments on such that
  ! it excludes the ghost zones.
  !
  imin = lbound(phi,1) - adjimin
  jmin = lbound(phi,2) - adjjmin
  kmin = lbound(phi,3) - adjkmin
  imax = ubound(phi,1) - adjimax
  jmax = ubound(phi,2) - adjjmax
  kmax = ubound(phi,3) - adjkmax
  !
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  end if

  ilow = imin
  klow = kmin

  !-----------------------------------------------------------------------------
  ! integrate
  !-----------------------------------------------------------------------------
  demenergy = 0.D0
  do k = klow, kmax
     do j = jmin, jmax
        do i = ilow, imax
           if (Symmetry == AXISYM) then
              temp = exp(10.d0*phi(i,j,k)) * &
                   ( gxx(i,j,k)*(Ex(i,j,k)**2 + Bx(i,j,k)**2) + &
                   2.d0*gxy(i,j,k)*(Ex(i,j,k)*Ey(i,j,k) + Bx(i,j,k)*By(i,j,k)) + &
                   2.d0*gxz(i,j,k)*(Ex(i,j,k)*Ez(i,j,k) + Bx(i,j,k)*Bz(i,j,k)) + &
                   gyy(i,j,k)*(Ey(i,j,k)**2 + By(i,j,k)**2) + &
                   2.d0*gyz(i,j,k)*(Ey(i,j,k)*Ez(i,j,k) + By(i,j,k)*Bz(i,j,k)) + &
                   gzz(i,j,k)*(Ez(i,j,k)**2 + Bz(i,j,k)**2) ) * abs(X(i,1,1))

              temp = temp*mskf(i,j,k)

              demenergy = demenergy + temp

           else
              temp = exp(10.d0*phi(i,j,k)) * &
                   ( gxx(i,j,k)*(Ex(i,j,k)**2 + Bx(i,j,k)**2) + &
                   2.d0*gxy(i,j,k)*(Ex(i,j,k)*Ey(i,j,k) + Bx(i,j,k)*By(i,j,k)) + &
                   2.d0*gxz(i,j,k)*(Ex(i,j,k)*Ez(i,j,k) + Bx(i,j,k)*Bz(i,j,k)) + &
                   gyy(i,j,k)*(Ey(i,j,k)**2 + By(i,j,k)**2) + &
                   2.d0*gyz(i,j,k)*(Ey(i,j,k)*Ez(i,j,k) + By(i,j,k)*Bz(i,j,k)) + &
                   gzz(i,j,k)*(Ez(i,j,k)**2 + Bz(i,j,k)**2) )

              demenergy = demenergy + temp
           end if
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  demenergy = demenergy * dV/(8.d0*PI)
  return
end subroutine em_energy_disk

