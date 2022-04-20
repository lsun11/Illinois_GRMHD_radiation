#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bbh_bondi_matter_local_2(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext

  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  integer :: vindex,handle,ierr,foundflag
  real*8  :: rho_max,tau_max
  real*8  :: ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)

  INTEGER :: i,j,k
  character :: filename*30,c2*2,c1
  CCTK_REAL :: reduction_value

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

   

  rho_b=max(rho_b,rho_b_atm)
  !u0=1.d0/sqrt((lapm1+1.d0)**2-psi**4*(vx**2+vy**2+vz**2))
  
  rho_star=rho_b*(lapm1+1.0)*u0*psi**6

  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
     if(rho_star(i,j,k).lt.0.D0) then
        write(*,*) "ERROR IN RHO_STAR:",i,j,k,rho_star(i,j,k)
        stop
     end if
  end do
  end do
  end do

  P=K_poly*rho_b**(gamma_th)
  h=1.0+gamma_th/(gamma_th-1.0)*P/rho_b
  w=rho_star*(lapm1+1.0)*u0
  tau=w*h-psi**(6)*P-rho_star
  
  st_x = rho_star*h*psi**4*vx*u0
  st_y = rho_star*h*psi**4*vy*u0
  st_z = rho_star*h*psi**4*vz*u0
  mhd_st_x=st_x
  mhd_st_y=st_y
  mhd_st_z=st_z

  if(fisheye_enable .eq. 1) then
     call trans_phys_fish_tensor_flat(cctkGH,cctk_lsh,Symmetry, &
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz)
     call trans_phys_fish_tensor_inv(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
     call trans_phys_fish_tensor(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
     call trans_phys_fish_phi(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,phi)
     call trans_phys_fish_gamt_flat(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,RadiusDerivative2, &
          Gammax,Gammay,Gammaz)
     call trans_phys_fish_matter(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative, &
          rho_star,tau,st_x,st_y,st_z,w,vx,vy,vz,Symmetry)
  end if

  mhd_st_x=st_x
  mhd_st_y=st_y
  mhd_st_z=st_z

  if(Symmetry==AXISYM) then
     write(*,*) "I PITY DA FOO TRIES TO USE AXISYMMETRY FOR ORBITING BLACK HOLES!"
     stop
  end if

  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry) 
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry) 
!  call Derivs(ext,X,Y,Z,dX,dY,dZ,shiftx,shiftxx,shiftxy,shiftxz,Symmetry) 
!  call Derivs(ext,X,Y,Z,dX,dY,dZ,shifty,shiftyx,shiftyy,shiftyz,Symmetry) 
!  call Derivs(ext,X,Y,Z,dX,dY,dZ,shiftz,shiftzx,shiftzy,shiftzz,Symmetry) 
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)
!  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,shiftx,shiftxx,shiftxy,shiftxz)
!  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,shifty,shiftyx,shiftyy,shiftyz)
!  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,shiftz,shiftzx,shiftzy,shiftzz)

  call CartSymGN(ierr,cctkGH,'BSSN::BSSN_vars')
!  call CCTK_SyncGroup(ierr,cctkGH,'BSSN::BSSN_vars') 

  !Initialize these variables:
  t_last_bbh_bondi = 0.d0
!!!!!!FIX THIS!!!!!!!!!!!!!
  !  xh1_last_bondi = xh0
  !  yh1_last_bondi = yh0
  !  zh1_last_bondi = zh0
  !  xh2_last_bondi = xh0
  !  yh2_last_bondi = yh0
  !  zh2_last_bondi = zh0
  ! ah_radii_last_bondi = r0_ah

  ah1_radii_last_bbh_bondi = rah1_initial
  ah2_radii_last_bbh_bondi = rah2_initial
!temporary
  xh1_last_bbh_bondi = xbh1_initial
  xh2_last_bbh_bondi = xbh2_initial
  yh1_last_bbh_bondi = 0.0
  yh2_last_bbh_bondi = 0.0
  zh1_last_bbh_bondi = 0.d0
  zh2_last_bbh_bondi = 0.d0
!!!!!!FIX THIS!!!!!!!!!!!!!

  surf1_radii_last_bbh_bondi = 0.d0
  surf2_radii_last_bbh_bondi = 0.d0
  surf3_radii_last_bbh_bondi = 0.d0
  surf4_radii_last_bbh_bondi = 0.d0


  gxx_rhs=0.d0
  gxy_rhs=0.d0
  gxz_rhs=0.d0
  gyy_rhs=0.d0
  gyz_rhs=0.d0
  gzz_rhs=0.d0
  shiftx_rhs=0.d0
  shifty_rhs=0.d0
  shiftz_rhs=0.d0
  lapm1_rhs=0.d0
  phi_rhs=0.d0

end subroutine bbh_bondi_matter_local_2
