#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine wdns_initialdata_local2(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext

  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  integer :: vindex,handle,ierr
  real*8  :: rho_max,tau_max
  real*8  :: ONE,ZERO
  real*8  :: cn2, cn3, an1, an2, an3
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

  ! the v's are really the u_i's
  u0=sqrt(1.0+psi**(-4)*(vx*vx+vy*vy+vz*vz))/(lapm1+1.0d0)

  rho_star=rho_b*(lapm1+1.0d0)*u0*psi**6

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

  an1 = 1.0d0/(gamma_tab(1)-1.0d0)
  an2 = 1.0d0/(gamma_tab(2)-1.0d0)
  an3 = 1.0d0/(gamma_tab(3)-1.0d0)
  cn2 = (an1-an2)*k_tab(1)*rho_tab(1)**(1.0d0/an1)
  cn3 = cn2 + (an2-an3)*k_tab(2)*rho_tab(2)**(1.0d0/an2)

!  print *,"an1,an2,an3 =",an1,an2,an3
!  print *,"Gam1,Gam2,Gam3 =",gamma_tab(1),gamma_tab(2),gamma_tab(3)
!  print *,"kappa1,kappa2,kappa3 =",k_tab(1),k_tab(2),k_tab(3)
!  print *,"rhoo1,rhoo2 =",rho_tab(1),rho_tab(2)
!  stop

  where ( (rho_b-rho_tab(1))/rho_tab(1) <= 1.d-15)
  P=k_tab(1)*rho_b**(gamma_tab(1))
  h=1.0d0+(an1+1.0d0)*P/rho_b
  elsewhere ( ( (rho_b-rho_tab(1))/rho_tab(1) > 1.d-15 ).AND.( (rho_b-rho_tab(2))/rho_tab(2) <= 1.d-15 ) )
  P=k_tab(2)*rho_b**(gamma_tab(2))
  h=1.0d0+(an2+1.0d0)*P/rho_b + cn2
  elsewhere
  P=k_tab(3)*rho_b**(gamma_tab(3))
  h=1.0d0+(an3+1.0d0)*P/rho_b + cn3
  endwhere


  w=rho_star*(lapm1+1.0d0)*u0
  tau=w*h-psi**(6)*P-rho_star

  st_x = rho_star*h*vx
  st_y = rho_star*h*vy
  st_z = rho_star*h*vz
  mhd_st_x=st_x
  mhd_st_y=st_y
  mhd_st_z=st_z
  
  ! and really the v's
  vx=psi**(-4)*vx/u0-shiftx
  vy=psi**(-4)*vy/u0-shifty
  vz=psi**(-4)*vz/u0-shiftz

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

  !still need to calculate basic matter variables

  rho = h*w*exp(-6.0*phi)-P
  Sx = st_x*exp(-6.0*phi)
  Sy = st_y*exp(-6.0*phi)
  Sz = st_z*exp(-6.0*phi)
  Sxx = st_x*st_x/w/h*exp(-6.0*phi) + P*gxx*exp(4.0*phi)
  Sxy = st_x*st_y/w/h*exp(-6.0*phi) + P*gxy*exp(4.0*phi)
  Sxz = st_x*st_z/w/h*exp(-6.0*phi) + P*gxz*exp(4.0*phi)
  Syy = st_y*st_y/w/h*exp(-6.0*phi) + P*gyy*exp(4.0*phi)
  Syz = st_y*st_z/w/h*exp(-6.0*phi) + P*gyz*exp(4.0*phi)
  Szz = st_z*st_z/w/h*exp(-6.0*phi) + P*gzz*exp(4.0*phi)
  S = exp(-4.0*phi)*(gupxx*Sxx+gupyy*Syy+gupzz*Szz+ &
       2.0*(gupxy*Sxy+gupxz*Sxz+gupyz*Syz))

  mhd_st_x=st_x
  mhd_st_y=st_y
  mhd_st_z=st_z

  if(Symmetry==AXISYM) then
     write(*,*) "I PITY DA FOO TRIES TO USE AXISYMMETRY FOR ORBITING BINARY STARS!"
     stop
  end if

  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry) 
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry) 
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)

  call CartSymGN(ierr,cctkGH,'BSSN::BSSN_vars')
!  call CCTK_SyncGroup(ierr,cctkGH,'BSSN::BSSN_vars') 

! Initialize t_last, M0dot_last and int_M0dot: variables necessary to compute the
!  time integrated M0 flux across BH
  t_last = 0.d0
  M0dot_last = 0.d0
  int_M0dot = 0.d0

! Initialize the horizon position and radius

  xh_last = xh0
  yh_last = yh0
  zh_last = zh0
  ah_radii_last = r0_ah

end subroutine wdns_initialdata_local2
