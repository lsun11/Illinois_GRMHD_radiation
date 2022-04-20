#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_initialdata_local2(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext

  real*8 			           :: dT,dX,dY,dZ,xL,yL,zL,rL,r_falloff_radius,r_falloff_dr,rhobatm_initial_value,rhobatm_final_value,rhobatm
  real*8 			           :: n1,n2,n3,mf

  integer :: vindex,handle,ierr
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

  ! the v's are really the u_i's
  u0=sqrt(1.0+psi**(-4)*(vx*vx+vy*vy+vz*vz))/(lapm1+1.0)

  rho_star=rho_b*(lapm1+1.0)*u0*psi**6

if(1==0) then
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           xL=X(i,j,k)
           yL=Y(i,j,k)
           zL=Z(i,j,k)
           rL=sqrt(xL*xL+yL*yL+zL*zL)

           r_falloff_radius = 15
           r_falloff_dr     = 3

           !UNPERT value: 1.5011293377E-001

           !rhobatm_initial_value = 1.2534469659822e-8
           !rhobatm_final_value = 1.2534469659822e-11
           !1.5538657586E-001 
           !rhobatm_initial_value = 1.D-7
           !rhobatm_final_value = 1.D-10

           !1.5053881771E-001
           !rhobatm_initial_value = 1.D-8
           !rhobatm_final_value = 1.D-11

           !1.5040410350E-001
           !rhobatm_initial_value = 7.D-9
           !rhobatm_final_value = 1.D-11

           !1.5013467510E-001
           rhobatm_initial_value = 1.D-9
           rhobatm_final_value = 1.D-11

           rhobatm=-erf((rL-r_falloff_radius-r_falloff_dr*2.0)/r_falloff_dr)*(rhobatm_initial_value-rhobatm_final_value)*0.5D0 + (rhobatm_initial_value-rhobatm_final_value)*0.5D0 + rhobatm_final_value

           rho_b(i,j,k) = max(rho_b(i,j,k),rhobatm)

           rho_star(i,j,k)=rho_b(i,j,k)*(lapm1(i,j,k)+1.0)*u0(i,j,k)*psi(i,j,k)**6

           if(rho_star(i,j,k).lt.0.D0) then
              write(*,*) "ERROR IN RHO_STAR:",i,j,k,rho_star(i,j,k)
              stop
           end if
        end do
     end do
  end do
end if

  P=K_poly*rho_b**(gamma_th)
  h=1.0+gamma_th/(gamma_th-1.0)*P/rho_b
  w=rho_star*(lapm1+1.0)*u0
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
     write(*,*) "I PITY DA FOO TRIES TO USE AXISYMMETRY FOR ORBITING BLACK HOLES!"
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

end subroutine bhns_initialdata_local2
