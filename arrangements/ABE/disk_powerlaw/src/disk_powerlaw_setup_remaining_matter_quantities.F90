#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine disk_powerlaw_setup_remaining_matter_quantities(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext

  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  integer :: vindex,handle,ierr,foundflag
  real*8  :: rho_max,tau_max
  real*8  :: alpha,e4phi,shift_x,shift_y,shift_z,beta2,u_t,psi6,u_x,u_y,u_z
  real*8  :: ONE,ZERO,TWO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0)

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
           e4phi = exp(4.d0*phi(i,j,k))
           shift_x = e4phi*(gxx(i,j,k)*shiftx(i,j,k) + & 
                gxy(i,j,k)*shifty(i,j,k) + gxz(i,j,k)*shiftz(i,j,k))
           shift_y = e4phi*(gxy(i,j,k)*shiftx(i,j,k) + & 
                gyy(i,j,k)*shifty(i,j,k) + gyz(i,j,k)*shiftz(i,j,k)) 
           shift_z = e4phi*(gxz(i,j,k)*shiftx(i,j,k) + & 
                gyz(i,j,k)*shifty(i,j,k) + gzz(i,j,k)*shiftz(i,j,k)) 
!!$
           !compute physical u_i
           u_x = u0(i,j,k)*(shift_x + &
                e4phi*gxx(i,j,k)*vx(i,j,k) + &
                e4phi*gxy(i,j,k)*vy(i,j,k) + &
                e4phi*gxz(i,j,k)*vz(i,j,k))
           
           u_y = u0(i,j,k)*(shift_y + &
                e4phi*gxy(i,j,k)*vx(i,j,k) + &
                e4phi*gyy(i,j,k)*vy(i,j,k) + &
                e4phi*gyz(i,j,k)*vz(i,j,k))
           
           u_z = u0(i,j,k)*(shift_z + &
                e4phi*gxz(i,j,k)*vx(i,j,k) + &
                e4phi*gyz(i,j,k)*vy(i,j,k) + &
                e4phi*gzz(i,j,k)*vz(i,j,k))

           
           alpha = lapm1(i,j,k)+1.d0
           P(i,j,k)=K_poly*rho_b(i,j,k)**gamma_th
           h(i,j,k)=1.d0+gamma_th/(gamma_th-1.0)*P(i,j,k)/rho_b(i,j,k)
           !w(i,j,k)=rho_star(i,j,k)*alpha*u0(i,j,k)
           
!!$           beta2 = shiftx(i,j,k)*shift_x + shifty(i,j,k)*shift_y + shiftz(i,j,k)*shift_z
!!$           u_t = u0(i,j,k)*(-alpha*alpha + beta2 + &
!!$                shift_x*vx(i,j,k) + shift_y*vy(i,j,k)+shift_z*vz(i,j,k))
         
           if ((i.eq.42).and.(j.eq.42).and.(k.eq.5).and.(1.eq.0)) then
              write(*,*) "*********************"
              write(*,*) "X: ",X(i,j,k)
              write(*,*) "Y: ",Y(i,j,k)
              write(*,*) "Z: ",Z(i,j,k)
              write(*,*) "vx: ",vx(i,j,k)
              write(*,*) "vy: ",vy(i,j,k)
              write(*,*) "vz: ",vz(i,j,k)
              write(*,*) "u_x: ",u_x
              write(*,*) "u_y: ",u_y
              write(*,*) "u_z: ",u_z
              write(*,*) "shiftx: ",shiftx(i,j,k)
              write(*,*) "shifty: ",shifty(i,j,k)
              write(*,*) "shiftz: ",shiftz(i,j,k)
              write(*,*) "gxx: ",gxx(i,j,k)
              write(*,*) "gxy: ",gxy(i,j,k)
              write(*,*) "gxz: ",gxz(i,j,k)
              write(*,*) "gyy: ",gyy(i,j,k)
              write(*,*) "gyz: ",gyz(i,j,k)
              write(*,*) "gzz: ",gzz(i,j,k)
              write(*,*) "**********************"
           endif
           psi6 = exp(6.d0*phi(i,j,k))
           w(i,j,k)=rho_star(i,j,k)*alpha*u0(i,j,k)
        !  tau(i,j,k) = -rho_star(i,j,k)*(1.d0+h(i,j,k)*u_t)-alpha*psi6*P(i,j,k)
           tau(i,j,k) = w(i,j,k)*h(i,j,k) - psi6*P(i,j,k) - rho_star(i,j,k)
           st_x(i,j,k)     = rho_star(i,j,k)*h(i,j,k)*u_x
           st_y(i,j,k)     = rho_star(i,j,k)*h(i,j,k)*u_y
           st_z(i,j,k)     = rho_star(i,j,k)*h(i,j,k)*u_z
	end do
     end do
  end do
  
 
  !mhd_st_x=st_x
  !mhd_st_y=st_y
  !mhd_st_z=st_z

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
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)

  call CartSymGN(ierr,cctkGH,'BSSN::BSSN_vars')
!  call CCTK_SyncGroup(ierr,cctkGH,'BSSN::BSSN_vars') 

  !Initialize these variables:
  t_last_disk_powerlaw = 0.d0
!!!!!!FIX THIS!!!!!!!!!!!!!
  !  xh1_last_bondi = xh0
  !  yh1_last_bondi = yh0
  !  zh1_last_bondi = zh0
  !  xh2_last_bondi = xh0
  !  yh2_last_bondi = yh0
  !  zh2_last_bondi = zh0
  ! ah_radii_last_bondi = r0_ah

  ah1_radii_last_disk_powerlaw = 1.d0
  ah2_radii_last_disk_powerlaw = 1.d0
!temporary
  xh1_last_disk_powerlaw = -10.0125d0
  xh2_last_disk_powerlaw = 10.0125d0
  yh1_last_disk_powerlaw = 0.0
  yh2_last_disk_powerlaw = 0.0
  zh1_last_disk_powerlaw = 0.d0
  zh2_last_disk_powerlaw = 0.d0
!!!!!!FIX THIS!!!!!!!!!!!!!

  surf1_radii_last_disk_powerlaw = 0.d0
  surf2_radii_last_disk_powerlaw = 0.d0
  surf3_radii_last_disk_powerlaw = 0.d0
  surf4_radii_last_disk_powerlaw = 0.d0

end subroutine disk_powerlaw_setup_remaining_matter_quantities
