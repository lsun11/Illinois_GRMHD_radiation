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

  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  integer :: vindex,handle,ierr
  real*8  :: rho_max,tau_max,valuetosetvz
  real*8  :: ONE,ZERO,W_L,U2, gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,psiL,psiL4
  real*8  :: cn2, cn3, an1, an2, an3
  real*8  :: gupxxL, gupxyL, gupxzL, gupyyL, gupyzL, gupzzL, u0L
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

!!$  ! ZACH SAYS: TESTING
!!$  !vz = 0.006 ! 0.1*vy
!!$  !vz = 0.012 ! 0.1*vy
!!$!  where(vy.ne.0.D0) 
!!$!     vz = 0.003*psi**4 ! 0.1*vy
!!$!  end where
!!$  
!!$  valuetosetvz = 0.006 ! <- this is about 5% kick
!!$!  valuetosetvz = 0.003 ! <- this is about 0.025 kick
!!$  !valuetosetvz = 0.012
!!$  where(vy.ne.0.D0)
!!$     vz = sqrt( (1+psi**(-4)*(vx**2+vy**2))*(valuetosetvz)**2 ) / sqrt(psi**(-8)*(lapm1+1.0)**2 - psi**(-4)*(valuetosetvz)**2)
!!$  end where


write(*,*) 'In bhns local2, rho_b_atm = ',rho_b_atm,dX


  ! the v's are really the u_i's when use_new_bhns_initial_data<2
  if (use_new_bhns_initial_data.ne.2) then
     u0=sqrt(1.0+psi**(-4)*(vx*vx+vy*vy+vz*vz))/(lapm1+1.0d0)
  else	
     ! The v's are the U's determined through u^a = W(n^a + U^a) when use_new_bhns_initial_data==2
     ! The Lorentz factor as measured by a normal observer is W = 1/Sqrt(1-g_ij U^i U^j)
     ! Knowing W we can find u^0 = W/lapse
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              psiL=psi(i,j,k)
              psiL4=psiL*psiL*psiL*psiL
              ! Full 3-metric
              gxxL=psiL4*gxx(i,j,k)
              gxyL=psiL4*gxy(i,j,k)
              gxzL=psiL4*gxz(i,j,k)
              gyyL=psiL4*gyy(i,j,k)
              gyzL=psiL4*gyz(i,j,k)
              gzzL=psiL4*gzz(i,j,k)
              if(rho_b(i,j,k) .le. 0.d0) then
                 vx(i,j,k) = 0.d0
                 vy(i,j,k) = 0.d0
                 vz(i,j,k) = 0.d0
              end if
              U2 = gxxL*vx(i,j,k)*vx(i,j,k) + gyyL*vy(i,j,k)*vy(i,j,k) + gzzL*vz(i,j,k)*vz(i,j,k) + &
                   2.d0*(gxyL*vx(i,j,k)*vy(i,j,k) + gxzL*vx(i,j,k)*vz(i,j,k) + gyzL*vy(i,j,k)*vz(i,j,k) )
              if (abs(U2-1.d0).le.1.0d-3) then
                 write(*,*) "Stopping: can't have near luminal initial data: check your data"
                 stop
              end if
              W_L = 1.d0/Sqrt(1.d0-U2)
              ! Compute u^0
              u0(i,j,k) =W_L/(lapm1(i,j,k)+1.0d0)
              ! And now redefine U^i to correspond to u_i = W * U_i
              vx(i,j,k) = W_L*( gxxL*vx(i,j,k)+gxyL*vy(i,j,k)+gxzL*vz(i,j,k) )
              vy(i,j,k) = W_L*( gxyL*vx(i,j,k)+gyyL*vy(i,j,k)+gyzL*vz(i,j,k) )
              vz(i,j,k) = W_L*( gxzL*vx(i,j,k)+gyzL*vy(i,j,k)+gzzL*vz(i,j,k) )
           enddo
        enddo
     enddo
  end if


  rho_b= max(rho_b,rho_b_atm)

  rho_star=rho_b*(lapm1+1.0d0)*u0*psi**6


  P=K_poly*rho_b**(gamma_th)

  h=1.0+gamma_th/(gamma_th-1.0)*P/rho_b
  w=rho_star*(lapm1+1.0d0)*u0
  tau=w*h-psi**(6)*P-rho_star


  st_x = rho_star*h*vx
  st_y = rho_star*h*vy
  st_z = rho_star*h*vz
  mhd_st_x=st_x
  mhd_st_y=st_y
  mhd_st_z=st_z

  ! and really the v's: v^i = u^i/u^0
  if (use_new_bhns_initial_data.ne.2) then
     vx=psi**(-4)*vx/u0-shiftx
     vy=psi**(-4)*vy/u0-shifty
     vz=psi**(-4)*vz/u0-shiftz
  else if (use_new_bhns_initial_data==2) then
     ! Remember now v^i are actually u_i
     ! and we want v^i = u^i/u^0
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              psiL=psi(i,j,k)
              psiL4=psiL*psiL*psiL*psiL
              ! Full inverse metric
              gupxxL=gupxx(i,j,k)/psiL4
              gupxyL=gupxy(i,j,k)/psiL4
              gupxzL=gupxz(i,j,k)/psiL4
              gupyyL=gupyy(i,j,k)/psiL4
              gupyzL=gupyz(i,j,k)/psiL4
              gupzzL=gupzz(i,j,k)/psiL4
              u0L = u0(i,j,k)
              ! And now set v^i = g^{ij}u_i/u^0 - beta^i
              vx(i,j,k) = ( gupxxL*vx(i,j,k)+gupxyL*vy(i,j,k)+gupxzL*vz(i,j,k) )/u0L - shiftx(i,j,k)
              vy(i,j,k) = ( gupxyL*vx(i,j,k)+gupyyL*vy(i,j,k)+gupyzL*vz(i,j,k) )/u0L - shifty(i,j,k)
              vz(i,j,k) = ( gupxzL*vx(i,j,k)+gupyzL*vy(i,j,k)+gupzzL*vz(i,j,k) )/u0L - shiftz(i,j,k)
           enddo
        enddo
     enddo
  end if

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

  !call CartSymGN(ierr,cctkGH,'BSSN::BSSN_vars')
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

  if(subtract_off_Omega_r==1) then
     vx = vx + Omega_value*y
     vy = vy - Omega_value*x
  end if

end subroutine bhns_initialdata_local2
