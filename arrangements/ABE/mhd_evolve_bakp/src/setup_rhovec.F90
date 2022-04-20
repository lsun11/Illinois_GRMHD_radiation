!-----------------------------------------------------------------------------
!
!$Id: correct_mhd.F90  $
!
!-----------------------------------------------------------------------------
!
! Correct & Correct RHS routines for mhd variables
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
!
! Corrector
!
!-----------------------------------------------------------------------------
subroutine Setup_Rhovec(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  interface
     REAL*8 FUNCTION find_rest_density(x1,x2,r,Gam,mhd_mdot,rho_inf,mhd_mbh,KPOLY)
       implicit none
       REAL*8 x1,x2,r,Gam,mhd_mdot,rho_inf,mhd_mbh,KPOLY
     end FUNCTION find_rest_density
  end interface

  integer                                  :: i,j,k
  integer, dimension(3)                    :: ext,global_ext
  real*8 :: pi,gam,u2,m0,u,a2,rho0_i
  real*8 :: rmax
  real*8 :: a2_inf,rho_inf,C_inf,dk,rr
  real*8 :: lower,upper,cl,cu,riso,kl,ku,fac,xi,rhorm,rhoeps,eps
  real*8, dimension(4) :: xa,ya
  integer :: fnr,err,jj,handle,index
  integer :: NO_SYMM,EQUATORIAL,OCTANT,PI_SYMM,AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: ZERO,ONE,TWO,FOUR
  parameter(ZERO=0.0d0,ONE=1.0d0, TWO=2.0d0, FOUR=4.0d0)

  integer :: nr,ict
  parameter (nr=10000)
  real*8, dimension(-nr:nr)  :: radius,rho0,smallM

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"fisheye::PhysicalRadius")
  call CCTK_Reduce (err,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,rmax,1,index)

  dxvec=2.0*rmax/(1.0*nrhovec)
  write(6,*)'nrhovec:',nrhovec,' dxvec:',dxvec,' rmax:',rmax

  pi=3.141592653
  Gam = gamma_th 
  u2 = mhd_mbh/(TWO*mhd_r_crit)
  u = sqrt(u2)
  !Eq. G17
  a2 = u2/(ONE - 3.D0*u2)
  !Eq. G28
  rho0_i = ( a2/K_Poly/(Gam - (Gam/(Gam-ONE))*a2) )**(ONE/(Gam-ONE))

  ! no need to rescale k_poly; it's already done

  !Eq. G30     
  a2_inf = Gam-ONE - sqrt(ONE + 3.D0*a2)*(Gam-ONE - a2)
  !Eq. G28
  rho_inf = ( a2_inf/K_POLY/(Gam - (Gam/(Gam-ONE))*a2_inf) )**(ONE/(Gam-ONE))
  !Eq. G29, this is the integration constant
  C_INF = ( ONE + K_POLY*(Gam/(Gam-ONE))*(rho_inf**(Gam-ONE)) )**2
  fnr = nr
  dk = (0.1d0/mhd_r_crit)**(ONE/nr)
  !This is the density at the sonic point
  rho0(1) = rho0_i
  radius(1) = mhd_r_crit
  rr=mhd_r_crit
  ! fill up 1D arrays
  do k = 2, nr
     !we work our way inward from the sonic point in geomtric increments
     rr = rr * dk
     lower = 0.9D0 * rho0(k-1)
     upper = 1.1D0 * rho0(k-1)
     radius(k) = rr
     ! these are the integration constants for densities increasing
     ! and decreasing away geometrically
     Cl = ( ( ONE + K_POLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*mhd_mbh/rr + (mhd_mdot/(4.D0*PI*lower*rr*rr))**2) - C_INF
     Cu = ( ( ONE + K_POLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*mhd_mbh/rr + (mhd_mdot/(4.D0*PI*upper*rr*rr))**2) - C_INF
     !now bracket
     ict=0
     do while (Cl*Cu.gt.ZERO .and. upper-lower.gt.1.d-3*rho0(k+1))
        lower = sqrt(lower*rho0(k-1))
        upper = sqrt(upper*rho0(k-1))
        Cl = ( ( ONE + K_POLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mhd_mbh/rr + (mhd_mdot/(4.D0*PI*lower*rr*rr))**2) - C_INF
        Cu = ( ( ONE + K_POLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mhd_mbh/rr + (mhd_mdot/(4.D0*PI*upper*rr*rr))**2) - C_INF
        ict=ict+1
        if(ict.gt.nr) then
           write(6,*)'Endless loop in setup_rhovec?',k,Cl,Cu,upper,lower,rho0(k+1)
           stop
        endif
     end do
     rho0(k) = find_rest_density(lower,upper,rr,Gam,mhd_mdot,rho_inf,mhd_mbh,K_POLY)
     rho0_i = rho0(k)
  end do
  rr = mhd_r_crit
  do k = 0, -nr+1, -1
     rr = rr / dk
     lower = 0.8D0 * rho0(k+1)
     upper = 1.2D0 * rho0(k+1)
     radius(k) = rr
     
     Cl = ( ( ONE + K_POLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*mhd_mbh/rr + (mhd_mdot/(4.D0*PI*lower*rr*rr))**2) - C_INF
     Cu = ( ( ONE + K_POLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*mhd_mbh/rr + (mhd_mdot/(4.D0*PI*upper*rr*rr))**2) - C_INF
     
     ict=0
     do while (Cl*Cu.gt.ZERO .and. upper-lower.gt.1.d-4*rho0(k+1))
        lower = sqrt(lower*rho0(k+1))
        upper = sqrt(upper*rho0(k+1))
        Cl = ( ( ONE + K_POLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mhd_mbh/rr + (mhd_mdot/(4.D0*PI*lower*rr*rr))**2) - C_INF
        Cu = ( ( ONE + K_POLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mhd_mbh/rr + (mhd_mdot/(4.D0*PI*upper*rr*rr))**2) - C_INF
        ict=ict+1
        if(ict.gt.nr) then
           write(6,*)'Endless loop in setup_rhovec?',k,Cl,Cu,upper,lower,rho0(k+1)
           stop
        endif
     end do
     
     rho0(k) = find_rest_density(lower,upper,rr,Gam,mhd_mdot,rho_inf,mhd_mbh,K_POLY)
     rho0_i = rho0(k)
     u2 = (mhd_mdot/(4.D0*PI*rho0_i*rr*rr))**2
  end do
  
  do i=1,nrhovec
     riso=dxvec*i
     rr=riso*(1.0+0.5/riso)**2
     kl = floor( fnr*log(mhd_r_crit/rr)/log(mhd_r_crit/0.1d0) ) + 1
     ku = kl + 1
     fac = (rr - radius(kl))/(radius(ku) - radius(kl))
     if(fac.lt.ZERO .or. fac.gt.ONE) then
        write(*,*) "Error: fac = ",fac
        write(*,*) "    ", kl, ku, rr, radius(kl), radius(ku)
     end if
     do jj=1,4
        xa(jj) = log( radius(jj-2+kl) )
        ya(jj) =  rho0(jj-2+kl)
     end do
     xi = log(rr)
     !returns the density at the given radial point to rho0_i
     call polint(xa,ya,4,xi,rhorm,err)
     
     rhovec(i)=rhorm
     rhoeps = K_POLY * (rhorm**Gam) / (Gam - ONE)
     eps = rhoeps / rhorm
     Pvec(i) = (Gam - ONE)*rhoeps
     u = mhd_mdot / (FOUR * PI * rr * rr * rhorm)
     !uiso from my notes
     vvec(i)= u/(1.0d0+mhd_mbh/2.0d0/riso)/(1.0d0-mhd_mbh/2.0d0/riso)
     
  enddo
  
end subroutine Setup_Rhovec
