#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
!
! Function to compute 1/8\pi \int \psi^6 (x A^m_y - yA^m_x) dS_m
!
!-----------------------------------------------------------------------------
subroutine setup_multiple_rhosurfaces(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  real*8                                   :: dX,dY,dZ
  real*8, dimension(1,3)                   :: pointcoords
  real*8                                   :: PI,costheta,sintheta,phiangle  
  integer                                  :: i,n,j,int_order
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  real*8 :: drhosurf,rr,getu0

  if(nsurfrho.gt.20) then
     write(6,*)'nsurfrho cannot be bigger than 20, since we cannot dynamically'
     write(6,*)'allocate the vectors in rhosurf_params. If you need more space'
     write(6,*)'edit bondi/interface.ccl'
     stop
  else if(nsurfrho.eq.0) then
     return
  else if (nsurfrho.ge.1) then
     rhosurfvec(1)=rhosurf
     u0vec(1)=getu0(rhosurf,gamma_th,mdot,r_crit,mbh)
     rplus(1)=0.
     rminus(1)=0.
     if(nsurfrho.eq.1)return

     write(6,*)'rhosurf 1:',rhosurfvec(1),u0vec(1)

  endif
  
  if(arithrhosurf==1) then
     drhosurf=(rhosurf2-rhosurf)/(nsurfrho-1.0d0)
  else
     drhosurf=exp((log(rhosurf2)-log(rhosurf))/(nsurfrho-1.0d0))
  endif
  
  do i=2,nsurfrho

     if(arithrhosurf==1) then
        rr=rhosurf+(i-1)*drhosurf
     else
        rr=rhosurf*drhosurf**(i-1)
     endif

     rhosurfvec(i)=rr
     u0vec(i)=getu0(rr,gamma_th,mdot,r_crit,mbh)
     rplus(i)=0.
     rminus(i)=0.
     write(6,*)'surf ',i,':rhosurf:',rhosurfvec(i),u0vec(i)
     
  enddo
  
  
end subroutine setup_multiple_rhosurfaces

real*8 function getu0(rhos,gam,md,rc,mh)
  implicit none
  real*8 :: rhos,gam,md,rc,mh
  
  real*8 ::  kpoly
  
  real*8 ::  ZERO,HALF,ONE,TWO,THREE,FOUR,PI
  parameter(ZERO=0.d0,HALF=0.5d0,ONE=1.d0,TWO=2.d0,THREE=3.d0,FOUR=4.d0)
  
  real*8 ::  rnew,u,r,cinf,rho0inf,a2inf,mdc,a2,rnewer
  real*8 ::  rho0c,a2c,uc,u2c,ur,unew,alph
  integer i,k,ii
  real*8 :: rhopnf,fnr,dk,lower,upper,cl,cu,rest_density,frac
  integer :: nr
  parameter(nr=10000)
  real*8, dimension(-nr:nr)  :: radius,rho0,smallM
  
  PI=acos(-ONE)
  
  kpoly=1.0d0
  
  !     Eq g17
  u2c=mh/(TWO*rc)
  uc=sqrt(u2c)
  a2c=u2c/(ONE-THREE*u2c)
  !     Eq g28
  rho0c=(a2c/Kpoly/(gam-(gam/(gam-ONE))*a2c))**(ONE/(Gam-ONE))
  !     Eq g33
  Mdc=Four*PI*uc*rho0c*rc**2
  
  write(6,*)'md,mdc:',md,mdc,u2c,uc,a2c,rho0c
  
  !     rescale mdot
  kpoly=(mdc/md)**(gam-ONE)
  rho0c=(a2c/Kpoly/(gam-(gam/(gam-ONE))*a2c))**(ONE/(Gam-ONE))
  Mdc=Four*PI*uc*rho0c*rc**2
  
   write(6,*)'new kpoly:',kpoly,rho0c,mdc

 !     eq g30
  a2inf=gam-ONE-sqrt(ONE+THREE*a2c)*(gam-ONE-a2c)
  !     eq g28
  rho0inf=(a2inf/Kpoly/(gam-(gam/(gam-ONE))*a2inf))**(ONE/(Gam-ONE))
  !     eq g29
  cinf=(ONE+kpoly*(gam/(gam-ONE))*(rho0inf**(gam-ONE)))**2
  
  write(6,*)'cinf:',cinf,a2inf,rho0inf,nr
  
  fnr = nr
  dk = (0.1d0/rc)**(ONE/nr)
  !This is the density at the sonic point
  rho0(1) = rho0c
  radius(1) = rc
  r=rc
  write(6,*)'fnr:',fnr,dk,rho0c,rc
  ! fill up 1D arrays
  do k = 2, nr
     !we work our way inward from the sonic point in geomtric increments
     r = r * dk
     lower = 0.9D0 * rho0(k-1)
     upper = 1.1D0 * rho0(k-1)
     radius(k) = r
     
     ! these are the integration constants for densities increasing
     ! and decreasing away geometrically
     Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*mh/r + (Md/(4.D0*PI*lower*r*r))**2) - CINF
     Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*mh/r + (Md/(4.D0*PI*upper*r*r))**2) - CINF
!     write(6,*)'k:',k,r,lower,upper,radius(k),cl,cu
     !now bracket
     do while (Cl*Cu.gt.ZERO .and. upper-lower.gt.1.d-3*rho0(k+1))
        lower = sqrt(lower*rho0(k-1))
        upper = sqrt(upper*rho0(k-1))
        Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mh/r + (Md/(4.D0*PI*lower*r*r))**2) - CINF
        Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mh/r + (Md/(4.D0*PI*upper*r*r))**2) - CINF
     end do
     
     rho0(k) = rest_density(lower,upper,r,Gam,Md,rho0inf,mh,KPOLY)
     rho0c = rho0(k)
  end do
  write(*,*) "critical soln: ",rho0(1),a2,CINF
  write(*,*) "KPOLY = ",KPOLY
  r=rc
  do k = 0, -nr+1, -1
     r = r / dk
     lower = 0.8D0 * rho0(k+1)
     upper = 1.2D0 * rho0(k+1)
     radius(k) = r
     
     Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*mh/r + (Md/(4.D0*PI*lower*r*r))**2) - CINF
     Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*mh/r + (Md/(4.D0*PI*upper*r*r))**2) - CINF
     
     do while (Cl*Cu.gt.ZERO .and. upper-lower.gt.1.d-4*rho0(k+1))
        lower = sqrt(lower*rho0(k+1))
        upper = sqrt(upper*rho0(k+1))
        Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mh/r + (Md/(4.D0*PI*lower*r*r))**2) - CINF
        Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mh/r + (Md/(4.D0*PI*upper*r*r))**2) - CINF
     end do
     
     rho0(k) = rest_density(lower,upper,r,Gam,Md,rho0inf,mh,KPOLY)
     rho0c = rho0(k)
     u2c = (Md/(4.D0*PI*rho0c*r*r))**2
  end do
  write(*,*) "rho limits: ",rho0(-nr+1),rho0(nr)
  write(*,*) "r limits: ",radius(-nr+1),radius(nr)

  do i=-nr+2,nr-2
     if(rho0(i).lt.rhos.and.rho0(i+1).gt.rhos) then
        ii=i
        goto 16
     endif
  enddo
16 continue
  frac=(rhos-rho0(i))/(rho0(i+1)-rho0(i))
  r=radius(i)+frac*(radius(i+1)-radius(i))
  u=md/4.0d0/pi/r**2/rhos
  write(6,*)'found the point:',i,r,u
     
  write(6,*)'check:',r,u,4.0*pi*r**2*rhos*u,md, &
       (rhopnf(rhos,gam,kpoly))**2*(one-two*mh/r+u**2),cinf

  alph=sqrt(1.0-2*mh/r)

  getu0=sqrt(1.0+u**2/alph**2)/alph

  return
end function getu0

real*8 function rhopnf(n,gam,k)
  implicit none
  real*8 n,gam,k
  
  rhopnf=1.0+gam/(gam-1.d0)*k*n**(gam-1.0)
  return
end function rhopnf

REAL*8 FUNCTION rest_density(x1,x2,r,Gamma,Mdot,rho0inf,MBH,KPOLY)
  implicit none
  INTEGER MAXIT
  REAL*8 x1,x2,r,Gamma,Mdot,rho0inf,MBH,KPOLY
  REAL*8 xacc, scale
  logical success
  !    EXTERNAL funcd
  PARAMETER (MAXIT=100)
  INTEGER j
  REAL*8 df,dx,dxold,f,fh,fl,temp,xh,xl
  xacc = 1.D-16
  success=.true.
  call transf2(x1,fl,df,r,Gamma,Mdot,rho0inf,MBH,KPOLY)
  call transf2(x2,fh,df,r,Gamma,Mdot,rho0inf,MBH,KPOLY)
  if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) then
     rest_density = x2 - 2.D0*fh/df
     call transf2(rest_density,f,df,r,Gamma,Mdot,rho0inf,MBH,KPOLY)
     if(abs(f).lt.xacc)then
        return
     else
        !          write(*,*) "bracket: ",r,fl,fh
        !       x2 = 1.1 * x1
     end if
  end if
  if(abs(fl).lt.xacc)then
     rest_density=x1
     return
  else if(abs(fh).lt.xacc)then
     rest_density=x2
     return
  else if(fl.lt.0.)then
     xl=x1
     xh=x2
  else
     xh=x1
     xl=x2
  endif
  rest_density=.5*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  call transf2(rest_density,f,df,r,Gamma,Mdot,rho0inf,MBH,KPOLY)
  do j=1,MAXIT
     if(((rest_density-xh)*df-f)*((rest_density-xl)*df-f).ge.0. &
          .or. abs(2.*f).gt.abs(dxold*df) ) then
        dxold=dx
        dx=0.5D0*(xh-xl)
        rest_density=xl+dx
        if(xl.eq.rest_density) return
     else
        dxold=dx
        dx=f/df
        temp=rest_density
        rest_density=rest_density-dx
        if(temp.eq.rest_density) return
     endif
     if(abs(dx).lt.xacc) return
     call transf2(rest_density,f,df,r,Gamma,Mdot,rho0inf,MBH,KPOLY)
     if(f.lt.0.) then
        xl=rest_density
     else
        xh=rest_density
     endif
  end do
  success=.false.
  return
END FUNCTION rest_density

subroutine transf2(rho0,f,df,r,Gamma,Mdot,rho_inf,MBH,KPOLY)
  implicit none
  real*8            :: rho0,f,df,r,Gamma,Mdot,MBH,KPOLY
  real*8            :: PI,ONE,TWO,u2,C_INF,rho_inf
  ONE = 1.D0
  TWO = 2.D0
  PI = acos(-ONE)
  C_INF = ( ONE + KPOLY*(Gamma/(Gamma-ONE))*(rho_inf**(Gamma-ONE)) )**2
  u2 = (Mdot/(4.D0*PI*rho0*r*r))**2
  f = ( ( ONE + KPOLY*(Gamma/(Gamma-ONE))*(rho0**(Gamma-ONE)) )**2 ) &
       * (ONE - TWO*MBH/r + u2) - C_INF
  df = TWO * (ONE + KPOLY*(Gamma/(Gamma-ONE)) * rho0**(Gamma-ONE)) &
       * KPOLY * Gamma * rho0**(Gamma-TWO) * (ONE - TWO*MBH/r + u2) &
       + u2 * (-TWO/rho0) * ((ONE + KPOLY*(Gamma/(Gamma-ONE)) &
       * rho0**(Gamma-ONE))**2)
  !    print '(1X,4e15.4)', rho0,u2,f,TWO/r
end subroutine transf2
