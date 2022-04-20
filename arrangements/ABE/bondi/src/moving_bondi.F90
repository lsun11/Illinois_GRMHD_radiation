subroutine moving_puncture_bondi(ex,X,Y,Z,RP,mbh, &
     rho_star,tau,st_x,st_y,st_z,w,h,u0,rho_b,P,vx,vy,vz, &
     gam,Mdot,r_crit, &
     Kpoly_init,phi,lap, &
     xbh1,zbh1,bigP,Symmetry)
  
  implicit none
  
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z,RP
  real*8, dimension(ex(1),ex(2),ex(3))  :: rho_star
  real*8, dimension(ex(1),ex(2),ex(3))  :: st_x,st_y,st_z,w,h,u0,rho_b,P
  real*8, dimension(ex(1),ex(2),ex(3))  :: phi,lap,tau,vx,vy,vz
  real*8                                   :: Mdot, r_crit, gam,vvx,vvy,vvz
  real*8                                   :: kpoly_init,mbh,xbh1,zbh1,bigP
  integer                                  :: Symmetry
  interface
     REAL*8 FUNCTION find_rest_density(x1,x2,r,Gam,Mdot,rho_inf,mbh,KPOLY)
       implicit none
       REAL*8 x1,x2,r,Gam,Mdot,rho_inf,mbh,KPOLY
     end FUNCTION find_rest_density
  end interface
!
! Other variables:
!
  real*8                     :: u_x,u_y,u_z,u_t,alpha,mfac,psi6,bigP2
  integer                    :: i, j, k, kl, ku, count
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8                     :: HALF, ONE, PI, ZERO, TWO, FOUR, Tiny
  real*8                     :: rho0_i,u,upper,lower,fac,alph,ut,dk
  real*8                     :: r,u_r,st_r,eps,rhoeps,Cl,Cu,riso,uiso,psi
  real*8                     :: delta,a2,u2,rho_inf,fnr,KPOLY,E_EM,sb2,sbr
  real*8                     :: a2_inf,D1,D2,f,df,C_INF,Press,lambda,Br
  integer :: nr
  parameter(nr=10000)
  real*8, dimension(-nr:nr)  :: radius,rho0,smallM
  real *8 :: nx,ny,nz,gamv,ux,uy,uz
  real*8,dimension(4)        :: xa,ya
  real*8                     :: xi,err,temp,rhorm,rhorm2,rhorm3,drhorm,r2
  integer                    :: jj

  real*8                     :: Mdot0,sol,drsmall,xnew,ynew,znew,rfish

  real*8:: rr,radj,radj2,rmadj,rmadj2,msradj,msradj2,small,rfunc,urm

  parameter(HALF = 0.5D0, ONE = 1.0D0, ZERO = 0.0D0, TWO = 2.0D0, FOUR=4.0D0)
  parameter(Tiny = 1.0d-8, small=0.1d0)

  parameter (rr=0.5d0,drsmall=0.01d0)

  ! radj is the half the radius at which we alter the fields to account for the 
  ! singular velocity at the horizon.  Don't set it smaller than 0.5,
  ! which represents alterations only within r=M.  rmadj scales up by the mass
  radj=2.0*rr
  rmadj=mbh*radj

  PI = acos(-ONE)
  KPOLY = ONE
  !
  imin = lbound(rho_star,1)
  jmin = lbound(rho_star,2)
  kmin = lbound(rho_star,3)
  imax = ubound(rho_star,1)
  jmax = ubound(rho_star,2)
  kmax = ubound(rho_star,3)

  if(Mdot > 0) then

     write(6,*)'mdot:',mdot,' r_crit:',r_crit,' KPOLY:',KPOLY,' gam:',gam

     ! critical solution
     r = r_crit
     u2 = mbh/(TWO*r_crit)
     u = sqrt(u2)
     !Eq. G17
     a2 = u2/(ONE - 3.D0*u2)
     !Eq. G28
     rho0_i = ( a2/KPOLY/(Gam - (Gam/(Gam-ONE))*a2) )**(ONE/(Gam-ONE))
     !Eq. G33
     Mdot0 = FOUR*PI*u*rho0_i*r*r

     ! rescale Mdot to correct value
     KPOLY = (Mdot0/Mdot)**((Gam-ONE))
     rho0_i = ( a2/KPOLY/(Gam - (Gam/(Gam-ONE))*a2) )**(ONE/(Gam-ONE))
     Mdot = FOUR*PI*u*rho0_i*r*r

     write(*,*) "critical soln: ",rho0(1),a2,C_INF
     write(*,*) "KPOLY = ",KPOLY
     write(6,*)'new, old mdot:',Mdot,Mdot0

     kpoly_init = KPOLY
     !Eq. G30
     a2_inf = Gam-ONE - sqrt(ONE + 3.D0*a2)*(Gam-ONE - a2)
     !Eq. G28
     rho_inf = ( a2_inf/KPOLY/(Gam - (Gam/(Gam-ONE))*a2_inf) )**(ONE/(Gam-ONE))
     !Eq. G29, this is the integration constant
     C_INF = ( ONE + KPOLY*(Gam/(Gam-ONE))*(rho_inf**(Gam-ONE)) )**2
     fnr = nr
     dk = (0.1d0/r_crit)**(ONE/nr)
     !This is the density at the sonic point
     rho0(1) = rho0_i
     radius(1) = r
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
             * (ONE - TWO*mbh/r + (Mdot/(4.D0*PI*lower*r*r))**2) - C_INF
        Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mbh/r + (Mdot/(4.D0*PI*upper*r*r))**2) - C_INF

        !now bracket
        do while (Cl*Cu.gt.ZERO .and. upper-lower.gt.1.d-3*rho0(k+1))
           lower = sqrt(lower*rho0(k-1))
           upper = sqrt(upper*rho0(k-1))
           Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
                * (ONE - TWO*mbh/r + (Mdot/(4.D0*PI*lower*r*r))**2) - C_INF
           Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
                * (ONE - TWO*mbh/r + (Mdot/(4.D0*PI*upper*r*r))**2) - C_INF
        end do

        rho0(k) = find_rest_density(lower,upper,r,Gam,Mdot,rho_inf,mbh,KPOLY)
        rho0_i = rho0(k)
!!$     u2 = (Mdot/(4.D0*PI*rho0_i*r*r))**2
!!$     if (modulo(k,10)==0 .and. k.lt.200) write(*,*) k,Cl,Cu, &
!!$          ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(rho0_i**(Gam-ONE)) )**2 ) &
!!$          * (ONE - TWO*M/r + u2) - C_INF
     end do
     write(*,*) "critical soln: ",rho0(1),a2,C_INF
     write(*,*) "KPOLY = ",KPOLY
     r = r_crit
     do k = 0, -nr+1, -1
        r = r / dk
        lower = 0.8D0 * rho0(k+1)
        upper = 1.2D0 * rho0(k+1)
        radius(k) = r

        Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mbh/r + (Mdot/(4.D0*PI*lower*r*r))**2) - C_INF
        Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*mbh/r + (Mdot/(4.D0*PI*upper*r*r))**2) - C_INF

        do while (Cl*Cu.gt.ZERO .and. upper-lower.gt.1.d-4*rho0(k+1))
           lower = sqrt(lower*rho0(k+1))
           upper = sqrt(upper*rho0(k+1))
           Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
                * (ONE - TWO*mbh/r + (Mdot/(4.D0*PI*lower*r*r))**2) - C_INF
           Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
                * (ONE - TWO*mbh/r + (Mdot/(4.D0*PI*upper*r*r))**2) - C_INF
        end do

        rho0(k) = find_rest_density(lower,upper,r,Gam,Mdot,rho_inf,mbh,KPOLY)
        rho0_i = rho0(k)
        u2 = (Mdot/(4.D0*PI*rho0_i*r*r))**2
     end do
     write(*,*) "rho limits: ",rho0(-nr+1),rho0(nr)
     write(*,*) "r limits: ",radius(-nr+1),radius(nr)
  endif

!msadj is the schwarzschild radius of the adjustment surface
  msradj=rmadj*(1.0+0.5/radj)**2

  kl = floor( fnr*log(r_crit/msradj)/log(r_crit/0.1d0) ) + 1
  ku = kl + 1
  fac = (msradj - radius(kl))/(radius(ku) - radius(kl))
  if(fac.lt.ZERO .or. fac.gt.ONE) then
     write(*,*) "Error: fac = ",fac
     write(*,*) "    ", kl, ku, msradj, radius(kl), radius(ku)
  end if
  do jj=1,4
     xa(jj) = log( radius(jj-2+kl) )
     ya(jj) =  rho0(jj-2+kl)
  end do
  xi = log(msradj)
  !returns the density at the given radial point to rho0_i

  call polint(xa,ya,4,xi,rhorm,err)

  write(6,*)'rhorm:',rhorm,xa,ya,xi,msradj,rmadj,radj

! finite difference
  radj=radj+drsmall
  rmadj2=mbh*radj
  msradj2=rmadj2*(1.0+0.5/radj)**2
  kl = floor( fnr*log(r_crit/msradj2)/log(r_crit/0.1d0) ) + 1
  ku = kl + 1
  fac = (msradj2 - radius(kl))/(radius(ku) - radius(kl))
  if(fac.lt.ZERO .or. fac.gt.ONE) then
     write(*,*) "Error: fac = ",fac
     write(*,*) "    ", kl, ku, msradj2, radius(kl), radius(ku)
  end if
  do jj=1,4
     xa(jj) = log( radius(jj-2+kl) )
     ya(jj) =  rho0(jj-2+kl)
  end do
  xi = log(msradj2)
  !returns the density just outside the adjustment surface to rhorm2
  call polint(xa,ya,4,xi,rhorm2,err)

  write(6,*)'rhorm:',rhorm2,rhorm

  ! finite difference to get rhorm
  drhorm=(rhorm-rhorm2)/drsmall
  ! rhorm3 is the coefficient to smooth the density field inside the adjustment surface
  rhorm3=rhorm*(1.0+drhorm/rhorm*rmadj/4.0)
  urm=Mdot / (FOUR * PI * msradj * msradj * rhorm)

  write(6,*)'rhorm3:',rhorm3,rhorm,drhorm,rmadj,mdot,msradj

  ! fill up 3D arrays
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           rfish=sqrt(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2)
           xnew=x(i,j,k)/rfish*RP(i,j,k)
           ynew=y(i,j,k)/rfish*RP(i,j,k)
           znew=z(i,j,k)/rfish*RP(i,j,k)
           
           riso = sqrt((xnew-xbh1)**2 + ynew**2 + (znew-zbh1)**2)
           if(i.eq.2.and.j.eq.2.and.k.eq.2)write(6,*)'riso:',riso,x(i,j,k)

           nx=(xnew-xbh1)/riso
           ny=ynew/riso
           nz=(znew-zbh1)/riso

           ! now we switch to schwarzschild
           r= riso*(1.0d0+mbh/2.0d0/riso)**2

!!$           if(r.gt.0.1d0 .and. r.lt.9.d0) then
!           if(i.eq.1.and.j.eq.2.and.k.eq.1)write(6,*)'loop:',riso,r,mbh,Mdot
           if(riso.gt.rmadj.and.Mdot.gt.0.0) then
              kl = floor( fnr*log(r_crit/r)/log(r_crit/0.1d0) ) + 1
              ku = kl + 1
              fac = (r - radius(kl))/(radius(ku) - radius(kl))
              if(fac.lt.ZERO .or. fac.gt.ONE) then
                 write(*,*) "Error: fac = ",fac
                 write(*,*) "    ", kl, ku, r, radius(kl), radius(ku)
              end if
!!$              rho0_i = fac * rho0(ku) + (ONE-fac) * rho0(kl)
              !
              ! Perform 3rd order interpolation
              !
              do jj=1,4
                 xa(jj) = log( radius(jj-2+kl) )
                 ya(jj) =  rho0(jj-2+kl)
              end do
              xi = log(r)
              !returns the density at the given radial point to rho0_i
              call polint(xa,ya,4,xi,rho0_i,err)
              !
              rhoeps = KPOLY * (rho0_i**Gam) / (Gam - ONE)
              eps = rhoeps / rho0_i
              h(i,j,k) = ONE + Gam * eps
              ! u is the value in Schwarzschild coords.
              u = Mdot / (FOUR * PI * r * r * rho0_i)
              !uiso from my notes
              uiso= u/(1.0d0+mbh/2.0d0/riso)/(1.0d0-mbh/2.0d0/riso)
              psi = exp(phi(i,j,k))
              alph = lap(i,j,k)+1.0
              ut = sqrt(1+uiso**2*psi**4)/alph
              gamv=1.0/sqrt(1.0-bigP**2)
              sol=alph/psi**2
              if(Symmetry==4) then
                 vvz=(sol*bigP-uiso*nz/ut)/(1.0-bigP*uiso*nz/ut/sol)
                 vvy=-1.0*uiso*ny/ut/gamv/(1.0-bigP*uiso*nz/ut/sol)
                 vvx=-1.0*uiso*nx/ut/gamv/(1.0-bigP*uiso*nz/ut/sol)
              else
                 vvx=(sol*bigP-uiso*nx/ut)/(1.0-bigP*uiso*nx/ut/sol)
                 vvy=-1.0*uiso*ny/ut/gamv/(1.0-bigP*uiso*nx/ut/sol)
                 vvz=-1.0*uiso*nz/ut/gamv/(1.0-bigP*uiso*nx/ut/sol)
              endif
              ut=(alph**2-psi**4*(vvx**2+vvy**2+vvz**2))**-0.5
              ux=vvx*ut
              uy=vvy*ut
              uz=vvz*ut

              u0(i,j,k)=ut
              rho_b(i,j,k) = rho0_i
              rho_star(i,j,k) = rho0_i * alph * ut * psi**6
              w(i,j,k) = rho_star(i,j,k) * alph * ut
              Press = (Gam - ONE)*rhoeps
              P(i,j,k) = Press
              tau(i,j,k)=w(i,j,k)*h(i,j,k)-psi**6*Press-rho_star(i,j,k)
              u_t=-1.0d0*alph**2*ut
              u_x=psi**4*ux
              u_y=psi**4*uy
              u_z=psi**4*uz
              st_x(i,j,k) = rho_star(i,j,k) * h(i,j,k) *u_x
              st_y(i,j,k) = rho_star(i,j,k) * h(i,j,k) *u_y
              st_z(i,j,k) = rho_star(i,j,k) * h(i,j,k) *u_z
              vx(i,j,k)=vvx
              vy(i,j,k)=vvy
              vz(i,j,k)=vvz

!              if(i.eq.60.and.j.eq.1.and.k.eq.1)write(6,*)'boundinit:',alph,psi,sol,bigP,uiso,nx,ut,vvx,vvy,vvz

            else if (riso.le.rmadj.and.riso.ge.rmadj/2.0.and.mdot.gt.0) then
              psi=exp(phi(i,j,k))
              rfunc=1.0-drhorm/rhorm*((rmadj-riso)-(rmadj**2-riso**2)/rmadj)

              rhorm2=rhorm*rfunc
              r2=msradj*rfunc

              rhoeps = KPOLY * (rhorm2**Gam) / (Gam - ONE)
              eps = rhoeps / rhorm2
              h(i,j,k) = ONE + Gam * eps
              alph = lap(i,j,k)+1.0
              sol=alph/psi**2
              u = urm*riso/rmadj
              uiso=u/(1.0-0.25*mbh**2/rmadj**2)
              ut = sqrt(1+uiso**2*psi**4)/alph
              bigP2=bigP*riso/rmadj
              gamv=1.0/sqrt(1.0-bigP2**2)
              if(Symmetry==4) then
                 vvz=(sol*bigP2-uiso*nz/ut)/(1.0-bigP2*uiso*nz/ut/sol)
                 vvy=-1.0*uiso*ny/ut/gamv/(1.0-bigP2*uiso*nz/ut/sol)
                 vvx=-1.0*uiso*nx/ut/gamv/(1.0-bigP2*uiso*nz/ut/sol)
              else
                 vvx=(sol*bigP2-uiso*nx/ut)/(1.0-bigP2*uiso*nx/ut/sol)
                 vvy=-1.0*uiso*ny/ut/gamv/(1.0-bigP2*uiso*nx/ut/sol)
                 vvz=-1.0*uiso*nz/ut/gamv/(1.0-bigP2*uiso*nx/ut/sol)
              endif
              ut=(alph**2-psi**4*(vvx**2+vvy**2+vvz**2))**-0.5
              ux=vvx*ut
              uy=vvy*ut
              uz=vvz*ut

              u0(i,j,k)=ut
              rho_b(i,j,k) = rhorm2
              rho_star(i,j,k) = rhorm2 * alph * ut * psi**6
              w(i,j,k) = rho_star(i,j,k) * alph * ut
              Press = (Gam - ONE)*rhoeps
              P(i,j,k) = Press
              u_t=-1.0d0*alph**2*ut
              u_x=psi**4*ux
              u_y=psi**4*uy
              u_z=psi**4*uz
              st_x(i,j,k) = rho_star(i,j,k) * h(i,j,k) *u_x
              st_y(i,j,k) = rho_star(i,j,k) * h(i,j,k) *u_y
              st_z(i,j,k) = rho_star(i,j,k) * h(i,j,k) *u_z
              vx(i,j,k)=vvx
              vy(i,j,k)=vvy
              vz(i,j,k)=vvz
              tau(i,j,k)=w(i,j,k)*h(i,j,k)-rho_star(i,j,k)-Press*psi**6

          else if (riso.le.rmadj/2.0.and.mdot.gt.0) then
              psi=exp(phi(i,j,k))

              rfunc=0.5*(1.0+small-(1.0-small)*cos(2.0*3.141592653*riso/rmadj))

              rhorm2=rhorm3*rfunc
              r2=msradj*rfunc
              rhoeps = KPOLY * (rhorm2**Gam) / (Gam - ONE)
              eps = rhoeps / rhorm2
              h(i,j,k) = ONE + Gam * eps
              alph = lap(i,j,k)+1.0
              sol=alph/psi**2
              if(i.eq.2.and.j.eq.2.and.k.eq.2)write(6,*)'a:',psi,rfunc,rhorm2,r2,rhoeps,eps,alph,sol
              u = urm*riso/rmadj
              uiso=u/(1.0-0.25*mbh**2/rmadj**2)
              ut = sqrt(1+uiso**2*psi**4)/alph
              bigP2=bigP*riso/rmadj
              gamv=1.0/sqrt(1.0-bigP2**2)
              if(Symmetry==4) then
                 vvz=(sol*bigP2-uiso*nz/ut)/(1.0-bigP2*uiso*nz/ut/sol)
                 vvy=-1.0*uiso*ny/ut/gamv/(1.0-bigP2*uiso*nz/ut/sol)
                 vvx=-1.0*uiso*nx/ut/gamv/(1.0-bigP2*uiso*nz/ut/sol)
              else
                 vvx=(sol*bigP2-uiso*nx/ut)/(1.0-bigP2*uiso*nx/ut/sol)
                 vvy=-1.0*uiso*ny/ut/gamv/(1.0-bigP2*uiso*nx/ut/sol)
                 vvz=-1.0*uiso*nz/ut/gamv/(1.0-bigP2*uiso*nx/ut/sol)
              endif

              ut=(alph**2-psi**4*(vvx**2+vvy**2+vvz**2))**-0.5
              ux=vvx*ut
              uy=vvy*ut
              uz=vvz*ut

              u0(i,j,k)=ut
              rho_b(i,j,k) = rhorm2
              rho_star(i,j,k) = rhorm2 * alph * ut * psi**6
              w(i,j,k) = rho_star(i,j,k) * alph * ut
              Press = (Gam - ONE)*rhoeps
              P(i,j,k) = Press
              u_t=-1.0d0*alph**2*ut
              u_x=psi**4*ux
              u_y=psi**4*uy
              u_z=psi**4*uz
              st_x(i,j,k) = rho_star(i,j,k) * h(i,j,k) *u_x
              st_y(i,j,k) = rho_star(i,j,k) * h(i,j,k) *u_y
              st_z(i,j,k) = rho_star(i,j,k) * h(i,j,k) *u_z
              vx(i,j,k)=vvx
              vy(i,j,k)=vvy
              vz(i,j,k)=vvz
              tau(i,j,k)=w(i,j,k)*h(i,j,k)-rho_star(i,j,k)-Press*psi**6
           else
              kpoly_init=1.0
              rho_star(i,j,k) = ZERO
              st_x(i,j,k)     = ZERO
              st_y(i,j,k)     = ZERO
              st_z(i,j,k)     = ZERO
              w(i,j,k)        = ZERO
              h(i,j,k)     = ONE
              u0(i,j,k)= ONE
              rho_b(i,j,k)     = ZERO
              P(i,j,k)  = ZERO
              vx(i,j,k)=ZERO
              vy(i,j,k)=ZERO
              vz(i,j,k)=ZERO
              tau(i,j,k)=ZERO
            end if
        end do
     end do
  end do

  write(6,*)'vxs:',st_x(2,2,2),st_x(3,2,2),st_x(4,2,2),st_x(3,2,3)
  write(6,*)'phis:',phi(2,2,2),phi(3,2,2),phi(4,2,2),phi(3,2,3)
  write(6,*)'rhos:',rho_star(2,2,2),rho_star(3,2,2),rho_star(4,2,2),rho_star(3,2,3)
  write(6,*)'rhob:',rho_b(2,2,2),rho_b(3,2,2),rho_b(4,2,2),rho_b(3,2,3)

end subroutine moving_puncture_bondi

REAL*8 FUNCTION find_rest_density(x1,x2,r,Gamma,Mdot,rho_inf,MBH,KPOLY)
  implicit none
  INTEGER MAXIT
  REAL*8 x1,x2,r,Gamma,Mdot,rho_inf,MBH,KPOLY
  REAL*8 xacc, scale
  logical success
  !    EXTERNAL funcd
  PARAMETER (MAXIT=100)
  INTEGER j
  REAL*8 df,dx,dxold,f,fh,fl,temp,xh,xl
  xacc = 1.D-16
  success=.true.
  call transf(x1,fl,df,r,Gamma,Mdot,rho_inf,MBH,KPOLY)
  call transf(x2,fh,df,r,Gamma,Mdot,rho_inf,MBH,KPOLY)
  if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) then
     find_rest_density = x2 - 2.D0*fh/df
     call transf(find_rest_density,f,df,r,Gamma,Mdot,rho_inf,MBH,KPOLY)
     if(abs(f).lt.xacc)then
        return
     else
        !          write(*,*) "bracket: ",r,fl,fh
        !       x2 = 1.1 * x1
     end if
  end if
  if(abs(fl).lt.xacc)then
     find_rest_density=x1
     return
  else if(abs(fh).lt.xacc)then
     find_rest_density=x2
     return
  else if(fl.lt.0.)then
     xl=x1
     xh=x2
  else
     xh=x1
     xl=x2
  endif
  find_rest_density=.5*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  call transf(find_rest_density,f,df,r,Gamma,Mdot,rho_inf,MBH,KPOLY)
  do j=1,MAXIT
     if(((find_rest_density-xh)*df-f)*((find_rest_density-xl)*df-f).ge.0. &
          .or. abs(2.*f).gt.abs(dxold*df) ) then
        dxold=dx
        dx=0.5D0*(xh-xl)
        find_rest_density=xl+dx
        if(xl.eq.find_rest_density) return
     else
        dxold=dx
        dx=f/df
        temp=find_rest_density
        find_rest_density=find_rest_density-dx
        if(temp.eq.find_rest_density) return
     endif
     if(abs(dx).lt.xacc) return
     call transf(find_rest_density,f,df,r,Gamma,Mdot,rho_inf,MBH,KPOLY)
     if(f.lt.0.) then
        xl=find_rest_density
     else
        xh=find_rest_density
     endif
  end do
  success=.false.
  return
END FUNCTION find_rest_density

subroutine transf(rho0,f,df,r,Gamma,Mdot,rho_inf,MBH,KPOLY)
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
end subroutine transf
