!-----------------------------------------------------------------------------
!
! $Id
!
!-----------------------------------------------------------------------------
!
! Set up matter source terms (spherical accretion onto a black hole)
!
!-----------------------------------------------------------------------------
subroutine mag_bondi_accreting_matter(ext, x, y, z, &
     neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, gamma_th, &
     rho_star, tau, st_x, st_y, st_z, mhd_st_x,mhd_st_y,mhd_st_z, &
     w, enth, rhob, P,u0,vx,vy,vz, n, r_crit, constrained_transport_scheme, &
     Ax,Ay,Az, Bx_stagger,By_stagger,Bz_stagger, &
     Bx,By,Bz,sbt,sbx,sby,sbz,Ex,Ey,Ez, Bxtilde,Bytilde,Bztilde, &
     Bstrength,r0, puncture_id, excision_enable, excision_radius,xmax,zmax)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ext
  integer                                   :: constrained_transport_scheme
  integer, parameter                       :: nr = 100
  real*8, dimension(ext(1),ext(2),ext(3))  :: x,y,z,rho_star,tau
  real*8, dimension(ext(1),ext(2),ext(3))  :: st_x,st_y,st_z,w,enth,rhob,P
  real*8, dimension(ext(1),ext(2),ext(3))  :: Bx,By,Bz, u0,vx,vy,vz
  real*8, dimension(ext(1),ext(2),ext(3))  :: Bxtilde,Bytilde,Bztilde
  real*8, dimension(ext(1),ext(2),ext(3))  :: sbt,sbx,sby,sbz
  real*8, dimension(ext(1),ext(2),ext(3))  :: Ex,Ey,Ez,Ax,Ay,Az
  real*8, dimension(ext(1),ext(2),ext(3))  :: mhd_st_x,mhd_st_y,mhd_st_z
  real*8, dimension(ext(1),ext(2),ext(3))  :: Bx_stagger,By_stagger,Bz_stagger
  real*8                                   :: n, Mdot, r_crit, dr,r0
  real*8                                   :: kpoly_init,Bstrength,gamma_th
  real*8                                   :: zmax,xmax,x0
  integer                                  :: neos
  real*8, dimension(10)                    :: rho_tab,P_tab,eps_tab
  real*8, dimension(11)                    :: k_tab,gamma_tab
  integer                                    :: excision_enable, puncture_id
  real*8                                    :: excision_radius,rmin
  interface
     REAL*8 FUNCTION find_rest_density(x1,x2,r,Gam,Mdot,rho_inf,M,KPOLY)
       implicit none
       REAL*8 x1,x2,r,Gam,Mdot,rho_inf,M,KPOLY
     end FUNCTION find_rest_density
  end interface
!
! Other variables:
!
  real*8                     :: u_x,u_y,u_z,u_t,alpha,mfac,psi6,Gam,psi6bar
  integer                    :: i, j, k, kl, ku, count, n_zout,n_xout,n_x0
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8                     :: HALF, ONE, PI, ZERO, TWO, FOUR, Tiny, delta
  real*8                     :: rho0_i,u,upper,lower,fac,alph,ut,dk,vr
  real*8                     :: r,rbar,u_r,st_r,H,eps,rhoeps,Cl,Cu
  real*8                     :: a2,u2,rho_inf,fnr,M,KPOLY,E_EM,sb2,sbr
  real*8                     :: a2_inf,D1,D2,f,df,C_INF,Press,lambda,Br
  real*8                     :: Ax0,Ay0,dX,dY,dZ,r_stagger,amp
  real*8                     :: x_stagger,y_stagger,z_stagger
  real*8                     :: sb_r,sb_th,sb_ph,sb_x,sb_y,sb_z,pm,sb_t
  real*8                     :: xup,yup,zup,zout,xout,xl,yl,zl
  real*8                     :: rbarxu, rbaryu,rbarzu
  real*8                     :: psi6bar_xs,psi6bar_ys,psi6bar_zs
  real*8                     :: Bxts,Byts,Bzts,Bxts0,Byts0,Bzts1,Byts1,Bxts1
  real*8                     :: b2oB02
  real*8, dimension(-nr:nr)  :: radius,rho0,smallM
  real*8,dimension(4)              :: xa,ya
  real*8                     :: xi,err,temp
  integer                    :: jj
  parameter(HALF = 0.5D0, ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0, FOUR = 4.D0)
  parameter(Tiny = 1.e-8)
!
  PI = acos(-ONE)
  Gam = ONE + ONE/n
  M = ONE
  KPOLY = ONE
!
  imin = lbound(rho_star,1)
  jmin = lbound(rho_star,2)
  kmin = lbound(rho_star,3)
  imax = ubound(rho_star,1)
  jmax = ubound(rho_star,2)
  kmax = ubound(rho_star,3)
  dX = x(2,1,1)-x(1,1,1)
  dY = y(1,2,1)-y(1,1,1)
  dZ = z(1,1,2)-z(1,1,1)
  delta = 0.1d0*min( dX, dY, dZ ) * r_crit
! critical solution
  r = r_crit
  u2 = M/(TWO*r_crit)
  u = sqrt(u2)
  a2 = u2/(ONE - 3.D0*u2)
  rho0_i = ( a2/KPOLY/(Gam - (Gam/(Gam-ONE))*a2) )**(ONE/(Gam-ONE))
  Mdot = FOUR*PI*u*rho0_i*r*r
! make Mdot = 1
  KPOLY = Mdot**(Gam-ONE)
  rho0_i = ( a2/KPOLY/(Gam - (Gam/(Gam-ONE))*a2) )**(ONE/(Gam-ONE))
  Mdot = FOUR*PI*u*rho0_i*r*r
  ! Reset K_init and the EOS table
  kpoly_init = KPOLY
  do i=1,2
     k_tab(i) = kpoly_init
     gamma_tab(i) = Gam
  end do
  neos = 1
  gamma_th = Gam
  rho_tab(1)=1.d0
  P_tab(1)=kpoly_init
  eps_tab(1)=n*P_tab(1)/rho_tab(1)
  a2_inf = Gam-ONE - sqrt(ONE + 3.D0*a2)*(Gam-ONE - a2)
  rho_inf = ( a2_inf/KPOLY/(Gam - (Gam/(Gam-ONE))*a2_inf) )**(ONE/(Gam-ONE))
  C_INF = ( ONE + KPOLY*(Gam/(Gam-ONE))*(rho_inf**(Gam-ONE)) )**2
  fnr = nr
  dr = 0.05d0
  dk = (dr/r_crit)**(ONE/nr)
  rho0(1) = rho0_i
  radius(1) = r
! fill up 1D arrays
  do k = 2, nr
     r = r * dk
     lower = 0.9D0 * rho0(k-1)
     upper = 1.1D0 * rho0(k-1)
     radius(k) = r
     Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*M/r + (Mdot/(4.D0*PI*lower*r*r))**2) - C_INF
     Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*M/r + (Mdot/(4.D0*PI*upper*r*r))**2) - C_INF
     do while (Cl*Cu.gt.ZERO .and. upper-lower.gt.1.e-3*rho0(k+1))
        lower = sqrt(lower*rho0(k-1))
        upper = sqrt(upper*rho0(k-1))
        Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*M/r + (Mdot/(4.D0*PI*lower*r*r))**2) - C_INF
        Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*M/r + (Mdot/(4.D0*PI*upper*r*r))**2) - C_INF
     end do
     rho0(k) = find_rest_density(lower,upper,r,Gam,Mdot,rho_inf,M,KPOLY)
     rho0_i = rho0(k)
  end do
  !write(*,*) critical soln: ,rho0(1),a2,C_INF
  !write(*,*) KPOLY = ,KPOLY
  r = r_crit
  do k = 0, -nr+1, -1
     r = r / dk
     lower = 0.8D0 * rho0(k+1)
     upper = 1.2D0 * rho0(k+1)
     radius(k) = r
     Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*M/r + (Mdot/(4.D0*PI*lower*r*r))**2) - C_INF
     Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
          * (ONE - TWO*M/r + (Mdot/(4.D0*PI*upper*r*r))**2) - C_INF
     do while (Cl*Cu.gt.ZERO .and. upper-lower.gt.1.e-4*rho0(k+1))
        lower = sqrt(lower*rho0(k+1))
        upper = sqrt(upper*rho0(k+1))
        Cl = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(lower**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*M/r + (Mdot/(4.D0*PI*lower*r*r))**2) - C_INF
        Cu = ( ( ONE + KPOLY*(Gam/(Gam-ONE))*(upper**(Gam-ONE)) )**2 ) &
             * (ONE - TWO*M/r + (Mdot/(4.D0*PI*upper*r*r))**2) - C_INF
     end do
     rho0(k) = find_rest_density(lower,upper,r,Gam,Mdot,rho_inf,M,KPOLY)
     rho0_i = rho0(k)
     u2 = (Mdot/(4.D0*PI*rho0_i*r*r))**2
  end do
  write(*,*) "rho limits: ",rho0(-nr+1),rho0(nr)
  write(*,*) "r limits: ",radius(-nr+1),radius(nr)
  ! Calculate b^2/B0^2 at r=2M
  ! Interpolate to find rho_0 at r=2M
  r = 2.d0*M
  kl = floor( fnr*log(r_crit/r)/log(r_crit/dr) ) + 1
  ku = kl + 1
  fac = (r - radius(kl))/(radius(ku) - radius(kl))
  do jj=1,4
     xa(jj) = log( radius(jj-2+kl) )
     ya(jj) =  rho0(jj-2+kl)
  end do
  xi = log(r)
  call polint(xa,ya,4,xi,rho0_i,err)
  u = Mdot / (FOUR * PI * r * r * rho0_i)
  H = TWO*M / r
  psi6 = sqrt(ONE + H)
  ! ut = (u*H - sqrt(ONE + u*u - H))/(-ONE + H) ! This formula doesnt 
  ! work at r=2M, cause you get 0/0, need to take limit.
  ut = u + 0.5d0/u
  ! B^r = B0/[ r^2 sqrt(1+2M/r)] for Kerr-Schild metric
  ! Hence B^2 = B0^2/r^4. 
  ! u_\mu B^\mu = u_r B^r  
  ! b^2 = [ B^2 + (u\mu B^\mu)^2 ] / [4 pi (alpha ut)^2] , alpha = 1/sqrt(1+2M/r) for 
  !            Kerr-Schild metric
  u_r = H * ut - (1.d0+H)*u
  b2oB02 = (1.d0/r**4 + u_r*u_r / r**4 ) / (4.d0*PI) / (ut/psi6)**2
  ! Determine B0 (amp here) so that b^2/rho_0 = Bstrength at r=2M
  amp = sqrt(Bstrength * rho0_i / b2oB02)
  !amp = sqrt(abs(Bstrength*rho0(nr)))
  Press = KPOLY*rho0_i**Gam
  write(*,*) ' '
  write(*,*) 'b^2/rho0 at r=2M: ',amp*amp*b2oB02/rho0_i
  write(*,*) 'b^2/P at r=2M: ',amp*amp*b2oB02/Press
! fill up 3D arrays
  rmin = dr
  if (Bstrength .lt. 0.d0) amp = -amp
  if (excision_enable==1) rmin = 0.8d0*excision_radius+r0
  do k = kmin,kmax
     do j = jmin,jmax
        do i = imin,imax
           xl = x(i,j,k)
           yl = y(i,j,k)
              zl = z(i,j,k)
            rbar = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
           if (puncture_id==1) then
              r = rbar*(1.d0+0.5d0/rbar)**2
           else
              r = rbar + r0
               end if
           if(r.gt.rmin) then
              kl = floor( fnr*log(r_crit/r)/log(r_crit/dr) ) + 1
              ku = kl + 1
              fac = (r - radius(kl))/(radius(ku) - radius(kl))
              if(fac.lt.ZERO .or. fac.gt.ONE) then
                 write(*,*) "Error: fac = ",fac
                 write(*,*) "    ", kl, ku, r, radius(kl), radius(ku)
                 stop
              end if
!
! Perform 3rd order interpolation 
! 
              do jj=1,4
                 xa(jj) = log( radius(jj-2+kl) )
                 ya(jj) =  rho0(jj-2+kl)
              end do
              xi = log(r)
              call polint(xa,ya,4,xi,rho0_i,err)
!
              rhoeps = KPOLY * (rho0_i**Gam) / (Gam - ONE)
              eps = rhoeps / rho0_i
              enth(i,j,k) = ONE + Gam * eps
              u = Mdot / (FOUR * PI * r * r * rho0_i)
              H = TWO*M / r
              psi6 = sqrt(ONE + H)
              if (puncture_id==1) then
                 psi6bar = (1.d0+0.5d0/rbar)**6
              else
                 psi6bar = psi6 * (r/rbar)**2
                end if
              if (ONE + u*u - H .gt. 0.d0 .and. H .lt. 1.d0) then
                 ut = (u*H - sqrt(ONE + u*u - H))/(-ONE + H)
                 vr = -u/ut
              else
                 ut = psi6
                 vr = -H/(1.d0+H)
              end if
              alph = ONE/psi6
              if (puncture_id==1) then
                 ! Assume puncture lapse alph = 1/psi^2
                 alph = 1.d0/(1.d0+0.5d0/rbar)**2
                 u = u/(1.d0-0.25d0/rbar**2)
                 ut = sqrt(1.d0+u*u*(1.d0+0.5d0/rbar)**4) / alph
                  ! Limit velocity: alpha ut <= 60
                 if (alph*ut .gt. 60.d0) then
                     u = sqrt(3599.d0)/(1.d0+0.5d0/rbar)**2
                    ut = 60.d0/alph
                  end if
                 vr = -u/ut
              end if
              rhob(i,j,k) = rho0_i
               u0(i,j,k) = ut
              vx(i,j,k) = vr * x(i,j,k)/rbar
              vy(i,j,k) = vr * y(i,j,k)/rbar
              vz(i,j,k) = vr * z(i,j,k)/rbar
              rho_star(i,j,k) = rho0_i * alph * ut * psi6bar
              w(i,j,k) = rho_star(i,j,k) * alph * ut
              Press = (Gam - ONE)*rhoeps
              P(i,j,k) = Press
              tau(i,j,k) = w(i,j,k)*enth(i,j,k) - psi6bar*Press &
                   - rho_star(i,j,k)
              u_r = -(ONE/(ONE-H))*(u - H * sqrt(ONE + u*u - H))
              u_t = -sqrt(ONE + u*u - H)
              if (puncture_id==1) then
                 u_r = ut*vr*(1.d0+0.5d0/rbar)**4
                  u_t = -alph*alph*ut
              end if
              u_x = u_r * x(i,j,k)/rbar
              u_y = u_r * y(i,j,k)/rbar
              u_z = u_r * z(i,j,k)/rbar
              st_r = rho_star(i,j,k) * enth(i,j,k) * u_r
              st_x(i,j,k) = st_r * x(i,j,k)/rbar
              st_y(i,j,k) = st_r * y(i,j,k)/rbar
              st_z(i,j,k) = st_r * z(i,j,k)/rbar
              Br = amp/psi6/r**2
              E_EM = (1+H)*Br**2/(8.d0*PI)
              if (puncture_id==1) then
                 Br = amp/psi6bar/rbar**2
              end if
              Bx(i,j,k) = Br * x(i,j,k)/rbar
              By(i,j,k) = Br * y(i,j,k)/rbar
              Bz(i,j,k) = Br * z(i,j,k)/rbar
                Bxtilde(i,j,k) = Bx(i,j,k) * psi6bar
              Bytilde(i,j,k) = By(i,j,k) * psi6bar
              Bztilde(i,j,k) = Bz(i,j,k) * psi6bar
              sbt(i,j,k) = u_r*Br/alph/sqrt(4.d0*PI)
              sbr = Br/(alph*ut)/sqrt(4.d0*PI) + sbt(i,j,k)*vr
              sb2 = -(ONE-H)*sbt(i,j,k)**2 + TWO*H*sbt(i,j,k)*sbr &
                   + (ONE+H)*sbr**2
              if (puncture_id==1) then
                 sb2 = sbr*sbr*(1.d0+0.5d0/rbar)**4 - (alph*sbt(i,j,k))**2
              end if
              sbx(i,j,k) = sbr * x(i,j,k)/rbar
              sby(i,j,k) = sbr * y(i,j,k)/rbar
              sbz(i,j,k) = sbr * z(i,j,k)/rbar
              Ex(i,j,k) = ZERO
              Ey(i,j,k) = ZERO
              Ez(i,j,k) = ZERO
              mhd_st_x(i,j,k) = st_x(i,j,k)
              mhd_st_y(i,j,k) = st_y(i,j,k)
              mhd_st_z(i,j,k) = st_z(i,j,k)
              if (constrained_transport_scheme==3) then
                 call compute_Ax_mag_bondi(xl,yl+0.5d0*dY,zl+0.5d0*dZ,amp,Ax(i,j,k))
                 call compute_Ay_mag_bondi(xl+0.5d0*dX,yl,zl+0.5d0*dZ,amp,Ay(i,j,k))
                 Az(i,j,k) = 0.d0
                 call compute_Bi_mag_bondi(xl,yl,zl,dX,dY,dZ,amp, &
                                Bx_stagger(i,j,k),By_stagger(i,j,k),Bz_stagger(i,j,k),&
                                Bx(i,j,k),By(i,j,k),Bz(i,j,k),r0,puncture_id)
                 Br = (Bx(i,j,k)*xl + By(i,j,k)*yl + Bz(i,j,k)*zl)/rbar
                 sbt(i,j,k) = u_r*Br/alph/sqrt(4.d0*PI)
                 sbr = Br/(alph*ut)/sqrt(4.d0*PI) + sbt(i,j,k)*vr
                 sbx(i,j,k) = Bx(i,j,k)/(alph*ut)/sqrt(4.d0*PI) + &
                                sbt(i,j,k)*vx(i,j,k)
                 sby(i,j,k) = By(i,j,k)/(alph*ut)/sqrt(4.d0*PI) + &
                                sbt(i,j,k)*vy(i,j,k)
                 sbz(i,j,k) = Bz(i,j,k)/(alph*ut)/sqrt(4.d0*PI) + &
                                sbt(i,j,k)*vz(i,j,k)
                 ! sb_theta = gamma_{theta theta} sb^theta ,
                 ! sb_phi =  gamma_{phi phi} sb^phi in both coordinate systems.
                 pm = sqrt(xl**2+yl**2)
                 sb_th = r*r/rbar/rbar*( xl*zl/pm*sbx(i,j,k) + &
                        yl*zl/pm*sby(i,j,k) - pm*sbz(i,j,k) )
                 sb_ph = r*r/rbar/rbar*(xl*sby(i,j,k)-yl*sbx(i,j,k))
                 if (puncture_id==1) then
                    sb_r = sbr*(1.d0+0.5d0/rbar)**4
                    sb_t = -sbt(i,j,k)*alph**2
                 else
                    sb_r = sbr*(1.d0+2.d0/r) + 2.d0/r*sbt(i,j,k)
                    sb_t = -(1.d0-2.d0/r)*sbt(i,j,k) + 2.d0/r*sbr
                 end if
                 sb_x = xl/rbar*sb_r + xl*zl/pm/rbar/rbar*sb_th &
                        - yl/pm/pm*sb_ph
                 sb_y = yl/rbar*sb_r + yl*zl/pm/rbar/rbar*sb_th &
                        + xl/pm/pm*sb_ph
                 sb_z = zl/rbar*sb_r - pm/rbar/rbar*sb_th
                 sb2 = sbt(i,j,k)*sb_t + sbx(i,j,k)*sb_x + sby(i,j,k)*sb_y + sbz(i,j,k)*sb_z
                 mhd_st_x(i,j,k) = st_x(i,j,k) + alph*psi6bar*( sb2*ut*u_x - sbt(i,j,k)*sb_x )
                 mhd_st_y(i,j,k) = st_y(i,j,k) + alph*psi6bar*( sb2*ut*u_y - sbt(i,j,k)*sb_y )
                 mhd_st_z(i,j,k) = st_z(i,j,k) + alph*psi6bar*( sb2*ut*u_z - sbt(i,j,k)*sb_z )
              end if
              tau(i,j,k) = tau(i,j,k) + psi6bar * &
                   (-HALF*sb2 - (alph*sbt(i,j,k))**2 + sb2*(alph*ut)**2)
              sbt(i,j,k) = sbt(i,j,k)/alph
              sbx(i,j,k) = sbx(i,j,k)/alph
              sby(i,j,k) = sby(i,j,k)/alph
              sbz(i,j,k) = sbz(i,j,k)/alph
           else
               H = TWO*M / r
              psi6 = sqrt(ONE + H)
              vr = -H/(1.d0+H)
              rho_star(i,j,k) = ZERO
              tau(i,j,k)   = ZERO
              st_x(i,j,k)     = ZERO
              st_y(i,j,k)     = ZERO
              st_z(i,j,k)     = ZERO
              w(i,j,k)        = ZERO
              enth(i,j,k)     = ONE
                P(i,j,k) = ZERO
              rhob(i,j,k) = ZERO
              Bx(i,j,k) = ZERO
              By(i,j,k) = ZERO
              Bz(i,j,k) = ZERO
              Bxtilde(i,j,k) = ZERO
              Bytilde(i,j,k) = ZERO
              Bztilde(i,j,k) = ZERO
              Ex(i,j,k) = ZERO
              Ey(i,j,k) = ZERO
              Ez(i,j,k) = ZERO
              sbt(i,j,k) = ZERO
              sbx(i,j,k) = ZERO
              sby(i,j,k) = ZERO
              sbz(i,j,k) = ZERO
              u0(i,j,k) = psi6
              vx(i,j,k) = vr*x(i,j,k)/rbar
              vy(i,j,k) = vr*y(i,j,k)/rbar
              vz(i,j,k) = vr*z(i,j,k)/rbar
              mhd_st_x(i,j,k) = ZERO
              mhd_st_y(i,j,k) = ZERO
              mhd_st_z(i,j,k) = ZERO
           end if
        end do
     end do
  end do
!!!  ! Setup magnetic field initial data for constrained_transport_scheme=3
!!!  if (constrained_transport_scheme==3) then 
!!!
!!!     x0 = x(imax,1,1)+0.5d0*dX - dX*int(x(imax,1,1)/dX + 0.5d0)
!!!     n_xout = int( (1.1d0*xmax-x(imax,1,1)-0.5d0*dX)/dX + 0.5d0 )
!!!     n_zout = int( (1.1d0*zmax-z(1,1,kmax)-0.5d0*dZ)/dZ + 0.5d0 )
!!!     xout = x(imax,1,1) + 0.5d0*dX + n_xout*dX
!!!     zout = z(1,1,kmax) + 0.5d0*dZ + n_zout*dZ
!!!     Bzts1 = 0.d0
!!!     do j=jmin,jmax
!!!        do i=imax,imin,-1
!!!           do k=kmax,kmin,-1
!!!              xl = x(i,j,k)
!!!              yl = y(i,j,k)
!!!              zl = z(i,j,k)
!!!              x_stagger = xl + 0.5d0*dX
!!!              y_stagger = yl + 0.5d0*dY
!!!              z_stagger = zl + 0.5d0*dZ
!!!              xup = xl + dX
!!!              yup = yl + dY
!!!              zup = zl + dZ
!!!              rbar = sqrt(xl**2 + yl**2 + zl**2)
!!!              rbarxu = sqrt(xup**2 + yl**2 + zl**2)
!!!              rbaryu = sqrt(xl**2 + yup**2 + zl**2)
!!!              rbarzu = sqrt(xl**2 + yl**2 + zup**2)
!!!              ! Calculate psi^6 on staggered grid
!!!              if (puncture_id==1) then 
!!!                 psi6bar_xs = ((1.d0+0.5d0/rbar)**3)*((1.d0+0.5d0/rbarxu)**3)
!!!                 psi6bar_ys = ((1.d0+0.5d0/rbar)**3)*((1.d0+0.5d0/rbaryu)**3)
!!!                 psi6bar_zs = ((1.d0+0.5d0/rbar)**3)*((1.d0+0.5d0/rbarzu)**3)
!!!              else
!!!                 psi6bar_xs = (1.d0+r0/rbar)*(1.d0+r0/rbarxu) * & 
!!!                              sqrt( sqrt( 1.d0+2.d0/(rbar+r0) ) ) * & 
!!!                              sqrt( sqrt( 1.d0+2.d0/(rbarxu+r0) ) )
!!!                 psi6bar_ys = (1.d0+r0/rbar)*(1.d0+r0/rbaryu) * &
!!!                              sqrt( sqrt( 1.d0+2.d0/(rbar+r0) ) ) * &
!!!                              sqrt( sqrt( 1.d0+2.d0/(rbaryu+r0) ) )
!!!                 psi6bar_zs = (1.d0+r0/rbar)*(1.d0+r0/rbarzu) * &
!!!                              sqrt( sqrt( 1.d0+2.d0/(rbar+r0) ) ) * &
!!!                              sqrt( sqrt( 1.d0+2.d0/(rbarzu+r0) ) )
!!!              end if
!!!
!!!              ! Compute Bxtilde, Bytilde, Bx and By on stagger grid
!!!              call compute_Bxtilde_mag_bondi(x_stagger,yl,zl,amp,Bxts,r0,puncture_id)
!!!              call compute_Bytilde_mag_bondi(xl,y_stagger,zl,amp,Byts,r0,puncture_id)
!!!              Bx_stagger(i,j,k) = Bxts/psi6bar_xs
!!!              By_stagger(i,j,k) = Byts/psi6bar_ys
!!!
!!!              ! Compute Bztilde, Bz, Ax, Ay and Az on stagger grid
!!!              if (k==kmax) then 
!!!                 call compute_Bztilde_mag_bondi(xl,yl,zout,n_zout,dX,dY,dZ,amp,Bzts,r0,puncture_id)
!!!                 call compute_Ax_mag_bondi(xl,y_stagger,zout,n_zout,dZ,Ax(i,j,k),amp,r0,puncture_id)
!!!                 !! n_xout = int( (xout-x_stagger)/dX + 0.1d0)
!!!                   !! call compute_Ay_mag_bondi(xout,yl,zout,n_xout,n_zout, & 
!!!                !!                dX,dY,dZ,Ay(i,j,k),amp,r0,puncture_id)
!!!                n_x0 = int( (x_stagger-x0)/dX + 0.1d0)
!!!                call compute_Ay_mag_bondi(x0,yl,zout,n_x0,n_zout, &
!!!                             dX,dY,dZ,Ay(i,j,k),amp,r0,puncture_id)
!!!              else
!!!                 call compute_Bxtilde_mag_bondi(x_stagger,yl,zl+dZ,amp,Bxts1,r0,puncture_id)
!!!                 call compute_Bxtilde_mag_bondi(x_stagger-dX,yl,zl+dZ,amp,Bxts0,r0,puncture_id)
!!!                 call compute_Bytilde_mag_bondi(xl,y_stagger,zl+dZ,amp,Byts1,r0,puncture_id)
!!!                 call compute_Bytilde_mag_bondi(xl,y_stagger-dY,zl+dZ,amp,Byts0,r0,puncture_id)
!!!                 Bzts = Bzts1 + dZ*( (Bxts1-Bxts0)/dX + (Byts1-Byts0)/dY )
!!!                 Ax(i,j,k) = Ax(i,j,k+1) - dZ*Byts1
!!!                  Ay(i,j,k) = Ay(i,j,k+1) + dZ*Bxts1
!!!              end if
!!!              Bz_stagger(i,j,k) = Bzts/psi6bar_zs
!!!              Bzts1 = Bzts
!!!              Az(i,j,k) = 0.d0
!!!           end do
!!!        end do
!!!     end do
!!!
!!!     ! Now compute Bx, By and Bz on unstagger grid
!!!     call compute_Bi_mag_bondi(ext,Bx,By,Bz,Bx_stagger,By_stagger, & 
!!!              Bz_stagger,x,y,z,dx,dy,dz,puncture_id,r0,amp)
!!!
!!!     ! Compute conservative variables and other source terms
!!!     do k=kmin,kmax
!!!        do j=jmin,jmax
!!!           do i=imin,imax
!!!              ! Compute sb_i
!!!              ! Note that the metric in spherical coordinates is 
!!!              ! gamma_{ij} = (1+M/2rbar)^4 (drbar^2 + rbar^2 d\Omega^2)
!!!              ! in isotropic (puncture) radial coordinate. 
!!!              ! ds^2 = -(1-2M/r)dt^2 + 4M/r dr dt + (1+2M/r)dr^2 +r^2 d\Omega^2 
!!!              ! in Kerr-Schild coordinates.
!!!              !
!!!              ! We know that B^i is radial analytically, but if you use Ai to 
!!!              ! compute B^i, the resultant B^i will be radial only up to 
!!!              ! truncation error. Hence we compute the other components here.
!!!              !
!!!                 xl = x(i,j,k)
!!!              yl = y(i,j,k)
!!!                zl = z(i,j,k)
!!!              pm = sqrt(xl**2+yl**2)
!!!              rbar = sqrt(xl**2+yl**2+zl**2)
!!!              if (puncture_id==1) then
!!!                 r = rbar*(1.d0+0.5d0/rbar)**2
!!!                 psi6bar = (1.d0+0.5d0/rbar)**6
!!!                 ! Assume puncture lapse alph = 1/psi^2
!!!                 alph = 1.d0/(1.d0+0.5d0/rbar)**2
!!!                 u = Mdot / (FOUR * PI * r * r * rhob(i,j,k))
!!!                 u = u/(1.d0-0.25d0/rbar**2)
!!!                 ut = sqrt(1.d0+u*u*(1.d0+0.5d0/rbar)**4) / alph
!!!                 ! Limit velocity: alpha ut <= 60
!!!                 if (alph*ut .gt. 60.d0) then
!!!                    u = sqrt(3599.d0)/(1.d0+0.5d0/rbar)**2
!!!                   ut = 60.d0/alph
!!!                 end if
!!!                 vr = -u/ut
!!!                 u_r = ut*vr*(1.d0+0.5d0/rbar)**4
!!!                 u_t = -alph*alph*ut
!!!              else
!!!                 r = rbar + r0
!!!                 psi6bar = psi6 * (r/rbar)**2
!!!                 alph = 1.d0/sqrt(ONE + H)
!!!                 H = 2.d0*M/r
!!!                 u = Mdot / (FOUR * PI * r * r * rhob(i,j,k))
!!!                 u_r = -(ONE/(ONE-H))*(u - H * sqrt(ONE + u*u - H))
!!!                 u_t = -sqrt(ONE + u*u - H)
!!!                 vr = -u/ut
!!!              end if
!!!              u_x = u_r * xl/rbar
!!!              u_y = u_r * yl/rbar
!!!              u_z = u_r * zl/rbar
!!!              Bxtilde(i,j,k) = Bx(i,j,k) * psi6bar
!!!              Bytilde(i,j,k) = By(i,j,k) * psi6bar
!!!              Bztilde(i,j,k) = Bz(i,j,k) * psi6bar
!!!              Br = (Bx(i,j,k)*xl + By(i,j,k)*yl + Bz(i,j,k)*zl)/rbar
!!!              sbt(i,j,k) = u_r*Br/alph/sqrt(4.d0*PI)
!!!              sbr = Br/(alph*ut)/sqrt(4.d0*PI) + sbt(i,j,k)*vr
!!!              sbx(i,j,k) = Bx(i,j,k)/(alph*ut)/sqrt(4.d0*PI) +  & 
!!!                           sbt(i,j,k)*vx(i,j,k)
!!!              sby(i,j,k) = By(i,j,k)/(alph*ut)/sqrt(4.d0*PI) + & 
!!!                           sbt(i,j,k)*vy(i,j,k)
!!!              sbz(i,j,k) = Bz(i,j,k)/(alph*ut)/sqrt(4.d0*PI) + & 
!!!                           sbt(i,j,k)*vz(i,j,k)
!!!
!!!              ! sb_theta = gamma_{theta theta} sb^theta , 
!!!              ! sb_phi =  gamma_{phi phi} sb^phi in both coordinate systems.
!!!              sb_th = r*r/rbar/rbar*( xl*zl/pm*sbx(i,j,k) + & 
!!!                                yl*zl/pm*sby(i,j,k) - pm*sbz(i,j,k) ) 
!!!              sb_ph = r*r/rbar/rbar*(xl*sby(i,j,k)-yl*sbx(i,j,k))
!!!              if (puncture_id==1) then
!!!                 sb_r = sbr*(1.d0+0.5d0/rbar)**4
!!!                 sb_t = -sbt(i,j,k)*alph**2
!!!              else
!!!                 sb_r = sbr*(1.d0+2.d0/r) + 2.d0/r*sbt(i,j,k)
!!!                 sb_t = -(1.d0-2.d0/r)*sbt(i,j,k) + 2.d0/r*sbr
!!!              end if
!!!
!!!              sb_x = xl/rbar*sb_r + xl*zl/pm/rbar/rbar*sb_th & 
!!!                     - yl/pm/pm*sb_ph
!!!              sb_y = yl/rbar*sb_r + yl*zl/pm/rbar/rbar*sb_th &
!!!                     + xl/pm/pm*sb_ph
!!!              sb_z = zl/rbar*sb_r - pm/rbar/rbar*sb_th
!!!              sb2 = sbt(i,j,k)*sb_t + sbx(i,j,k)*sb_x + sby(i,j,k)*sb_y + sbz(i,j,k)*sb_z
!!!
!!!              mhd_st_x(i,j,k) = st_x(i,j,k) + alph*psi6bar*( sb2*ut*u_x - sbt(i,j,k)*sb_x ) 
!!!              mhd_st_y(i,j,k) = st_y(i,j,k) + alph*psi6bar*( sb2*ut*u_y - sbt(i,j,k)*sb_y )
!!!              mhd_st_y(i,j,k) = st_z(i,j,k) + alph*psi6bar*( sb2*ut*u_z - sbt(i,j,k)*sb_z )
!!!
!!!           end do
!!!        end do
!!!     end do
!!!
!!!  end if
end subroutine mag_bondi_accreting_matter
!-----------------------------------------------------------
! Find rho_0
!-----------------------------------------------------------
  REAL*8 FUNCTION find_rest_density(x1,x2,r,Gamma,Mdot,rho_inf,M,KPOLY)
    implicit none
    INTEGER MAXIT
    REAL*8 x1,x2,r,Gamma,Mdot,rho_inf,M,KPOLY
    REAL*8 xacc, scale
    logical success
!    EXTERNAL funcd
    PARAMETER (MAXIT=100)
    INTEGER j
    REAL*8 df,dx,dxold,f,fh,fl,temp,xh,xl
    xacc = 1.D-16
    success=.true.
    call transf(x1,fl,df,r,Gamma,Mdot,rho_inf,M,KPOLY)
    call transf(x2,fh,df,r,Gamma,Mdot,rho_inf,M,KPOLY)
    if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) then
       find_rest_density = x2 - 2.D0*fh/df
       call transf(find_rest_density,f,df,r,Gamma,Mdot,rho_inf,M,KPOLY)
       if(abs(f).lt.xacc)then
          return
       else
!          write(*,*) bracket: ,r,fl,fh
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
    call transf(find_rest_density,f,df,r,Gamma,Mdot,rho_inf,M,KPOLY)
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
       call transf(find_rest_density,f,df,r,Gamma,Mdot,rho_inf,M,KPOLY)
       if(f.lt.0.) then
          xl=find_rest_density
       else
          xh=find_rest_density
       endif
    end do
    success=.false.
    return
  END FUNCTION find_rest_density
!
  subroutine transf(rho0,f,df,r,Gamma,Mdot,rho_inf,M,KPOLY)
    implicit none
    real*8            :: rho0,f,df,r,Gamma,Mdot,M,KPOLY
    real*8            :: PI,ONE,TWO,u2,C_INF,rho_inf
    ONE = 1.D0
    TWO = 2.D0
    PI = acos(-ONE)
    C_INF = ( ONE + KPOLY*(Gamma/(Gamma-ONE))*(rho_inf**(Gamma-ONE)) )**2
    u2 = (Mdot/(4.D0*PI*rho0*r*r))**2
    f = ( ( ONE + KPOLY*(Gamma/(Gamma-ONE))*(rho0**(Gamma-ONE)) )**2 ) &
         * (ONE - TWO*M/r + u2) - C_INF
    df = TWO * (ONE + KPOLY*(Gamma/(Gamma-ONE)) * rho0**(Gamma-ONE)) &
         * KPOLY * Gamma * rho0**(Gamma-TWO) * (ONE - TWO*M/r + u2) &
         + u2 * (-TWO/rho0) * ((ONE + KPOLY*(Gamma/(Gamma-ONE)) &
         * rho0**(Gamma-ONE))**2)
!    print (1X,4e15.4), rho0,u2,f,TWO/r
  end subroutine transf
! Ax = -amp y / [ r(r+z) ]
subroutine compute_Ax_mag_bondi(x,y,z,amp,Ax)
  implicit none
  real*8 :: x,y,z,amp,Ax,r
!
  r = sqrt(x*x+y*y+z*z)
  Ax = -amp*y/( r * (r+z) )
end subroutine compute_Ax_mag_bondi
! Ay = amp x / [ r(r+z) ]
subroutine compute_Ay_mag_bondi(x,y,z,amp,Ay)
  implicit none
  real*8 :: x,y,z,amp,Ay,r
!
  r = sqrt(x*x+y*y+z*z)
  Ay = amp*x/( r * (r+z) )
end subroutine compute_Ay_mag_bondi
!
! Compute B^i
!
subroutine compute_Bi_mag_bondi(x,y,z,dx,dy,dz,amp,Bxs,Bys,Bzs,Bx,By,Bz, &
                                r0,puncture_id)
  implicit none
  integer :: puncture_id
  real*8 :: x,y,z,dx,dy,dz,amp,Bxs,Bys,Bzs,r0,rbar
  real*8 :: xu,yu,zu,rbar_xu,rbar_yu,rbar_zu
  real*8 :: psi6bar_xs,psi6bar_ys,psi6bar_zs,spsi6bar
  real*8 :: Bxt,Byt,Bzt,Bx,By,Bz,xd,yd,zd,Bxs0,Bys0,Bzs0
  real*8 :: rbar_xd,rbar_yd,rbar_zd,psi6bar_xs0,psi6bar_ys0,psi6bar_zs0
!
  xu = x + dx
  yu = y + dy
  zu = z + dz
  xd = x - dx
  yd = y - dy
  zd = z - dz
  rbar = sqrt(x*x + y*y + z*z)
  rbar_xu = sqrt(xu*xu + y*y + z*z)
  rbar_yu = sqrt(x*x + yu*yu + z*z)
  rbar_zu = sqrt(x*x + y*y + zu*zu)
  rbar_xd = sqrt(xd*xd + y*y + z*z)
  rbar_yd = sqrt(x*x + yd*yd + z*z)
  rbar_zd = sqrt(x*x + y*y + zd*zd)
  ! Compute psi^6 on staggered grid
  if (puncture_id==1) then
     spsi6bar = (1.d0+0.5d0/rbar)**3
     psi6bar_xs = spsi6bar*((1.d0+0.5d0/rbar_xu)**3)
     psi6bar_ys = spsi6bar*((1.d0+0.5d0/rbar_yu)**3)
     psi6bar_zs = spsi6bar*((1.d0+0.5d0/rbar_zu)**3)
     psi6bar_xs0 = spsi6bar*((1.d0+0.5d0/rbar_xd)**3)
     psi6bar_ys0 = spsi6bar*((1.d0+0.5d0/rbar_yd)**3)
     psi6bar_zs0 = spsi6bar*((1.d0+0.5d0/rbar_zd)**3)
  else
     spsi6bar = (1.d0+r0/rbar) * sqrt( sqrt( 1.d0+2.d0/(rbar+r0) ) )
     psi6bar_xs = spsi6bar*(1.d0+r0/rbar_xu) * sqrt( sqrt( 1.d0+2.d0/(rbar_xu+r0) ) )
     psi6bar_ys = spsi6bar*(1.d0+r0/rbar_yu) * sqrt( sqrt( 1.d0+2.d0/(rbar_yu+r0) ) )
     psi6bar_zs = spsi6bar*(1.d0+r0/rbar_zu) * sqrt( sqrt( 1.d0+2.d0/(rbar_zu+r0) ) )
     psi6bar_xs0 = spsi6bar*(1.d0+r0/rbar_xd) * sqrt( sqrt( 1.d0+2.d0/(rbar_xd+r0) ) )
     psi6bar_ys0 = spsi6bar*(1.d0+r0/rbar_yd) * sqrt( sqrt( 1.d0+2.d0/(rbar_yd+r0) ) )
     psi6bar_zs0 = spsi6bar*(1.d0+r0/rbar_zd) * sqrt( sqrt( 1.d0+2.d0/(rbar_zd+r0) ) )
  end if
  ! Compute Bitilde on staggered grid
  call compute_Bxtilde_mag_bondi(x,y,z,dx,dy,dz,amp,Bxt)
  call compute_Bytilde_mag_bondi(x,y,z,dx,dy,dz,amp,Byt)
  call compute_Bztilde_mag_bondi(x,y,z,dx,dy,dz,amp,Bzt)
  ! Compute B^i on staggered grid
  Bxs = Bxt / psi6bar_xs
  Bys = Byt / psi6bar_ys
  Bzs = Bzt / psi6bar_zs
  ! Compute B^i on unstaggered grid
  call compute_Bxtilde_mag_bondi(xd,y,z,dx,dy,dz,amp,Bxt)
  call compute_Bytilde_mag_bondi(x,yd,z,dx,dy,dz,amp,Byt)
  call compute_Bztilde_mag_bondi(x,y,zd,dx,dy,dz,amp,Bzt)
  Bxs0 = Bxt / psi6bar_xs0
  Bys0 = Byt / psi6bar_ys0
  Bzs0 = Bzt / psi6bar_zs0
  Bx = 0.5d0*(Bxs + Bxs0)
  By = 0.5d0*(Bys + Bys0)
  Bz = 0.5d0*(Bzs + Bzs0)
end subroutine compute_Bi_mag_bondi
!
! Compute Bxtilde at (x+dx/2, y, z) using Ai
!
subroutine compute_Bxtilde_mag_bondi(x,y,z,dx,dy,dz,amp,Bxtilde)
  implicit none
  real*8 :: x,y,z,dx,dy,dz,amp,Bxtilde,Ay1,Ay0
  real*8 :: x_stagger,z_stagger
!
  x_stagger = x+0.5d0*dx
  z_stagger = z+0.5d0*dz
  call compute_Ay_mag_bondi(x_stagger,y,z_stagger,amp,Ay1)
  call compute_Ay_mag_bondi(x_stagger,y,z_stagger-dz,amp,Ay0)
  Bxtilde = (Ay0-Ay1)/dz
end subroutine compute_Bxtilde_mag_bondi
!
! Compute Bytilde at (x, y+dy/2, z) using Ai
!
subroutine compute_Bytilde_mag_bondi(x,y,z,dx,dy,dz,amp,Bytilde)
  implicit none
  real*8 :: x,y,z,dx,dy,dz,amp,Bytilde,Ax1,Ax0
  real*8 :: y_stagger,z_stagger
!
  y_stagger = y+0.5d0*dy
  z_stagger = z+0.5d0*dz
  call compute_Ax_mag_bondi(x,y_stagger,z_stagger,amp,Ax1)
  call compute_Ax_mag_bondi(x,y_stagger,z_stagger-dz,amp,Ax0)
  Bytilde = (Ax1-Ax0)/dz
end subroutine compute_Bytilde_mag_bondi
!
! Compute Bztilde at (x, y, z+dz/2) using Ai
!
subroutine compute_Bztilde_mag_bondi(x,y,z,dx,dy,dz,amp,Bztilde)
  implicit none
  real*8 :: x,y,z,dx,dy,dz,amp,Bztilde,Ax1,Ax0,Ay1,Ay0
  real*8 :: y_stagger,x_stagger,z_stagger
!
  x_stagger = x+0.5d0*dx
  y_stagger = y+0.5d0*dy
  z_stagger = z+0.5d0*dz
  call compute_Ax_mag_bondi(x,y_stagger,z_stagger,amp,Ax1)
  call compute_Ax_mag_bondi(x,y_stagger-dy,z_stagger,amp,Ax0)
  call compute_Ay_mag_bondi(x_stagger,y,z_stagger,amp,Ay1)
  call compute_Ay_mag_bondi(x_stagger-dx,y,z_stagger,amp,Ay0)
  Bztilde = (Ay1-Ay0)/dx - (Ax1-Ax0)/dy
end subroutine compute_Bztilde_mag_bondi
!!!!! Compute Bx at (x+dx/2,y,z) using the vector potential
!!!!  subroutine compute_Bx_stagger_mag_bondi(x,y,z,dx,Bxs,puncture_id,r0,amp)
!!!!  implicit none
!!!!  integer :: puncture_id
!!!!  real*8 :: x,y,z,Bxs,r0,amp,dx
!!!!  real*8 :: xs,rs0,rs1
!!!!  real*8 :: Bxtildes,rbar0,rbar1,psi6bar
!!!!! 
!!!!  xs = x+0.5d0*dx
!!!!  call compute_Bxtilde_mag_bondi(xs,y,z,amp,Bxtildes,r0,puncture_id)
!!!!  
!!!!  rbar1 = sqrt((x+dx)**2 + y*y + z*z)
!!!!  rbar0 = sqrt(x*x + y*y + z*z)
!!!!  ! psi^6 at (x+dx/2,y,z) are given by exp[ 3 (phi0+phi1)]=psi0^3 * psi1^3, 
!!!!  ! where psi0 = psi(x,y,z), psi1 = psi(x+dx,y,z)
!!!!  if (puncture_id==1) then
!!!!     psi6bar = ((1.d0+0.5d0/rbar0)**3) * ((1.d0+0.5d0/rbar1)**3)
!!!!  else
!!!!     rs0 = rbar0 + r0
!!!!     rs1 = rbar1 + r0
!!!!     psi6bar = (sqrt(sqrt(1.d0+2.d0/rs0))*(rs0/rbar0)) * & 
!!!!               (sqrt(sqrt(1.d0+2.d0/rs1))*(rs1/rbar1))
!!!!  end if
!!!!
!!!!  Bxs = Bxtildes/psi6bar
!!!!
!!!!  end subroutine compute_Bx_stagger_mag_bondi
!!!!
!!!!! Compute By at (x,y+dy/2,z) using the vector potential
!!!!  subroutine compute_By_stagger_mag_bondi(x,y,z,dy,Bys,puncture_id,r0,amp)
!!!!  implicit none
!!!!  integer :: puncture_id
!!!!  real*8 :: x,y,z,Bys,r0,amp,dy
!!!!  real*8 :: ys,rs0,rs1
!!!!  real*8 :: Bytildes,rbar0,rbar1,psi6bar
!!!!!
!!!!  ys = y+0.5d0*dy
!!!!  call compute_Bytilde_mag_bondi(x,ys,z,amp,Bytildes,r0,puncture_id)
!!!!
!!!!  rbar1 = sqrt(x*x + (y+dy)**2 + z*z)
!!!!  rbar0 = sqrt(x*x + y*y + z*z)
!!!!  ! psi^6 at (x,y+dy/2,z) are given by exp[ 3 (phi0+phi1)]=psi0^3 * psi1^3,
!!!!  ! where psi0 = psi(x,y,z), psi1 = psi(x,y+dy,z)
!!!!  if (puncture_id==1) then
!!!!     psi6bar = ((1.d0+0.5d0/rbar0)**3) * ((1.d0+0.5d0/rbar1)**3)
!!!!  else
!!!!     rs0 = rbar0 + r0
!!!!     rs1 = rbar1 + r0
!!!!     psi6bar = (sqrt(sqrt(1.d0+2.d0/rs0))*(rs0/rbar0)) * &
!!!!               (sqrt(sqrt(1.d0+2.d0/rs1))*(rs1/rbar1))
!!!!  end if
!!!!
!!!!  Bys = Bytildes/psi6bar
!!!!
!!!!  end subroutine compute_By_stagger_mag_bondi
!!!!
!!!!! Set Bxtilde = amp x/r^3
!!!!  subroutine compute_Bxtilde_mag_bondi(x,y,z,amp,Bxtilde,r0,puncture_id)
!!!!  implicit none
!!!!  real*8 :: x,y,z,amp,Bxtilde,r0,rbar,r
!!!!  integer :: puncture_id
!!!!!
!!!!  rbar = sqrt(x*x+y*y+z*z)
!!!!  if (puncture_id==1) then 
!!!!     Bxtilde = amp * x/rbar**3
!!!!  else
!!!!     r = rbar+r0
!!!!     Bxtilde = amp * x/rbar/(r*r)
!!!!  end if
!!!!  end subroutine compute_Bxtilde_mag_bondi
!!!!
!!!!! Set Bytilde = amp y/r^3
!!!!  subroutine compute_Bytilde_mag_bondi(x,y,z,amp,Bytilde,r0,puncture_id)
!!!!  implicit none
!!!!  real*8 :: x,y,z,amp,Bytilde,r0,rbar,r
!!!!  integer :: puncture_id
!!!!!
!!!!  rbar = sqrt(x*x+y*y+z*z)
!!!!  if (puncture_id==1) then
!!!!     Bytilde = amp * y/rbar**3
!!!!  else
!!!!     r = rbar+r0
!!!!     Bytilde = amp * y/rbar/(r*r)
!!!!  end if
!!!!  end subroutine compute_Bytilde_mag_bondi
!!!!
!!!!! Compute Dx(Bxtilde) + Dy(Bytilde)
!!!!  subroutine compute_Dx_Bxt_plus_Dy_Byt(x,y,z,amp,dx,dy,Dx_Bxt_plus_Dy_Byt, & 
!!!!                                        r0,puncture_id)
!!!!  implicit none
!!!!  real*8 :: x,y,z,amp,dx,dy,Dx_Bxt_plus_Dy_Byt,Bxt1,Bxt0,Byt1,Byt0,r0
!!!!  integer :: puncture_id
!!!!!
!!!!  call compute_Bxtilde_mag_bondi(x+0.5d0*dx,y,z,amp,Bxt1,r0,puncture_id)
!!!!  call compute_Bxtilde_mag_bondi(x-0.5d0*dx,y,z,amp,Bxt0,r0,puncture_id)
!!!!  call compute_Bytilde_mag_bondi(x,y+0.5d0*dy,z,amp,Byt1,r0,puncture_id)
!!!!  call compute_Bytilde_mag_bondi(x,y-0.5d0*dy,z,amp,Byt0,r0,puncture_id)
!!!!  Dx_Bxt_plus_Dy_Byt = (Bxt1-Bxt0)/dx + (Byt1-Byt0)/dy
!!!!  end subroutine compute_Dx_Bxt_plus_Dy_Byt
!!!!
!!!!! Compute Bztilde at point (x,y,zout - n dz) using div(B)=0 and a freely 
!!!!! specifying function Bztilde(x,y,zout), where zout is a constant parameter. 
!!!!! 
!!!!  subroutine compute_Bztilde_mag_bondi(x,y,zout,n,dx,dy,dz,amp,Bztilde,r0,puncture_id)
!!!!  implicit none
!!!!  integer :: n,k,puncture_id
!!!!  real*8 :: x,y,z,zout,dx,dy,dz,amp,Bztilde,Bzt_zout
!!!!  real*8 :: Dx_Bxt_plus_Dy_Byt,r0,rbar
!!!!!
!!!!  ! Set Bztilde(x,y,zout) = zout/r^3
!!!!  rbar = sqrt(x*x+y*y+zout*zout)
!!!!  if (puncture_id==1) then 
!!!!     Bzt_zout = amp*zout/rbar**3
!!!!  else
!!!!     Bzt_zout = amp*zout/rbar/(rbar+r0)**2
!!!!  end if
!!!!  
!!!!  Bztilde = Bzt_zout
!!!!  z = zout - n*dz
!!!!  do k=1,n
!!!!     call compute_Dx_Bxt_plus_Dy_Byt(x,y,z+(k-0.5d0)*dz,amp,dx,dy,Dx_Bxt_plus_Dy_Byt,r0,puncture_id)
!!!!     Bztilde = Bztilde + dz*Dx_Bxt_plus_Dy_Byt
!!!!  end do
!!!!
!!!!  end subroutine compute_Bztilde_mag_bondi 
!!!!
!!!!! Specify Ax on z=zout plane
!!!!  subroutine compute_Ax_zout(x,y,zout,Ax,amp)
!!!!  implicit none
!!!!  real*8 :: x,y,zout,Ax,amp
!!!!  real*8 :: rbar
!!!!  real*8, parameter :: eps = 0.1d0
!!!!!
!!!!  rbar = sqrt(x*x + y*y + zout*zout)
!!!!  Ax = amp*y*zout/(x*x + y*y + eps)/rbar
!!!!  end subroutine compute_Ax_zout
!!!!
!!!!! Specify Ay on the line x=xout, z=zout
!!!!  subroutine compute_Ay_xout_zout(xout,y,zout,Ay,amp)
!!!!  implicit none
!!!!  real*8 :: xout,y,zout,Ay,amp
!!!!  real*8 :: rbar
!!!!  real*8, parameter :: eps = 0.1d0
!!!!!
!!!!  rbar = sqrt(xout*xout + y*y + zout*zout)
!!!!  Ay = -amp*xout*zout/(xout*xout + y*y + eps)/rbar
!!!!  end subroutine compute_Ay_xout_zout
!!!!
!!!!! Compute Ax at point (x,y,zout - n dz) by solving a finite difference equation
!!!!!
!!!!  subroutine compute_Ax_mag_bondi(x,y,zout,n,dz,Ax,amp,r0,puncture_id)
!!!!  implicit none
!!!!  real*8 :: x,y,zout,z,dz,Ax,amp,Bytilde,r0
!!!!  integer :: n,k,puncture_id
!!!!!
!!!!  call compute_Ax_zout(x,y,zout,Ax,amp)
!!!!  z = zout - n*dz
!!!!  do k=1,n
!!!!     call compute_Bytilde_mag_bondi(x,y,z+(k-0.5d0)*dz,amp,Bytilde,r0,puncture_id)
!!!!     Ax = Ax - dz*Bytilde
!!!!  end do
!!!!  end subroutine compute_Ax_mag_bondi
!!!!
!!!!! Compute Ay at point (xout - n dx, y, zout) by solving 
!!!!!  a finite difference equation
!!!!!!  subroutine compute_Ay_zout(xout,y,zout,n,dx,dy,Ay,amp,r0,puncture_id)
!!!!!!  implicit none
!!!!!!  real*8 :: xout,y,zout,dx,dy,Ay,amp,x,Bzt_zout,xi,rbar
!!!!!!  real*8 :: Ax1,Ax0,dy_Ax,r0
!!!!!!  integer :: n,i,puncture_id
!!!!!!!
!!!!!!  call compute_Ay_xout_zout(xout,y,zout,Ay,amp)
!!!!!!  x = xout - n*dx
!!!!!!  do i=1,n
!!!!!!     xi = x + (i-0.5d0)*dx
!!!!!!     rbar = sqrt(xi*xi+y*y+zout*zout)
!!!!!!     if (puncture_id==1) then
!!!!!!        Bzt_zout = amp*zout/rbar**3
!!!!!!     else
!!!!!!        Bzt_zout = amp*zout/rbar/(rbar+r0)**2
!!!!!!     end if
!!!!!!     call compute_Ax_zout(xi,y+0.5d0*dy,zout,Ax1,amp)
!!!!!!     call compute_Ax_zout(xi,y-0.5d0*dy,zout,Ax0,amp)
!!!!!!     dy_Ax = (Ax1-Ax0)/dy
!!!!!!     Ay = Ay - dx*(Bzt_zout + dy_Ax)
!!!!!!  end do
!!!!!!  end subroutine compute_Ay_zout
!!!!
!!!!! Compute Ay at point (xout - nx dx, y, zout - nz dz) by solving
!!!!!  a finite difference equation
!!!!!!  subroutine compute_Ay_mag_bondi(xout,y,zout,nx,nz,dx,dy,dz,Ay,amp,r0,puncture_id)
!!!!!!  implicit none
!!!!!!  real*8 :: xout,y,zout,dx,dy,dz,Ay,amp,x,z,Bxtilde,r0
!!!!!!  integer :: k,nx,nz,puncture_id
!!!!!!!
!!!!!!  call compute_Ay_zout(xout,y,zout,nx,dx,dy,Ay,amp,r0,puncture_id)
!!!!!!  x = xout - nx*dx
!!!!!!  z = zout - nz*dz
!!!!!!  do k=1,nz
!!!!!!     call compute_Bxtilde_mag_bondi(x,y,z+(k-0.5d0)*dz,amp,Bxtilde,r0,puncture_id)
!!!!!!     Ay = Ay + dz*Bxtilde
!!!!!!  end do
!!!!!!  end subroutine compute_Ay_mag_bondi
!!!!
!!!!! Compute Ay at point (x0 + n dx, y, zout) by solving
!!!!!  a finite difference equation
!!!!  subroutine compute_Ay_zout(x0,y,zout,n,dx,dy,Ay,amp,r0,puncture_id)
!!!!  implicit none
!!!!  real*8 :: x0,y,zout,dx,dy,Ay,amp,x,Bzt_zout,xi,rbar
!!!!  real*8 :: Ax1,Ax0,dy_Ax,r0
!!!!  integer :: n,i,puncture_id
!!!!!
!!!!  call compute_Ay_xout_zout(x0,y,zout,Ay,amp)
!!!!  x = x0 + n*dx
!!!!
!!!!  if (n==0) return
!!!!
!!!!  if (n .gt. 0) then 
!!!!     do i=1,n
!!!!        xi = x0 + (i-0.5d0)*dx
!!!!        rbar = sqrt(xi*xi+y*y+zout*zout)
!!!!        if (puncture_id==1) then
!!!!           Bzt_zout = amp*zout/rbar**3
!!!!        else
!!!!           Bzt_zout = amp*zout/rbar/(rbar+r0)**2
!!!!        end if
!!!!        call compute_Ax_zout(xi,y+0.5d0*dy,zout,Ax1,amp)
!!!!        call compute_Ax_zout(xi,y-0.5d0*dy,zout,Ax0,amp)
!!!!        dy_Ax = (Ax1-Ax0)/dy
!!!!        Ay = Ay + dx*(Bzt_zout + dy_Ax)
!!!!     end do
!!!!     return
!!!!  end if
!!!!
!!!!  if (n .lt. 0) then 
!!!!     do i=1,-n
!!!!         xi = x + (i-0.5d0)*dx
!!!!        rbar = sqrt(xi*xi+y*y+zout*zout)
!!!!        if (puncture_id==1) then
!!!!           Bzt_zout = amp*zout/rbar**3
!!!!        else
!!!!           Bzt_zout = amp*zout/rbar/(rbar+r0)**2
!!!!        end if
!!!!        call compute_Ax_zout(xi,y+0.5d0*dy,zout,Ax1,amp)
!!!!        call compute_Ax_zout(xi,y-0.5d0*dy,zout,Ax0,amp)
!!!!        dy_Ax = (Ax1-Ax0)/dy
!!!!        Ay = Ay - dx*(Bzt_zout + dy_Ax)
!!!!     end do
!!!!     return
!!!!  end if
!!!!
!!!!  end subroutine compute_Ay_zout
!!!!
!!!!! Compute Ay at point (x0 + nx dx, y, zout - nz dz) by solving
!!!!!  a finite difference equation
!!!!  subroutine compute_Ay_mag_bondi(x0,y,zout,nx,nz,dx,dy,dz,Ay,amp,r0,puncture_id)
!!!!  implicit none
!!!!  real*8 :: x0,y,zout,dx,dy,dz,Ay,amp,x,z,Bxtilde,r0
!!!!  integer :: k,nx,nz,puncture_id
!!!!!
!!!!  call compute_Ay_zout(x0,y,zout,nx,dx,dy,Ay,amp,r0,puncture_id)
!!!!  x = x0 + nx*dx
!!!!  z = zout - nz*dz
!!!!  do k=1,nz
!!!!     call compute_Bxtilde_mag_bondi(x,y,z+(k-0.5d0)*dz,amp,Bxtilde,r0,puncture_id)
!!!!     Ay = Ay + dz*Bxtilde
!!!!  end do
!!!!  end subroutine compute_Ay_mag_bondi
!!!!
!!!!  subroutine compute_Bz_kmin(x,y,z,Bz_kmin,Bzs_kmin,dx,dy,dz, & 
!!!!               puncture_id,r0,amp)
!!!!  implicit none
!!!!  integer :: puncture_id
!!!!  real*8 :: x,y,z,Bz_kmin,Bzs_kmin,dx,dy,dz,r0,amp,rbar
!!!!  real*8 :: z1,rbar1,psi6bar_zs,Bzts_kmin,Bzts_kminm1
!!!!  real*8 :: z0,rbar0,psi6bar_zs0,Bxts0,Byts0,Bzs_kminm1
!!!!  real*8 :: Bxts,Byts
!!!!!
!!!!  rbar = sqrt(x*x + y*y + z*z)
!!!!  z1 = z + dz
!!!!  z0 = z - dz
!!!!  rbar1 = sqrt(x*x + y*y + z1*z1)
!!!!  rbar0 = sqrt(x*x + y*y + z0*z0)
!!!!  if (puncture_id==1) then 
!!!!     psi6bar_zs = ((1.d0+0.5d0/rbar)**3)*((1.d0+0.5d0/rbar1)**3)
!!!!     psi6bar_zs0 = ((1.d0+0.5d0/rbar)**3)*((1.d0+0.5d0/rbar0)**3)
!!!!  else
!!!!     psi6bar_zs = (1.d0+r0/rbar)*(1.d0+r0/rbar1) * &
!!!!                              sqrt( sqrt( 1.d0+2.d0/(rbar+r0) ) ) * &
!!!!                              sqrt( sqrt( 1.d0+2.d0/(rbar1+r0) ) )
!!!!     psi6bar_zs0 = (1.d0+r0/rbar)*(1.d0+r0/rbar0) * &
!!!!                              sqrt( sqrt( 1.d0+2.d0/(rbar+r0) ) ) * &
!!!!                              sqrt( sqrt( 1.d0+2.d0/(rbar0+r0) ) )
!!!!  end if
!!!!  Bzts_kmin = Bzs_kmin * psi6bar_zs
!!!!
!!!!  call compute_Bxtilde_mag_bondi(x+0.5d0*dx,y,z,amp,Bxts,r0,puncture_id)
!!!!  call compute_Bytilde_mag_bondi(x,y+0.5d0*dy,z,amp,Byts,r0,puncture_id)
!!!!  call compute_Bxtilde_mag_bondi(x-0.5d0*dx,y,z,amp,Bxts0,r0,puncture_id)
!!!!  call compute_Bytilde_mag_bondi(x,y-0.5d0*dy,z,amp,Byts0,r0,puncture_id)
!!!!  Bzts_kminm1 = Bzts_kmin + dz*( (Bxts-Bxts0)/dx + (Byts-Byts0)/dy )
!!!!  Bzs_kminm1 = Bzts_kminm1/psi6bar_zs0
!!!!  Bz_kmin = 0.5d0*(Bzs_kmin + Bzs_kminm1)
!!!!  end subroutine compute_Bz_kmin
!!!!
!!!!  subroutine compute_Bi_mag_bondi(ext,Bx,By,Bz,Bxs,Bys,Bzs,x,y,z, & 
!!!!              dx,dy,dz,puncture_id,r0,amp)
!!!!  implicit none
!!!!  integer, dimension(3) :: ext
!!!!  real*8, dimension(ext(1),ext(2),ext(3)) :: Bx,By,Bz,Bxs,Bys,Bzs,x,y,z
!!!!  real*8 :: dx,dy,dz,r0,amp
!!!!  integer :: puncture_id,i,j,k
!!!!  real*8 :: Bxs0,Bys0,xl,yl,zl
!!!!!
!!!!  do k=1,ext(3)
!!!!     do j=1,ext(2)
!!!!        do i=1,ext(1)
!!!!
!!!!           if (i==1) then 
!!!!              xl = x(i,j,k) - dx
!!!!              yl = y(i,j,k)
!!!!              zl = z(i,j,k)
!!!!                   call compute_Bx_stagger_mag_bondi(xl,yl,zl,dx,Bxs0,puncture_id,r0,amp)
!!!!                   Bx(i,j,k) = 0.5d0*(Bxs0 + Bxs(i,j,k))
!!!!           else
!!!!                   Bx(i,j,k) = 0.5d0*(Bxs(i-1,j,k)+Bxs(i,j,k))
!!!!           end if
!!!!
!!!!           if (j==1) then 
!!!!              xl = x(i,j,k)
!!!!              yl = y(i,j,k) - dy
!!!!              zl = z(i,j,k)
!!!!                   call compute_By_stagger_mag_bondi(xl,yl,zl,dy,Bys0,puncture_id,r0,amp)
!!!!                   By(i,j,k) = 0.5d0*(Bys0+Bys(i,j,k))
!!!!           else
!!!!              By(i,j,k) = 0.5d0*(Bys(i,j-1,k)+Bys(i,j,k))
!!!!           end if
!!!!
!!!!           if (k==1) then 
!!!!              xl = x(i,j,k)
!!!!              yl = y(i,j,k)
!!!!              zl = z(i,j,k)
!!!!              call compute_Bz_kmin(xl,yl,zl,Bz(i,j,k),Bzs(i,j,k),dx,dy,dz, &
!!!!                   puncture_id,r0,amp)
!!!!           else
!!!!              Bz(i,j,k) = 0.5d0*(Bzs(i,j,k-1)+Bzs(i,j,k))
!!!!           end if
!!!!
!!!!        end do
!!!!     end do
!!!!  end do
!!!!
!!!!  end subroutine compute_Bi_mag_bondi
