!-----------------------------------------
! Set up magnetar in alt. atmosphere case
!-----------------------------------------
subroutine compute_magnetar_hybrid_alt_atmosphere(ex, DEPLETION, RADEQUAT, DELTA, &
     X, Y, Z, &
     neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
     lapse, phi, shiftx, shifty, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     rho, S, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     rho_star, tau, st_x, st_y, st_z, P, w, vx, vy, vz, &
     rho_b, u0, h, rho_b_atm_gf, pfloor_gf, Omega_Frame,PhysR,eps_flag,K_poly,n)
  implicit none
   interface
    subroutine compute_pcold_epscold(rhob, P_cold, eps_cold, &
                neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
     implicit none
     integer :: neos
     real*8  :: rhob, P_cold, eps_cold
     real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
     real*8, dimension(neos+1) :: k_tab, gamma_tab
    end subroutine compute_pcold_epscold
   end interface
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  integer                                   :: neos
  real*8, dimension(neos)                  :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1)                :: k_tab, gamma_tab
  real*8                                   :: rho_b_atm, gamma_th
  real*8                                   :: DEPLETION,Omega_Frame
  real*8                                   :: RADEQUAT, DELTA
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: lapse,phi,shiftx,shifty
  real*8, dimension(ex(1),ex(2),ex(3))     :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3))     :: rho,S,Sx,Sy,Sz,PhysR
  real*8, dimension(ex(1),ex(2),ex(3))     :: Sxx,Sxy,Sxz,Syy,Syz,Szz
  real*8, dimension(ex(1),ex(2),ex(3))     :: rho_star,tau,st_x,st_y,st_z
  real*8, dimension(ex(1),ex(2),ex(3))     :: P,w,vx,vy,vz,rho_b,u0,h
  real*8, dimension(ex(1),ex(2),ex(3))     :: rho_b_atm_gf, pfloor_gf
!
! Other variables:
!
  real*8                     :: rho0,Omega_R,Press,U,eps
  real*8                     :: u_x,u_y,u_z,alpha,mfac,psi6
  real*8                     :: pgxx,pgxy,pgxz,pgyy,pgyz,pgzz
  real*8                     :: fac2,w_scalar,p_met,betax,betay,betaz
  real*8                     :: gamma,v2,gamma2,v,ut,ux,uy
  integer                    :: i,j,k,l,m,MAT,eps_flag
  real*8                     :: R,cost,sint,R1,R2,phi1,phi2,u_r,t
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8                     :: theta,enth,alpi,gami,K_poly,n,f1on
  real*8                     :: r_Cook,mu_Cook,v_Cook
  real*8                     :: HALF, ONE, PI
  real*8                     :: x1,x2,y1,y2
  real*8                     :: riso,A,B,lower,upper
  real*8                     :: P_cold, eps_cold, hm1,au0m1
  parameter(HALF = 0.5D0, ONE = 1.D0)
  PI = acos(-ONE)
!
! compute
!
  imin = lbound(rho,1)
  jmin = lbound(rho,2)
  kmin = lbound(rho,3)
  imax = ubound(rho,1)
  jmax = ubound(rho,2)
  kmax = ubound(rho,3)
  f1on = 1.d0/n
  do i = imin, imax
     do j = jmin, jmax
        do k = kmin, kmax
           alpha = lapse(i,j,k) + ONE
           mfac  = exp(4.D0*phi(i,j,k))
           psi6  = exp(6.D0*phi(i,j,k))
           pgxx  = mfac*gxx(i,j,k)
           pgxy  = mfac*gxy(i,j,k)
           pgxz  = mfac*gxz(i,j,k)
           pgyy  = mfac*gyy(i,j,k)
           pgyz  = mfac*gyz(i,j,k)
           pgzz  = mfac*gzz(i,j,k)
           R    = sqrt(X(i,1,1)**2 + Y(1,j,1)**2 + Z(1,1,k)**2)
           riso = PhysR(i,j,k)
           cost = Z(1,1,k)/R
           sint = sqrt(1.D0 - cost**2)
           theta = acos(cost)
           p_met = P(i,j,k)
           rho0 = max(rho_b(i,j,k),rho_b_atm_gf(i,j,k))
!  Add the m=2 perturbation when checking for stability to the bar mode.
!           rho0 = rho0 * (1.D0 + DELTA * (X(i,1,1)**2 - Y(1,j,1)**2)/RADEQUAT**2)
!           if (rho0 .lt. 0.d0) rho0=0.d0
           Omega_R = vy(i,j,k)
           w_scalar = sqrt( (shiftx(i,j,k)**2 + shifty(i,j,k)**2) &
                / (X(i,1,1)**2 + Y(1,j,1)**2) )
           shiftx(i,j,k) = shiftx(i,j,k) - Omega_Frame*Y(1,j,1)
           shifty(i,j,k) = shifty(i,j,k) + Omega_Frame*X(i,1,1)
           betax = shiftx(i,j,k)
           betay = shifty(i,j,k)
           betaz = 0.d0
           v     = (Omega_R - w_scalar)*riso*sint*exp(-p_met)
           if (rho0>1.000001*rho_b_atm_gf(i,j,k)) then
              gamma = 1.D0/sqrt(1.D0 - v*v)
              au0m1 = v*v/(1.d0-v*v+sqrt(1.d0-v*v))
           else
              Omega_R = 0.D0
              v = 0.D0
              gamma = 1.D0
              au0m1 = 0.d0
           end if
           Omega_R = Omega_R - Omega_Frame
           ut    = gamma/alpha
           U     = gamma
           ux    = -ut*Omega_R*Y(1,j,1)
           uy    = ut*Omega_R*X(i,1,1)
           u_x   = pgxx*(ux + ut*betax) + pgxy*(uy + ut*betay)
           u_y   = pgxy*(ux + ut*betax) + pgyy*(uy + ut*betay)
           u_z   = pgxz*(ux + ut*betax) + pgyz*(uy + ut*betay)
           call compute_pcold_epscold(rho0, P_cold, eps_cold, &
                     neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
           Press = max(P_cold,pfloor_gf(i,j,k))
           eps = eps_cold + (Press-P_cold)/(gamma_th-1.d0)/rho0
!!$           ! This may not be right in the atmosphere, but it will be reset anyway.
!!$           Press = P_cold * DEPLETION
!!$           eps = eps_cold * DEPLETION
!!$           ! Do you want to keep the original epsilon? 
!!$           if (eps_flag == 1) then 
!!$              eps = K_poly*n*rho0**f1on
!!$              Press = P_cold + (gamma_th-1.d0)*rho0*(eps-eps_cold)
!!$           end if
           hm1 = eps + Press/rho0
           enth = 1.D0 + hm1
           fac2       = rho0*enth*U*U
           h(i,j,k)         = enth
           rho(i,j,k)       = fac2 - Press
           S(i,j,k)         = 3.D0*Press + enth*rho0*(U + 1.D0)*au0m1
           Sx(i,j,k)        = rho0*enth*U*u_x
           Sy(i,j,k)        = rho0*enth*U*u_y
           Sz(i,j,k)        = rho0*enth*U*u_z
           Sxx(i,j,k)       = Press*pgxx + rho0*enth*u_x*u_x
           Sxy(i,j,k)       = Press*pgxy + rho0*enth*u_x*u_y
           Sxz(i,j,k)       = Press*pgxz + rho0*enth*u_x*u_z
           Syy(i,j,k)       = Press*pgyy + rho0*enth*u_y*u_y
           Syz(i,j,k)       = Press*pgyz + rho0*enth*u_y*u_z
           Szz(i,j,k)       = Press*pgzz + rho0*enth*u_z*u_z
           rho_star(i,j,k)  = ut*alpha*psi6*rho0
           st_x(i,j,k)      = rho_star(i,j,k)*enth*u_x
           st_y(i,j,k)      = rho_star(i,j,k)*enth*u_y
           st_z(i,j,k)      = rho_star(i,j,k)*enth*u_z
           vx(i,j,k)        = ux/ut
           vy(i,j,k)        = uy/ut
           vz(i,j,k)        = 0.D0
           P(i,j,k)         = Press
           w(i,j,k)         = rho_star(i,j,k)*U
! Use this for Fonts energy variable 
           tau(i,j,k)       = (au0m1 + U*hm1)*rho_star(i,j,k) - psi6*Press
           u0(i,j,k)        = ut
           rho_b(i,j,k)     = rho0
         end do
     end do
  end do
end subroutine compute_magnetar_hybrid_alt_atmosphere