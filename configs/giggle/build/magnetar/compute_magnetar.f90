!-----------------------------------------------------------------------------
!
! Set up matter field variables for hybrid EOS
!
!-----------------------------------------------------------------------------
subroutine compute_magnetar_hybrid(ex, DEPLETION, RADEQUAT, DELTA, &
     X, Y, Z, &
     neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
     lapse, phi, shiftx, shifty, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     rho, S, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     rho_star, tau, st_x, st_y, st_z, P, w, vx, vy, vz, &
     rho_b, u0, h, rho_b_atm, Omega_Frame,PhysR,eps_flag,K_poly,n)
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
           rho0 = max(rho_b(i,j,k),rho_b_atm)
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
           if (rho0>rho_b_atm) then
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
           Press = P_cold * DEPLETION
           eps = eps_cold * DEPLETION
           ! Do you want to keep the original epsilon? 
           if (eps_flag == 1) then
              eps = K_poly*n*rho0**f1on
              Press = P_cold + (gamma_th-1.d0)*rho0*(eps-eps_cold)
           end if
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
! Add a strong perturbation to the star, to emulate the induced perturbation in a BHNS run, 
!    to check for magnetic field growth:
!!$           vx(i,j,k) = vx(i,j,k) + 0.2D0*sin(theta)
!!$           vy(i,j,k) = vy(i,j,k) + 0.2D0*sin(theta)
!!$           vz(i,j,k) = vz(i,j,k) - 0.2D0*cos(theta)
! Use this for Fonts energy variable 
           tau(i,j,k)       = (au0m1 + U*hm1)*rho_star(i,j,k) - psi6*Press
           u0(i,j,k)        = ut
           rho_b(i,j,k)     = rho0
         end do
     end do
  end do
  !write(*,*) comp_mag_hyb: ,gxx(2,2,2)
end subroutine compute_magnetar_hybrid
!
! Introduce perturbations to seed the MRI
!
subroutine seed_mri(ext,X,Y,Z,rho_b,P,vx,vy,vz,Bx,By,Bz, &
                    rho_star,mhd_st_x,mhd_st_y,mhd_st_z, tau, &
                    alpham1,betax,betay,betaz,phi, &
                    gxx,gxy,gxz,gyy,gyz,gzz,rho_b_atm,n)
 implicit none
 integer, dimension(3)                                  :: ext
 real*8, dimension(ext(1),ext(2),ext(3))                :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))                :: rho_b,P,vx,vy,vz, cs
 real*8, dimension(ext(1),ext(2),ext(3))                :: sbt,sbx,sby,sbz
 real*8, dimension(ext(1),ext(2),ext(3))                :: Bx,By,Bz
 real*8, dimension(ext(1),ext(2),ext(3))                :: rho_star,tau
 real*8, dimension(ext(1),ext(2),ext(3))                :: mhd_st_x,mhd_st_y,mhd_st_z
 real*8, dimension(ext(1),ext(2),ext(3))                :: alpham1,phi, er,el,au0m1
 real*8, dimension(ext(1),ext(2),ext(3))                :: betax,betay, betaz
 real*8, dimension(ext(1),ext(2),ext(3))                :: gxx,gxy,gxz
 real*8, dimension(ext(1),ext(2),ext(3))                :: gyy,gyz,gzz
 real*8, dimension(ext(1),ext(2),ext(3))                :: u_x,u_y,u_z, u0
 real*8, dimension(ext(1),ext(2),ext(3))                :: sb_x,sb_y,sb_z,b2, h
 real*8                                                 :: k1, k2, k3, k4
 real*8                                                        :: n, PI, rho_b_atm, gamma
 integer                                                :: i,j,k,imin,jmin,kmin
 integer                                                :: imax,jmax,kmax
 integer                                                :: idum
 real*8                                                 :: temp
 real*8, parameter :: TINY = 1.d-14
 real*8, parameter :: perturb = 1.d-3
 real*8, parameter :: lam1 = 0.5
 real*8, parameter :: lam2 = 1.0
 real*8, parameter :: lam3 = 1.5
 real*8, parameter :: lam4 = 2.0
 real*8, parameter :: fac = 0.99
 real                                                   :: ran1
 PI = acos(-1.0)
 gamma = 1.d0 + 1.d0/n
!
 imin = lbound(P,1)
 imax = ubound(P,1)
 jmin = lbound(P,2)
 jmax = ubound(P,2)
 kmin = lbound(P,3)
 kmax = ubound(P,3)
 k1 = 2.d0*PI/lam1
 k2 = 2.d0*PI/lam2
 k3 = 2.d0*PI/lam3
 k4 = 2.d0*PI/lam4
 where (rho_b .gt. 1.1d0*rho_b_atm)
    cs = sqrt(gamma*P/(rho_b + (n+1.0)*P) )
 elsewhere
    cs = 0.d0
 end where
 idum = -10
 temp = ran1(idum)
 do k = kmin,kmax
    do j = jmin,jmax
       do i = imin, imax
          vx(i,j,k) = vx(i,j,k) + perturb*cs(i,j,k)*(cos(k1*Z(1,1,k)) + &
                      cos(k2*Z(1,1,k)) + cos(k3*Z(1,1,k)) + cos(k4*Z(1,1,k)))
          vy(i,j,k) = vy(i,j,k) + perturb*cs(i,j,k)*(cos(k1*Z(1,1,k)) + &
                      cos(k2*Z(1,1,k)) + cos(k3*Z(1,1,k)) + cos(k4*Z(1,1,k)))
          vz(i,j,k) = vz(i,j,k) + perturb*cs(i,j,k)*(cos(k1*Z(1,1,k)) + &
                      cos(k2*Z(1,1,k)) + cos(k3*Z(1,1,k)) + cos(k4*Z(1,1,k)))
!!$
!!$          temp = ran1(idum) 
!!$          temp = temp-0.5
!!$          if (i==imin) temp = 0.d0
!!$
!!$          vx(i,j,k) = vx(i,j,k) + temp*perturb*cs(i,j,k)
!!$
!!$          temp = ran1(idum)
!!$          temp = temp-0.5
!!$          vy(i,j,k) = vy(i,j,k) + temp*perturb*cs(i,j,k)
!!$
!!$          temp = ran1(idum)
!!$          temp = temp-0.5
!!$          if (k==kmin) temp = 0.d0
!!$          vz(i,j,k) = vz(i,j,k) + temp*perturb*cs(i,j,k)
       end do
    end do
 end do
 !now recalculate stuff.
! Compute al*u0-1
  er = exp(4.d0*phi)*(gxx*(vx + betax)**2 + &
       2.d0*gxy*(vx + betax)*(vy + betay) +         &
       2.d0*gxz*(vx + betax)*(vz + betaz) +         &
       gyy*(vy + betay)**2 +                           &
       2.d0*gyz*(vy + betay)*(vz + betaz) +         &
       gzz*(vz + betaz)**2 )/(alpham1+1.d0)**2
! *** Check for superluminal velocity ***
  where (er .gt. 1.d0)
    vx = (vx + betax)*sqrt(fac/er)-betax
    vy = (vy + betay)*sqrt(fac/er)-betay
    vz = (vz + betaz)*sqrt(fac/er)-betaz
    er = fac
  end where
! ***************************************
  el = sqrt(1.d0-er)
  au0m1 = er/el/(1.d0+el)
  u0 = (au0m1+1.d0)/(alpham1+1.d0)
  u_x = u0*exp(4.d0*phi)*(gxx*(betax + vx) + gxy*(betay+vy) + gxz*(betaz+vz))
  u_y = u0*exp(4.d0*phi)*(gxy*(betax + vx) + gyy*(betay+vy) + gyz*(betaz+vz))
  u_z = u0*exp(4.d0*phi)*(gxz*(betax + vx) + gyz*(betay+vy) + gzz*(betaz+vz))
! Calculate small b^{\mu}
  sbt = u_x*Bx + u_y*By + u_z*Bz
  sbx = Bx/u0 + vx*sbt
  sby = By/u0 + vy*sbt
  sbz = Bz/u0 + vz*sbt
  sb_x = exp(4.d0*phi)*(gxx*(betax*sbt+sbx) + gxy*(betay*sbt+sby) + gxz*(betaz*sbt+sbz))
  sb_y = exp(4.d0*phi)*(gxy*(betax*sbt+sbx) + gyy*(betay*sbt+sby) + gyz*(betaz*sbt+sbz))
  sb_z = exp(4.d0*phi)*(gxz*(betax*sbt+sbx) + gyz*(betay*sbt+sby) + gzz*(betaz*sbt+sbz))
  b2 = -(alpham1+1.d0)**2*(sbt)**2 + exp(4.d0*phi)*(gxx*(betax*sbt+sbx)**2 + &
       2.d0*gxy*(betax*sbt+sbx)*(betay*sbt+sby) + 2.d0*gxz*(betax*sbt+sbx)*(betaz*sbt+sbz) + &
       gyy*(betay*sbt+sby)**2 + 2.d0*gyz*(betay*sbt+sby)*(betaz*sbt+sbz) + &
       gzz*(betaz*sbt+sbz)**2)
  rho_star = (alpham1+1.d0)*u0*exp(6.d0*phi)*rho_b
  h = 1.d0+(n+1.d0)*P/rho_b
  mhd_st_x = rho_star*h*u_x + exp(6.d0*phi)*(b2*u0*u_x - sbt*sb_x)/(4.d0*PI)/(alpham1+1.d0)
  mhd_st_y = rho_star*h*u_y + exp(6.d0*phi)*(b2*u0*u_y - sbt*sb_y)/(4.d0*PI)/(alpham1+1.d0)
  mhd_st_z = rho_star*h*u_z + exp(6.d0*phi)*(b2*u0*u_z - sbt*sb_z)/(4.d0*PI)/(alpham1+1.d0)
  tau = (au0m1 + (alpham1+1.0)*u0*(n+1.d0)*P/rho_b)*rho_star + &
       exp(6.d0*phi)*b2*u0**2/(4.d0*PI) - exp(6.d0*phi)*(P+b2/(8.d0*PI)/(alpham1+1.d0)**2) - &
       exp(6.d0*phi)*sbt**2/(4.d0*PI)
end subroutine seed_mri
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
      NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
