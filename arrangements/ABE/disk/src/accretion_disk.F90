!-----------------------------------------------------------------------------
!
! Set up matter source terms (disk accretion onto a black hole)
! This routine assumes axisymmetry
!
!-----------------------------------------------------------------------------
subroutine setup_fishbone_moncrief_disk_spin(ext, x, y, z, &
     rho_star, tau, st_x, st_y, st_z, w_rho, enth, rhob, P, &
     vx, vy, vz, u0, n, ell, W_in, rho_b_atm_gf, K_init, sam, Rp, dRdr, &
     pfloor_gf, lapse, phi, shiftx, shifty, shiftz, &
     gxx, gxy, gxz, gyy, gyz, gzz, enable_HARM_energyvariable, proc_imin, proc_kmin)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ext
  integer                                  :: nr
  real*8, dimension(ext(1),ext(2),ext(3))  :: rho_star,tau,rho_b_atm_gf, pfloor_gf
  real*8, dimension(ext(1),ext(2),ext(3))  :: st_x,st_y,st_z,w_rho,enth,rhob,P
  real*8, dimension(ext(1),ext(2),ext(3))  :: vx,vy,vz,u0, Rp, dRdr, x, y, z
  real*8, dimension(ext(1),ext(2),ext(3))  :: lapse, phi, shiftx, shifty, shiftz
  real*8, dimension(ext(1),ext(2),ext(3))  :: gxx, gxy, gxz, gyy, gyz, gzz
  integer                                  :: enable_HARM_energyvariable, proc_imin, proc_kmin
  real*8                                   :: n,ell,W_in,K_init, sam
  !
  ! Other variables:
  !
  integer                    :: i, j, k, kl, ku, glob_imin, glob_kmin
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8                     :: HALF, ONE, PI, ZERO, TWO, FOUR
  real*8                     :: rho0,epsilon,alpha,ut,u_x,u_y,u_z,log_h,h,M
  real*8                     :: varpi,sin2theta,psi6,Gam, r2, rp2, xp, yp, zp
  real*8                     :: cos2theta, costheta, fac_v
  real*8                     :: Delta, Delta_in, Sigma, Sigma_in
  real*8                     :: A, A_in, W, w2, hh, vyp
  real*8                     :: xi, xi_in, zeta, zeta_in, gphiphi, gphit
  real*8                     :: e4phi, shift_x, shift_y, shift_z, beta2, u_t
  logical                    :: fish_to_phys
  real, external             :: ran3
  parameter(HALF = 0.5D0, ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0, FOUR = 4.D0)
  PI = acos(-ONE)
  Gam = ONE + ONE/n 
  M = ONE
  !
  imin = lbound(rho_star,1)
  jmin = lbound(rho_star,2)
  kmin = lbound(rho_star,3)
  imax = ubound(rho_star,1)
  jmax = ubound(rho_star,2)
  kmax = ubound(rho_star,3)
  !
  glob_kmin = 0
  glob_imin = 0

  if (proc_imin == glob_imin) then
     imin = imin + 1
  end if

  if (proc_kmin == glob_kmin) then
     kmin = kmin + 1
  end if

  ! a little modification that assumes axisymmetry:
  jmin = jmin + 1
  jmax = jmax - 1

  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           r2 = x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2
           rp2 = Rp(i,j,k)*Rp(i,j,k)
           w2 = rp2 - sam * sam
           xp = x(i,j,k)*sqrt(rp2/r2)
           yp = y(i,j,k)*sqrt(rp2/r2)
           zp = z(i,j,k)*sqrt(rp2/r2)
           w2 = HALF *( w2 + sqrt( w2 * w2 + FOUR * sam*sam*zp*zp ))
           W= sqrt(w2)
           costheta = zp/W
           cos2theta = costheta*costheta
           sin2theta = ONE - cos2theta
           hh = M * w2 * W/( w2 * w2 + sam*sam*zp*zp )

           ! this is not the physical varpi, but it is only used
           ! for calculating v^y.
           varpi = sqrt(x(i,j,k)**2 + y(i,j,k)**2)

           if(W .gt. 0.9d0*W_in) then
              Delta    = W**2    - TWO*M*W   + sam*sam
              Delta_in = W_in**2 - TWO*M*W_in + sam*sam

              Sigma    = W**2    + sam*sam*cos2theta
              Sigma_in = W_in**2

              A    = (W*W + sam*sam)**2 - Delta*sam*sam*sin2theta
              A_in = W_in*(W_in**3 + W_in*sam*sam + TWO*M*sam*sam)

              zeta = ONE + sqrt(ONE + FOUR*ell**2*Sigma**2*Delta/(A**2*sin2theta))
              zeta = zeta * A / (Sigma * Delta)
              zeta_in = ONE + sqrt(ONE + FOUR*ell**2*Sigma_in**2*Delta_in/A_in**2)
              zeta_in = zeta_in * A_in / (Sigma_in * Delta_in)

              xi = -HALF * sqrt(ONE + FOUR*ell**2*Sigma**2*Delta/(A**2*sin2theta))
              xi = xi - TWO*sam*M*W*ell/A
              xi_in = -HALF * sqrt(ONE + FOUR*ell**2*Sigma_in**2*Delta_in/(A_in**2))
              xi_in = xi_in - TWO*sam*M*W_in*ell/A_in

              log_h = HALF*log(zeta) + xi - HALF*log(zeta_in) - xi_in 

              if(log_h.gt.ZERO) then
                 h = exp(log_h)
                 epsilon = (h - ONE)/Gam
                 rho0 = ( (Gam-ONE)*epsilon/(K_init) )**n
              else
                 h = ONE
                 epsilon = ZERO
                 rho0 = ZERO
              end if

              rhob(i,j,k) = rho0

              ! introduce a perturbation, then recalculate h
              epsilon = epsilon * (1.d0 + 4.e-2*(ran3(0)-0.5d0))
              P(i,j,k) = (Gam - ONE)*rho0*epsilon
              h = 1.d0 + Gam*epsilon
              rhob(i,j,k) = rho0
              if (P(i,j,k) .lt. pfloor_gf(i,j,k)) then
                 P(i,j,k) = pfloor_gf(i,j,k)
                 epsilon = P(i,j,k)/((Gam-ONE)*rho0)
                 h = 1.d0 + Gam*epsilon 
              end if
              enth(i,j,k) = h

              gphiphi = sin2theta*(Sigma + sam*sam*(1+TWO*M*W/Sigma)*sin2theta)
              gphit = -TWO*M*sam*W*sin2theta/Sigma
              ut = gphiphi+sqrt(gphiphi**2 + FOUR*ell**2*Delta*sin2theta)
              ut = sqrt(ut/(TWO*Delta*sin2theta)) 
              vx(i,j,k) = ZERO
              vy(i,j,k) = varpi/ut * (ell - gphit*ut**2)/(gphiphi*ut)
              vyp = vy(i,j,k)*sqrt(rp2/r2)
              vz(i,j,k) = ZERO

              !compute physical u_i
              u_x = TWO*hh*W*xp/(w2+sam*sam)*ut & 
                   *(ONE-sam*xp/(w2+sam*sam)*vyp)
              u_y = ut*(vyp-TWO*hh*sam*xp/(w2+sam*sam) &
                   *(ONE-sam*xp/(w2+sam*sam)*vyp))
              u_z = TWO*hh*zp/W*ut & 
                   *(ONE-sam*xp/(w2+sam*sam)*vyp)
              ! transform (locally) to phisheye.
              fish_to_phys = .true.
              call  trans_fish_phys_vj2(x(i,j,k), y(i,j,k), z(i,j,k), Rp(i,j,k), dRdr(i,j,k), &
                   u_x,u_y,u_z,fish_to_phys)
              alpha = 1.d0 + lapse(i,j,k)

              e4phi = exp(4.d0*phi(i,j,k))

              shift_x = e4phi*(gxx(i,j,k)*shiftx(i,j,k) + &
                   gxy(i,j,k)*shifty(i,j,k) + gxz(i,j,k)*shiftz(i,j,k))
              shift_y = e4phi*(gxy(i,j,k)*shiftx(i,j,k) + &
                   gyy(i,j,k)*shifty(i,j,k) + gyz(i,j,k)*shiftz(i,j,k))
              shift_z = e4phi*(gxz(i,j,k)*shiftx(i,j,k) + &
                   gyz(i,j,k)*shifty(i,j,k) + gzz(i,j,k)*shiftz(i,j,k))

              beta2 = shiftx(i,j,k)*shift_x + shifty(i,j,k)*shift_y + shiftz(i,j,k)*shift_z
              u_t = ut*(-alpha*alpha + beta2 + shift_y*vy(i,j,k))

              psi6 = exp(6.d0*phi(i,j,k))
              rho_star(i,j,k) = rho0*ut*alpha*psi6
              st_x(i,j,k)     = rho_star(i,j,k)*h*u_x
              st_y(i,j,k)     = rho_star(i,j,k)*h*u_y
              st_z(i,j,k)     = rho_star(i,j,k)*h*u_z
              w_rho(i,j,k)    = rho_star(i,j,k)*alpha*ut

              if (enable_HARM_energyvariable == 0) then
                 tau(i,j,k)      = w_rho(i,j,k)*h - psi6*P(i,j,k) - rho_star(i,j,k)
              else
                 tau(i,j,k) = -rho_star(i,j,k)*(1.d0+h*u_t)-alpha*psi6*P(i,j,k)
              end if

           else
              rho0 = ZERO
              rhob(i,j,k) = rho0
              vx(i,j,k) = ZERO
              vy(i,j,k) = ZERO
              vz(i,j,k) = ZERO
              psi6  = sqrt(ONE + TWO*hh)
              alpha = ONE/psi6
           end if

           if(rho0.lt.rho_b_atm_gf(i,j,k)) then
              rho0 = rho_b_atm_gf(i,j,k)
              rhob(i,j,k) = rho0
              alpha = 1.d0+lapse(i,j,k)
              psi6 = exp(6.d0*phi(i,j,k))
              ut = ONE/alpha
              u_t = -alpha
              rho_star(i,j,k) = rho0*ut*alpha*psi6

              ! compute the physical u_i
              u_x             = 0.d0
              u_y             = 0.d0
              u_z             = 0.d0
              vx(i,j,k) = -shiftx(i,j,k)
              vy(i,j,k) = -shifty(i,j,k)
              vz(i,j,k) = -shiftz(i,j,k)
              epsilon         = n*K_init*rho0**(Gam-1)
              h               = ONE + Gam*epsilon
              enth(i,j,k)     = h
              P(i,j,k)        = (Gam - ONE)*rho0*epsilon

              if (P(i,j,k) .lt. pfloor_gf(i,j,k)) then
                 P(i,j,k) = pfloor_gf(i,j,k)
                 epsilon = P(i,j,k)/(rho0*(Gam-1))
                 h = ONE + Gam*epsilon
              end if

              w_rho(i,j,k)    = rho_star(i,j,k)
              if (enable_HARM_energyvariable == 0) then
                 tau(i,j,k)   = Gam*psi6*rho0*epsilon - psi6*P(i,j,k)
              else
                 tau(i,j,k)   = -rho_star(i,j,k)*(1.d0+h*u_t) - alpha*psi6*P(i,j,k)
              end if

              st_x(i,j,k)     = rho_star(i,j,k)*h*u_x
              st_y(i,j,k)     = rho_star(i,j,k)*h*u_y
              st_z(i,j,k)     = rho_star(i,j,k)*h*u_z
           end if

           u0(i,j,k) = ut

        end do
     end do
  end do

end subroutine setup_fishbone_moncrief_disk_spin

!-----------------------------------------------------------------------------
!
! Set up matter source terms (disk accretion onto a black hole)
! This routine assumes axisymmetry
! this version is when you use "w" for the code coordinate radius
!
!-----------------------------------------------------------------------------
subroutine setup_fishbone_moncrief_disk_spin_bl(ext, x, y, z, &
     rho_star, tau, st_x, st_y, st_z, w_rho, enth, rhob, P, &
     vx, vy, vz, u0, n, ell, W_in, rho_b_atm_gf, K_init, sam, & 
     Rp, dRdr, pfloor_gf, lapse, phi, shiftx, shifty, shiftz, &
     gxx, gxy, gxz, gyy, gyz, gzz,enable_HARM_energyvariable,proc_imin,proc_kmin)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ext
  integer                                  :: nr
  real*8, dimension(ext(1),ext(2),ext(3))  :: rho_star,tau,rho_b_atm_gf, pfloor_gf
  real*8, dimension(ext(1),ext(2),ext(3))  :: st_x,st_y,st_z,w_rho,enth,rhob,P
  real*8, dimension(ext(1),ext(2),ext(3))  :: vx,vy,vz,u0, Rp, dRdr, x, y, z
  real*8, dimension(ext(1),ext(2),ext(3))  :: lapse, phi, shiftx, shifty, shiftz
  real*8, dimension(ext(1),ext(2),ext(3))  :: gxx, gxy, gxz, gyy, gyz, gzz
  real*8                                   :: n,ell,W_in,K_init, sam
  integer                                  :: enable_HARM_energyvariable
  integer                                  :: proc_imin, proc_kmin
  !
  ! Other variables:
  !
  integer                    :: i, j, k, kl, ku
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: glob_kmin, glob_imin
  real*8                     :: HALF, ONE, PI, ZERO, TWO, FOUR
  real*8                     :: rho0,epsilon,alpha,ut,u_x,u_y,u_z,log_h,h,M
  real*8                     :: varpi,sin2theta,psi6,Gam, r2, rp2, xp, yp, zp
  real*8                     :: cos2theta, costheta, fac_v
  real*8                     :: Delta, Delta_in, Sigma, Sigma_in
  real*8                     :: A, A_in, W, w2, hh, vyp
  real*8                     :: xi, xi_in, zeta, zeta_in, gphiphi, gphit
  real*8                     :: shift_x, shift_y, shift_z, e4phi
  real*8                     :: beta2, u_t
  logical                    :: fish_to_phys
  parameter(HALF = 0.5D0, ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0, FOUR = 4.D0)
  PI = acos(-ONE)
  Gam = ONE + ONE/n 
  M = ONE
  !
  imin = lbound(rho_star,1)
  jmin = lbound(rho_star,2)
  kmin = lbound(rho_star,3)
  imax = ubound(rho_star,1)
  jmax = ubound(rho_star,2)
  kmax = ubound(rho_star,3)

  glob_imin = 0
  glob_kmin = 0

  !
  if (proc_imin == glob_imin) then
     imin = imin + 1
  end if

  if (proc_kmin == glob_kmin) then
     kmin = kmin + 1
  end if

  ! a little modification that assumes axisymmetry:
  jmin = jmin + 1
  jmax = jmax - 1

  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           r2 = x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2
           rp2 = Rp(i,j,k)*Rp(i,j,k)
           w2 = rp2 
           xp = x(i,j,k)*sqrt(rp2/r2)
           yp = y(i,j,k)*sqrt(rp2/r2)
           zp = z(i,j,k)*sqrt(rp2/r2)
           W = sqrt(w2)
           costheta = zp/W
           cos2theta = costheta*costheta
           sin2theta = ONE - cos2theta

           ! this is not the physical varpi, but it is only used
           ! for calculating v^y.
           varpi = sqrt(x(i,j,k)**2 + y(i,j,k)**2)

           if(W .gt. 0.9d0*W_in) then
              Delta    = W**2    - TWO*M*W    + sam*sam
              Delta_in = W_in**2 - TWO*M*W_in + sam*sam

              Sigma    = W**2    + sam*sam*cos2theta
              Sigma_in = W_in**2

              A    = (W*W + sam*sam)**2 - Delta*sam*sam*sin2theta
              A_in = W_in*(W_in**3 + W_in*sam*sam + TWO*M*sam*sam)

              zeta = ONE + sqrt(ONE + FOUR*ell**2*Sigma**2*Delta/(A**2*sin2theta))
              zeta = zeta * A / (Sigma * Delta)
              zeta_in = ONE + sqrt(ONE + FOUR*ell**2*Sigma_in**2*Delta_in/A_in**2)
              zeta_in = zeta_in * A_in / (Sigma_in * Delta_in)

              xi = -HALF * sqrt(ONE + FOUR*ell**2*Sigma**2*Delta/(A**2*sin2theta))
              xi = xi - TWO*sam*M*W*ell/A
              xi_in = -HALF * sqrt(ONE + FOUR*ell**2*Sigma_in**2*Delta_in/(A_in**2))
              xi_in = xi_in - TWO*sam*M*W_in*ell/A_in

              log_h = HALF*log(zeta) + xi - HALF*log(zeta_in) - xi_in 

              if(log_h.gt.ZERO) then
                 h = exp(log_h)
                 epsilon = (h - ONE)/Gam
                 rho0 = ( (Gam-ONE)*epsilon/(K_init) )**n
              else
                 h = ONE
                 epsilon = ZERO
                 rho0 = ZERO
              end if

              enth(i,j,k) = h
              rhob(i,j,k) = rho0
              P(i,j,k) = (Gam - ONE)*rho0*epsilon
              if (P(i,j,k) .lt. pfloor_gf(i,j,k)) then
                 P(i,j,k) = pfloor_gf(i,j,k)
                 epsilon = P(i,j,k)/((Gam-ONE)*rho0)
                 h = 1.d0 + Gam*epsilon 
              end if

              gphiphi = sin2theta*(Sigma + sam*sam*(1+TWO*M*W/Sigma)*sin2theta)
              gphit = -TWO*M*sam*W*sin2theta/Sigma
              ut = gphiphi+sqrt(gphiphi**2 + FOUR*ell**2*Delta*sin2theta)
              ut = sqrt(ut/(TWO*Delta*sin2theta)) 
              vx(i,j,k) = ZERO
              vy(i,j,k) = varpi/ut * (ell - gphit*ut**2)/(gphiphi*ut)
              vyp = vy(i,j,k)*sqrt(rp2/r2)
              vz(i,j,k) = ZERO

              !compute fisheye coordinate values of  u_i
              e4phi = exp(4.d0*phi(i,j,k))

              shift_x = e4phi*(gxx(i,j,k)*shiftx(i,j,k) + &
                   gxy(i,j,k)*shifty(i,j,k) + gxz(i,j,k)*shiftz(i,j,k))
              shift_y = e4phi*(gxy(i,j,k)*shiftx(i,j,k) + &
                   gyy(i,j,k)*shifty(i,j,k) + gyz(i,j,k)*shiftz(i,j,k))
              shift_z = e4phi*(gxz(i,j,k)*shiftx(i,j,k) + &
                   gyz(i,j,k)*shifty(i,j,k) + gzz(i,j,k)*shiftz(i,j,k))

              u_x = ut*(shift_x + e4phi*gxy(i,j,k)*vy(i,j,k))
              u_y = ut*(shift_y + e4phi*gyy(i,j,k)*vy(i,j,k))
              u_z = ut*(shift_z + e4phi*gyz(i,j,k)*vy(i,j,k))

              alpha = 1.d0 + lapse(i,j,k)
              beta2 = shiftx(i,j,k)*shift_x + shifty(i,j,k)*shift_y + shiftz(i,j,k)*shift_z
              u_t = ut*(-alpha*alpha + beta2 + shift_y*vy(i,j,k))
              psi6 = exp(6.d0*phi(i,j,k))
              rho_star(i,j,k) = rho0*ut*alpha*psi6
              st_x(i,j,k)     = rho_star(i,j,k)*h*u_x
              st_y(i,j,k)     = rho_star(i,j,k)*h*u_y
              st_z(i,j,k)     = rho_star(i,j,k)*h*u_z
              w_rho(i,j,k)    = rho_star(i,j,k)*alpha*ut
              if (enable_HARM_energyvariable == 0) then
                 tau(i,j,k) = w_rho(i,j,k)*h - psi6*P(i,j,k) - rho_star(i,j,k)
                 if(i==2 .and. j==2 .and. k==2) write(*,*) "HIIITAU1:",tau(i,j,k)
              else
                 tau(i,j,k) = -rho_star(i,j,k)*(1.d0+h*u_t)-alpha*psi6*P(i,j,k)
                 if(i==2 .and. j==2 .and. k==2) write(*,*) "HIIITAU1BAD:",tau(i,j,k)
              end if

           else
              rho0 = ZERO
              rhob(i,j,k) = rho0
              vx(i,j,k) = ZERO
              vy(i,j,k) = ZERO
              vz(i,j,k) = ZERO
              psi6  = exp(6.d0*phi(i,j,k))
              alpha = 1.d0 + lapse(i,j,k)
           end if

           ! Is the correct metric being used in the low density regions?

           if(rho0.lt.rho_b_atm_gf(i,j,k)) then
              rho0 = rho_b_atm_gf(i,j,k)
              rhob(i,j,k) = rho0
              alpha = 1.d0+lapse(i,j,k)
              psi6 = exp(6.d0*phi(i,j,k))
              ut = ONE/alpha
              u_t = -alpha
              rho_star(i,j,k) = rho0*ut*alpha*psi6

              ! compute the fisheye coordinate values of u_i
              u_x             = 0.d0
              u_y             = 0.d0
              u_z             = 0.d0
              vx(i,j,k) = -shiftx(i,j,k)
              vy(i,j,k) = -shifty(i,j,k)
              vz(i,j,k) = -shiftz(i,j,k)

              epsilon         = n*K_init*rho0**(Gam-1)
              h               = ONE + Gam*epsilon
              enth(i,j,k)     = h
              P(i,j,k)        = (Gam - ONE)*rho0*epsilon

              if (P(i,j,k) .lt. pfloor_gf(i,j,k)) then
                 P(i,j,k) = pfloor_gf(i,j,k)
                 epsilon = P(i,j,k)/(rho0*(Gam-1))
                 h = ONE + Gam*epsilon
              end if

              w_rho(i,j,k)    = rho_star(i,j,k)
              if (enable_HARM_energyvariable == 0) then
                 tau(i,j,k)   = Gam*psi6*rho0*epsilon - psi6*P(i,j,k)
                 if(i==2 .and. j==2 .and. k==2) write(*,*) "HIIITAU2:",tau(i,j,k)
              else
                 tau(i,j,k)   = -rho_star(i,j,k)*(1.d0+h*u_t) - alpha*psi6*P(i,j,k)
                 if(i==2 .and. j==2 .and. k==2) write(*,*) "HIIITAU2BAD:",tau(i,j,k)
              end if
              st_x(i,j,k)     = rho_star(i,j,k)*h*u_x
              st_y(i,j,k)     = rho_star(i,j,k)*h*u_y
              st_z(i,j,k)     = rho_star(i,j,k)*h*u_z
           end if

           u0(i,j,k) = ut

        end do
     end do
  end do

end subroutine setup_fishbone_moncrief_disk_spin_bl

! I wouldn't call this if I were you.  It hasn't really been
! maintained.
subroutine setup_hsw_disk(ext, x, y, z, &
     rho_star, tau, st_x, st_y, st_z, w, enth, rhob, P, &
     vx, vy, vz, u0, n, ell_in, r_in, rho_tiny, K_init,p_deplete, &
     proc_imin, proc_kmin)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ext
  integer                                  :: nr, proc_imin, proc_kmin
  real*8, dimension(ext(1),ext(2),ext(3))  :: rho_star,tau, x, y, z
  real*8, dimension(ext(1),ext(2),ext(3))  :: st_x,st_y,st_z,w,enth,rhob,P
  real*8, dimension(ext(1),ext(2),ext(3))  :: vx,vy,vz,u0
  real*8                                   :: n,ell_in,r_in,rho_tiny,K_init
  !
  ! Other variables:
  !
  integer                    :: i, j, k, kl, ku, glob_imin, glob_kmin
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8                     :: HALF, ONE, PI, ZERO, TWO, FOUR
  real*8                     :: rho0,epsilon,alpha,ut,u_x,u_y,u_z,log_h,h,M
  real*8                     :: r,varpi,sin2theta,A,A_in,bind,bind_in,psi6,Gam
  real*8                     :: ell,Omega,Omega_in,q,lambda,lambda_in,f,f_in
  real*8                     :: p_deplete
  parameter(HALF = 0.5D0, ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0, FOUR = 4.D0)
  PI = acos(-ONE)
  Gam = ONE + ONE/n 
  M = ONE
  q = 1.68d0
  !
  imin = lbound(rho_star,1)
  jmin = lbound(rho_star,2)
  kmin = lbound(rho_star,3)
  imax = ubound(rho_star,1)
  jmax = ubound(rho_star,2)
  kmax = ubound(rho_star,3)

  glob_imin = 0
  glob_kmin = 0

  !
  if (proc_imin == glob_imin) then
     imin = imin + 1
  end if

  if (proc_kmin == glob_kmin) then
     kmin = kmin + 1
  end if

  ! a little modification that assumes axisymmetry:
  jmin = jmin + 1
  jmax = jmax - 1

  !
  !  write(*,*) r_in,ell_in,K_init

  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           r = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
           if(r.gt.0.9d0*r_in) then
              varpi = sqrt(x(i,j,k)**2 + y(i,j,k)**2)
              sin2theta = (varpi/r)**2
              A = ONE - TWO*M/r
              A_in = ONE - TWO*M/r_in
              lambda = varpi/sqrt(A)
              lambda_in = r_in/sqrt(A_in)
              Omega_in = ell_in/lambda_in**2
              Omega = Omega_in * ( sqrt(A/A_in)/(varpi/r_in) )**q
              ell = Omega * lambda**2
              bind = ONE/sqrt(ONE/A - ell**2/(varpi**2))
              bind_in = ONE/sqrt(ONE/A_in - ell_in**2/(r_in**2))
              alpha = q/(q-TWO)
              f = abs(ONE - ell*Omega)**(ONE/(alpha + ONE))
              f_in = abs(ONE - ell_in*Omega_in)**(ONE/(alpha + ONE))

              h = bind_in * f_in / ( bind * f )

              if(h.gt.ONE) then
                 epsilon = (h - ONE)/Gam
                 rho0 = ( (Gam-ONE)*epsilon/(K_init) )**n
                 ! Now add pressure perturbation 
		 epsilon = epsilon*p_deplete
	         h = ONE + Gam*epsilon
              else
                 h = ONE
                 epsilon = ZERO
                 rho0 = ZERO
              end if

              enth(i,j,k) = h
              rhob(i,j,k) = rho0
              P(i,j,k) = (Gam - ONE)*rho0*epsilon
              vx(i,j,k) = ZERO
              vy(i,j,k) = varpi * Omega
              vz(i,j,k) = ZERO
              ut = ONE/sqrt(ONE - TWO*M/r - vy(i,j,k)**2)
              psi6  = sqrt(ONE + TWO*M/r)
              alpha = ONE/psi6
              u_x = TWO*M*ut*x(i,j,k)/(r*r)
              u_y = ut*vy(i,j,k)
              u_z = TWO*M*ut*z(i,j,k)/(r*r)

              rho_star(i,j,k) = rho0*ut
              st_x(i,j,k)     = rho_star(i,j,k)*h*u_x
              st_y(i,j,k)     = rho_star(i,j,k)*h*u_y
              st_z(i,j,k)     = rho_star(i,j,k)*h*u_z
              w(i,j,k)        = rho_star(i,j,k)*alpha*ut
              tau(i,j,k)      = w(i,j,k)*h - psi6*P(i,j,k) - rho_star(i,j,k)

           else
              rho0 = ZERO
              rhob(i,j,k) = rho0
              vx(i,j,k) = ZERO
              vy(i,j,k) = ZERO
              vz(i,j,k) = ZERO
              psi6  = sqrt(ONE + TWO*M/r)
              alpha = ONE/psi6
           end if

           if(rho0.lt.rho_tiny) then
              rho0 = rho_tiny 
              rhob(i,j,k) = rho0
              psi6  = sqrt(ONE + TWO*M/r)
              alpha = ONE/psi6
              ut = ONE/sqrt(max(ONE - TWO*M/r,1.d-10))
              u_x = TWO*M*ut*x(i,j,k)/(r*r)
              u_y = ZERO
              u_z = TWO*M*ut*z(i,j,k)/(r*r)
              rho_star(i,j,k) = rho0*ut
              epsilon         = n*K_init*rho0**(Gam-1)
              h               = ONE + Gam*epsilon
              enth(i,j,k)     = h
              P(i,j,k)        = (Gam - ONE)*rho0*epsilon
              w(i,j,k)        = alpha*rho_star(i,j,k)*ut
              tau(i,j,k)      = (alpha*ut*h - 1.d0)*rho_star(i,j,k) - &
                   psi6*P(i,j,k)
              st_x(i,j,k)     = rho_star(i,j,k)*h*u_x
              st_y(i,j,k)     = rho_star(i,j,k)*h*u_y
              st_z(i,j,k)     = rho_star(i,j,k)*h*u_z

           end if

           u0(i,j,k) = ut

        end do
     end do
  end do

end subroutine setup_hsw_disk

!------------------------------------------------------------
! Setup the initial EM fields from a vector potential 
!------------------------------------------------------------
!
subroutine initial_emfields_disk(ext,r_in,rhob_max,rho_cut,betam1,X,Y,Z,phi, &
     alpham1,shiftx,shifty,shiftz, & 
     gxx,gxy,gxz,gyy,gyz,gzz,rho_b, P,Bx,By,Bz, sb2, &
     st_x,st_y,st_z,mhd_st_x,mhd_st_y,mhd_st_z,tau,u0, &
     vx,vy,vz,Symmetry, Sym_Bz, fac, enable_HARM_energyvariable,Rp,dRdr, &
     A_phi,A_phix,A_phiy,A_phiz, &
     Bx_phys,By_phys,Bz_phys)
  implicit none
  integer, dimension(3)                   :: ext
  real*8, dimension(ext(1),ext(2),ext(3)) :: X, Y, Z
  real*8, dimension(ext(1),ext(2),ext(3)) :: Rp,dRdr
  real*8, dimension(ext(1),ext(2),ext(3)) :: rho_b,P,Bx,By,Bz,sb2
  real*8, dimension(ext(1),ext(2),ext(3)) :: mhd_st_x,mhd_st_y,mhd_st_z,tau
  real*8, dimension(ext(1),ext(2),ext(3)) :: st_x,st_y,st_z
  real*8, dimension(ext(1),ext(2),ext(3)) :: u0,vx,vy,vz
  real*8, dimension(ext(1),ext(2),ext(3)) :: u_x,u_y,u_z
  real*8                                  :: B_xs, B_ys, B_zs, B2s
  real*8, dimension(ext(1),ext(2),ext(3)) :: phi,psi_n,alpham1
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: shiftx,shifty,shiftz
  real*8                                  :: r_in,rho_cut,betam1,rhob_max
  real*8                                  :: dX,dY,dZ,al,sqrtg4,sqrtg
  real*8                                  :: sb0,sb2_l,sb_x,sb_y,sb_z,fac
  real*8                                  :: u0u0, b0b0, sbx, sby, sbz
  integer                                 :: Symmetry,flag,enable_HARM_energyvariable
  real*8, dimension(ext(1),ext(2),ext(3)) :: A_phi,A_phix,A_phiy,A_phiz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Bx_phys,By_phys,Bz_phys
  !
  integer                                        :: AXISYM,EQUATORIAL
  integer                                        :: OCTANT
  integer                                        :: i,j,k,imin,imax
  integer                                        :: jmin,jmax,kmin,kmax
  real*8, parameter 				:: SYM = 1.d0, ANTI = -1.d0
  real*8						:: pomega2,psim6, fs4pi,psin
  real*8                                         :: Sym_Bz
  real*8                                         :: lambda, omega
  real*8                                         :: PI
  parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
  !
  PI = acos(-1.d0)
  fs4pi = sqrt(4.d0*acos(-1.d0))
  imin = lbound(P,1)
  imax = ubound(P,1)
  jmin = lbound(P,2)
  jmax = ubound(P,2)
  kmin = lbound(P,3)
  kmax = ubound(P,3)

  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)

  psi_n = exp(4.d0*phi)*u0
  u_x = (gxx*(shiftx+vx) + gxy*(shifty+vy) + gxz*(shiftz+vz))*psi_n
  u_y = (gxy*(shiftx+vx) + gyy*(shifty+vy) + gyz*(shiftz+vz))*psi_n
  u_z = (gxz*(shiftx+vx) + gyz*(shifty+vy) + gzz*(shiftz+vz))*psi_n
  !
  ! Compute the vector potential A_phi 
  !
  call compute_Aphi_disk(ext,Z,r_in,rho_b,A_phi,rho_cut,rhob_max,Sym_Bz)
  !  call compute_Aphi_omega(ext,X,Z,rho_b,vy,A_phi,r_in,rho_cut,rhob_max,Sym_Bz)
  ! 
  ! Compute the derivatives of A_phi
  !
  if (Symmetry==OCTANT) then
     call gderivs_oct(ext,A_phi,A_phix,A_phiy,A_phiz,dX,dY,dZ,SYM,SYM,Sym_Bz)
  elseif (Symmetry==EQUATORIAL) then
     call gderivs_eq(ext,A_phi,A_phix,A_phiy,A_phiz,dX,dY,dZ,SYM,SYM ,Sym_Bz )
  elseif (Symmetry==AXISYM) then
     call gderivs_axi(ext,A_phi,A_phix,A_phiy,A_phiz,dX,dY,dZ,SYM,SYM,Sym_Bz)
  else
     call gderivs(ext,A_phi,A_phix,A_phiy,A_phiz,dX,dY,dZ,SYM,SYM,SYM)
  end if
  !
  ! Now compute B^i according to (exercise for the readers)
  !  B^x = (-x/pomega^2) e^(-6 phi) * A_{phi,z}; 
  !  B^y = (-y/pomega^2) e^(-6 phi) * A_{phi,z};
  !  B^z = e^(-6 phi) * (x A_{phi,x} + y A_{phi,y})/pomega^2; 
  !  pomega^2 = x^2 + y^2
  !
  ! and then calculate mhd_st_i and tau
  !
  do k = kmin,kmax
     do j = jmin,jmax
        do i=imin,imax
           psim6 = exp(-6.d0*phi(i,j,k))

           pomega2 = x(i,j,k)**2 + y(i,j,k)**2 
           Bx(i,j,k) = -x(i,j,k)/pomega2 * psim6 * A_phiz(i,j,k)
           By(i,j,k) = -y(i,j,k)/pomega2 * psim6 * A_phiz(i,j,k)
           Bz(i,j,k) = psim6/pomega2 * (x(i,j,k)*A_phix(i,j,k) +  & 
                y(i,j,k)*A_phiy(i,j,k))

           ! Renormalize the B-field
           Bx(i,j,k) = fac*Bx(i,j,k)
           By(i,j,k) = fac*By(i,j,k)
           Bz(i,j,k) = fac*Bz(i,j,k)


           ! Compute b^0 and b_i
           al = 1.d0 + alpham1(i,j,k)
           sqrtg = 1.d0/psim6
           sqrtg4 = al * sqrtg
           B2s = exp(4.d0*phi(i,j,k))*(gxx(i,j,k)*Bx(i,j,k)**2 + & 
                2.d0*gxy(i,j,k)*Bx(i,j,k)*By(i,j,k) + & 
                2.d0*gxz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) + &
                gyy(i,j,k)*By(i,j,k)**2 + 2.d0*gyz(i,j,k)*By(i,j,k)*Bz(i,j,k) + & 
                gzz(i,j,k)*Bz(i,j,k)**2)/(fs4pi*al)**2
           psin = exp(4.d0*phi(i,j,k))/al/fs4pi
           B_xs  = psin * (gxx(i,j,k) * Bx(i,j,k) + gxy(i,j,k) * By(i,j,k) + & 
                gxz(i,j,k) * Bz(i,j,k))
           B_ys  = psin * (gxy(i,j,k) * Bx(i,j,k) + gyy(i,j,k) * By(i,j,k) + & 
                gyz(i,j,k) * Bz(i,j,k))
           B_zs  = psin * (gxz(i,j,k) * Bx(i,j,k) + gyz(i,j,k) * By(i,j,k) + & 
                gzz(i,j,k) * Bz(i,j,k))
           sb0 = (u_x(i,j,k)*Bx(i,j,k) + u_y(i,j,k)*By(i,j,k) + &
                u_z(i,j,k)*Bz(i,j,k))/fs4pi/al
           sb2_l = (B2s + sb0**2)/u0(i,j,k)**2

           ! Depending on how you normalize the B-field, you may need to change this.
           !             sb2(i,j,k) = sb2_l/P(i,j,k)
           sb2(i,j,k) = sb2_l

           sb_x = (B_xs + u_x(i,j,k)*sb0)/u0(i,j,k)
           sb_y = (B_ys + u_y(i,j,k)*sb0)/u0(i,j,k)
           sb_z = (B_zs + u_z(i,j,k)*sb0)/u0(i,j,k)
           ! Now compute mhd_st_i and tau
           mhd_st_x(i,j,k) = st_x(i,j,k)+sqrtg4*(u0(i,j,k)*sb2_l*u_x(i,j,k)-sb0*sb_x)
           mhd_st_y(i,j,k) = st_y(i,j,k)+sqrtg4*(u0(i,j,k)*sb2_l*u_y(i,j,k)-sb0*sb_y)
           mhd_st_z(i,j,k) = st_z(i,j,k)+sqrtg4*(u0(i,j,k)*sb2_l*u_z(i,j,k)-sb0*sb_z)

           if (enable_HARM_energyvariable == 1) then
              u0u0 = -1.d0-u0(i,j,k)*(vx(i,j,k)*u_x(i,j,k) + &
                   vy(i,j,k)*u_y(i,j,k) + vz(i,j,k)*u_z(i,j,k) )
              sbx = Bx(i,j,k)/(al*fs4pi*u0(i,j,k)) + sb0*vx(i,j,k)
              sby = By(i,j,k)/(al*fs4pi*u0(i,j,k)) + sb0*vy(i,j,k)
              sbz = Bz(i,j,k)/(al*fs4pi*u0(i,j,k)) + sb0*vz(i,j,k)
              b0b0 = sb2_l - (sbx*sb_x + sby*sb_y + sbz*sb_z)
              tau(i,j,k) = tau(i,j,k) - sqrtg4*(sb2_l*(u0u0+0.5d0)-b0b0) 
           else 
              tau(i,j,k) = tau(i,j,k) + sqrtg*( sb2_l*(al*u0(i,j,k))**2 &
                   - sb2_l*0.5d0 - (al*sb0)**2 )
           end if

        end do
     end do
  end do


  if (Symmetry==AXISYM) then
     call axibc_vector(ext,X,Y,Z,Bx,By,Bz)
     call axibc_vector(ext,X,Y,Z,mhd_st_x,mhd_st_y,mhd_st_z)
     call axibc_scalar(ext,X,Y,Z,tau)
     call axibc_scalar(ext,X,Y,Z,sb2)
  end if

  ! compute physical B-fields
  call bi_phys(ext,X,Y,Z,Rp,dRdr,Bx,By,Bz,Bx_phys,By_phys,Bz_phys,Symmetry)

  ! try to compute the MRI wavelength.
  k = kmin + 1
  j = 2

  do i = imin+1,imax
     omega = vy(i,j,k)/X(i,j,k)

     lambda = 2.d0*sqrt(4*PI/15)*Bz_phys(i,j,k)/(omega*sqrt(rho_b(i,j,k) + 4.d0*P(i,j,k)))

! Uncomment below line to ENABLE output of MRI wavelength
!     write(*,*) 'x, \lambda = ', x(i,j,k), ' ', lambda
  end do


end subroutine initial_emfields_disk

!------------------------------------------------------------
! Compute the vector potential A_phi: 
!   A_phi = max(rhob/rhob_max - rho_cut, 0)
!------------------------------------------------------------
subroutine compute_Aphi_disk(ext,Z,r_in,rho_b,A_phi,rho_cut,rhob_max,Sym_Bz)
  implicit none
  integer, dimension(3)                          :: ext
  real*8, dimension(ext(1),ext(2),ext(3))        :: A_phi, rho_b, Z
  real*8						:: rho_cut
  real*8                                         :: rhob_max,r_in
  integer					:: i,j,k
  integer                                        :: kmin,kmax
  real*8                                         :: Sym_Bz
  !
  kmin = lbound(rho_b,3)
  kmax = ubound(rho_b,3)
  do k = kmin,kmax
     A_phi(:,:,k) = max(rho_b(:,:,k)/rhob_max - rho_cut, 0.d0)
     if (Sym_Bz .lt. 0.d0) A_phi(:,:,k) = A_phi(:,:,k)*z(:,:,k)/r_in
  end do

end subroutine compute_Aphi_disk

!------------------------------------------------------------
! Compute the vector potential A_phi:
!   A_phi = Omega                     (rho_b > rho_cut) 
!   A_phi = Omega * (rho_b/rho_cut)   (rho_b < rho_cut)
!------------------------------------------------------------
subroutine compute_Aphi_omega(ext,X,Z,rho_b,vy,A_phi,r_in,rho_cut,rhob_max,Sym_Bz)
  implicit none
  integer, dimension(3)                          :: ext
  real*8, dimension(ext(1),ext(2),ext(3))        :: A_phi, rho_b, vy, X, Z
  real*8                                         :: rho_cut
  real*8                                         :: rhob_max,r_in
  integer                                        :: i,j,k
  integer                                        :: kmin,kmax,imin,imax
  real*8                                         :: Sym_Bz
  !
  imin = lbound(rho_b,1)
  imax = ubound(rho_b,1)
  kmin = lbound(rho_b,3)
  kmax = ubound(rho_b,3)
  j = 2
  A_phi = 0.d0
  do k = kmin,kmax
     do i = imin,imax
        A_phi(i,j,k) = vy(i,j,k)/x(i,j,k)
        if (rho_b(i,j,k) .lt. rho_cut*rhob_max) &
             A_phi(i,j,k) = A_phi(i,j,k)*rho_b(i,j,k)/(rhob_max*rho_cut)
        if (Sym_Bz .lt. 0.d0) A_phi(i,j,k) = A_phi(i,j,k)*z(i,j,k)/r_in
     end do
  end do

end subroutine compute_Aphi_omega

subroutine trans_fish_phys_vj2(x, y, z, Rp, dR,  &
     vx,vy,vz,fish_to_phys)
  implicit none
  real*8                                   :: x, y, z
  real*8                                   :: Rp, dR
  real*8                                   :: vx,vy,vz
  real*8                                   :: r,fac, vx_t, vy_t, vz_t
  real*8                                   :: Jxx,Jxy,Jxz,Jyy,Jyz,Jzz
  logical				   :: fish_to_phys
  !
  r = sqrt(x**2 + y**2 + z**2)

  ! transformation matrix
  if (fish_to_phys) then
     fac = dR - Rp/r
     Jxx = Rp/r + x/r * x/r * fac
     Jxy = x/r * y/r * fac
     Jxz = x/r * z/r * fac;
     Jyy = Rp/r + y/r * y/r * fac
     Jyz = y/r * z/r * fac
     Jzz = Rp/r + z/r * z/r * fac
  else 
     fac = 1.d0/dR - r/Rp
     Jxx = r/Rp + x/r * x/r * fac
     Jxy = x/r * y/r * fac
     Jxz = x/r * z/r * fac;
     Jyy = r/Rp + y/r * y/r * fac
     Jyz = y/r * z/r * fac
     Jzz = r/Rp + z/r * z/r * fac
  end if

  vx_t = vx
  vy_t = vy
  vz_t = vz

  vx = Jxx*vx_t + Jxy*vy_t + Jxz*vz_t
  vy = Jxy*vx_t + Jyy*vy_t + Jyz*vz_t
  vz = Jxz*vx_t + Jyz*vy_t + Jzz*vz_t

end subroutine trans_fish_phys_vj2

!------------------------------------------------------------

FUNCTION ran3(idum)
  INTEGER idum
  INTEGER MBIG,MSEED,MZ
  !C     REAL MBIG,MSEED,MZ
  REAL ran3,FAC
  PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
  !C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
  INTEGER i,iff,ii,inext,inextp,k
  INTEGER mj,mk,ma(55)
  !C     REAL mj,mk,ma(55)
  SAVE iff,inext,inextp,ma
  DATA iff /0/
  if(idum.lt.0.or.iff.eq.0)then
     iff=1
     mj=MSEED-iabs(idum)
     mj=mod(mj,MBIG)
     ma(55)=mj
     mk=1
     do 11 i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if(mk.lt.MZ)mk=mk+MBIG
        mj=ma(ii)
11      continue
        do 13 k=1,4
           do 12 i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12            continue
13            continue
              inext=0
              inextp=31
              idum=1
           endif
           inext=inext+1
           if(inext.eq.56)inext=1
           inextp=inextp+1
           if(inextp.eq.56)inextp=1
           mj=ma(inext)-ma(inextp)
           if(mj.lt.MZ)mj=mj+MBIG
           ma(inext)=mj
           ran3=mj*FAC
           return
   END


