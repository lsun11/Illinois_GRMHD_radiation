!-----------------------------------------------------------------------------
! reconstruct primitive variables, compute sources for hybrid EOS
! This version supports a spatially varying density floor
! Also, this version allows for a hotter atmosphere.  Notice the
!               Pmax = 50.d0*P_cold
!  line of code below.
! Notice also that the initial guess for u0 is NOT 1.0, since we can do better
!
! USE THIS PRIMITIVES SOLVER FOR DISK RUNS
!-----------------------------------------------------------------------------
subroutine primitive_vars_alt_disk(ext, X, Y, Z, rho_star, tau, &
     st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
     neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
     w, rho_b, rho, P, h, Sx, Sy, Sz, &
     Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     phi, alpha, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
     h_old, u0, rho_max, rho_b_atm_gf, rho_fail_max_step, M_fail_step, &
     Bx, By, Bz, Ex, Ey, Ez, vx, vy, vz, sbt, sbx, sby, sbz, &
     proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
     glob_imax,glob_jmax,glob_kmax,Symmetry,Fontfix_tracker_gf,pfloor_gf,Nfont, &
     enable_HARM_energyvariable, excision_zone_gf, force_font_fix_fail, excision_enable)
  implicit none
! Hydro Input
  integer, dimension(3),intent(in)                    :: ext
  integer                                              :: neos,ergo_star
  real*8, dimension(neos)                              :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1)                              :: k_tab, gamma_tab
  real*8                                              :: gamma_th,ergo_sigma
  real*8, dimension(ext(1),ext(2),ext(3)), intent(in) :: X,Y,Z
  integer, dimension(ext(1),ext(2),ext(3)),intent(in) :: excision_zone_gf
  real*8, dimension(ext(1),ext(2),ext(3))             :: h_old
  real*8, dimension(ext(1),ext(2),ext(3))             :: rho_star,tau, Fontfix_tracker_gf, pfloor_gf
  real*8, dimension(ext(1),ext(2),ext(3))             :: st_x,st_y,st_z,rho_b_atm_gf
  real*8, dimension(ext(1),ext(2),ext(3))             :: mhd_st_x,mhd_st_y,mhd_st_z
  real*8                                              :: rho_max
  integer                                             :: proc_imin,proc_jmin,proc_kmin,Nfont
  integer                                             :: proc_imax,proc_jmax,proc_kmax
  integer                                             :: glob_imax,glob_jmax,glob_kmax
  integer, intent(in)                                      :: enable_HARM_energyvariable
  integer                                             :: force_font_fix_fail, excision_enable
! Output
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: w, rho_b, rho
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: P, h
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Sx, Sy, Sz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Sxx, Sxy, Sxz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Syy, Syz, Szz
  real*8, dimension(ext(1),ext(2),ext(3))             :: u0,vx,vy,vz
! Metric
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: phi,alpha
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: shiftx,shifty,shiftz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: gupxx, gupxy, gupxz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: gupyy, gupyz, gupzz
! MHD
  real*8, dimension(ext(1),ext(2),ext(3))              :: Bx,By,Bz,Ex,Ey,Ez
  real*8, dimension(ext(1),ext(2),ext(3))              :: sbt,sbx,sby,sbz
! Internal variables:
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  real*8                             :: P_l, eps, eps_star,u2,Psi6,Psi4,Psi2
  real*8                             :: rho_s, st_x_l, st_y_l,st_z_l, Tiny
  real*8                             :: deltE,B_xl,B_yl,B_zl
  real*8                             :: E_x_l,E_y_l,E_z_l,alpn1,Fx_l,Fy_l,Fz_l
  real*8                             :: ux_l,uy_l,uz_l, u_xl,u_yl,u_zl
  real*8                             :: Pmax,Pmin,rho_bl,h_l,P_tiny,BB
  real*8                             :: gijuiuj,au0m1,sb0,sb2,sb_x,sb_y,sb_z
  real*8                             :: fac, B2, dX, dY, dZ, dV
  real*8                             :: rho_fail_max_step, M_fail_step
  real*8                             :: P_cold, eps_cold, eps_tiny
  real*8, parameter                  :: TWO  = 2.D0
  real*8, parameter                  :: ZERO = 0.D0
  real*8, parameter                  :: ONE  = 1.D0
  real*8, parameter                  :: FOUR = 4.D0
  real*8, parameter                  :: SIX  = 6.D0
  real*8                             :: u_scal, eps_scal, sti_scal, tau_scal
  real*8                             :: f, rho_cutoff, factor
  real*8, parameter                     :: max_gamma = 60.d0
  real*8                             :: E_xl, E_yl, E_zl
  real*8                             :: u0l,temp,u_xll,u_yll,u_zll
!  real*8, dimension(ext(1),ext(2),ext(3)) :: E_x, E_y, E_z, B_x, B_y, B_z
!  real*8, dimension(ext(1),ext(2),ext(3)) :: u_x, u_y, u_z, u0
!  real*8, dimension(ext(1),ext(2),ext(3)) :: temp, psi_6, psi_4
  real*8 :: u_0n1,u_0,u_0hn1,u0u0,sb_0,b0b0
  real*8 :: PI,f1o4p, f1o8p,f1o4pa
  real*8, dimension(4) :: UU
  real*8, dimension(3) :: UU_font_fix
  integer, parameter :: m=31
  real*8, dimension(m) :: AUX
  integer :: Symmetry
  integer, parameter :: AXISYM = 4
  integer :: count,nn
  logical :: check, recom,exit_do, repairs_needed, ex_kmax, ex_imax
  logical, dimension(ext(1),ext(2),ext(3)) :: failures
  logical :: outside_excision, have_inner_boundary
  external funcv_hybrid_disk,fdjac_hybrid_disk
  external funcv_hybrid_font_fix2,fdjac_hybrid_font_fix2
  integer :: i_lower, k_upper, i_upper, k_lower
!
  count = 0
  PI = acos(-ONE)
  f1o4p = ONE/(FOUR*PI)
  f1o8p = f1o4p/TWO
!
! Input translation
!
  imin = lbound(gxx,1)
  jmin = lbound(gxx,2)
  kmin = lbound(gxx,3)
  imax = ubound(gxx,1)
  jmax = ubound(gxx,2)
  kmax = ubound(gxx,3)
  dX = X(imin+1,jmin,kmin)-X(imin,jmin,kmin)
  dY = Y(imin,jmin+1,kmin)-Y(imin,jmin,kmin)
  dZ = Z(imin,jmin,kmin+1)-Z(imin,jmin,kmin)
  i_upper = imax
  k_upper = kmax
  ! Lets go ahead and do the inversion on the outer boundary.
  ! testing
  if (proc_kmax .eq. glob_kmax) k_upper = kmax - 1
  if (proc_imax .eq. glob_imax) i_upper = imax - 1
  if (Symmetry ==AXISYM) then
     dV = 4.d0*acos(-1.d0)*dX*dZ
     ! only compute primitives and sources on y=0 plane
     jmin = 2
     jmax = 2
     rho = ZERO
     Sx  = ZERO
     Sy  = ZERO
     Sz  = ZERO
     Sxx = ZERO
     Sxy = ZERO
     Sxz = ZERO
     Syy = ZERO
     Syz = ZERO
     Szz = ZERO
  else
     dV = dX * dY * dZ
  end if
  rho_fail_max_step = ZERO
  M_fail_step =ZERO
  Fontfix_tracker_gf = ZERO
  Nfont = 0
!-----------------------------------------------------------------------------
! Funny parameters...  (See Shibata, Oohara & Nakamura...)
!-----------------------------------------------------------------------------
  f = 4.D-5
  rho_cutoff = 1D-4
  Tiny = ZERO
!  psi_4 = exp(4.d0*phi)
!  psi_6 = exp(6.d0*phi)
!  u_x = (gxx*(shiftx+vx) + gxy*(shifty+vy) + gxz*(shiftz+vz))*psi_4*u0
!  u_y = (gxy*(shiftx+vx) + gyy*(shifty+vy) + gyz*(shiftz+vz))*psi_4*u0
!  u_z = (gxz*(shiftx+vx) + gyz*(shifty+vy) + gzz*(shiftz+vz))*psi_4*u0
!  B_x  = psi_4 * (gxx * Bx + gxy * By + gxz * Bz)
!  B_y  = psi_4 * (gxy * Bx + gyy * By + gyz * Bz)
!  B_z  = psi_4 * (gxz * Bx + gyz * By + gzz * Bz)
! Initialize the failures array
  failures = .FALSE.
  repairs_needed = .FALSE.
  write(*,*) "HI1:",rho_max
!  write(*,*) HI2:,P_cold
!  write(*,*) HI3:,eps_scal
  write(*,*) "HI4:",neos
  write(*,*) "HI5:",rho_tab
! Set eps_scal
  call compute_pcold_epscold(rho_max, P_cold, eps_scal, &
       neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
!
! Now go to each gridpoint
!
  i_lower = imin+1
  k_lower = kmin+1
  do k = k_lower, k_upper
     do j = jmin, jmax
        do i = i_lower, i_upper
           ! Set eps_tiny and P_tiny
           call compute_pcold_epscold(rho_b_atm_gf(i,j,k), P_tiny, eps_tiny, &
                neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
           Psi2 = exp(2.d0*phi(i,j,k))
           Psi4 = Psi2*Psi2
           Psi6 = Psi2*Psi4
           alpn1 = alpha(i,j,k) + ONE
           u0l = u0(i,j,k)
           u_xll = (gxx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) + &
                    gxy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                    gxz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*Psi4*u0l
           u_yll = (gxy(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) + &
                    gyy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                    gyz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)))*Psi4*u0l
           u_zll = (gxz(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) + &
                    gyz(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                    gzz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)))*Psi4*u0l
           B_xl  = Psi4 * (gxx(i,j,k) * Bx(i,j,k) + gxy(i,j,k) * By(i,j,k) + &
                           gxz(i,j,k) * Bz(i,j,k) )
           B_yl  = Psi4 * (gxy(i,j,k) * Bx(i,j,k) + gyy(i,j,k) * By(i,j,k) + &
                           gyz(i,j,k) * Bz(i,j,k))
           B_zl  = Psi4 * (gxz(i,j,k) * Bx(i,j,k) + gyz(i,j,k) * By(i,j,k) + &
                           gzz(i,j,k) * Bz(i,j,k))
           if (rho_b_atm_gf(i,j,k)*Psi6*alpn1*u0(i,j,k) .lt. 1.001*rho_star(i,j,k)) then
              rho_s = rho_star(i,j,k)
           else
              rho_s = -1.d0
           end if
           f1o4pa = sqrt(f1o4p)/alpn1
           if (rho_s .gt. Tiny) then
              rho_bl = rho_b(i,j,k)
              eps = h_old(i,j,k)-1.d0-P(i,j,k)/rho_bl
              nn = 4
              AUX(1) = rho_s
              AUX(2) = tau(i,j,k)
              AUX(3) = gamma_th
              AUX(4) = Bx(i,j,k)*f1o4pa
              AUX(5) = By(i,j,k)*f1o4pa
              AUX(6) = Bz(i,j,k)*f1o4pa
              B2 = Psi4*( gxx(i,j,k)*AUX(4)**2 + &
                   2.d0*gxy(i,j,k)*AUX(4)*AUX(5) + 2.d0*gxz(i,j,k)*AUX(4)*AUX(6) + &
                   gyy(i,j,k)*AUX(5)**2 + 2.d0*gyz(i,j,k)*AUX(5)*AUX(6) + &
                   gzz(i,j,k)*AUX(6)**2 )
              AUX(7) = mhd_st_x(i,j,k)
              AUX(8) = mhd_st_y(i,j,k)
              AUX(9) = mhd_st_z(i,j,k)
              AUX(10) = alpha(i,j,k)
              AUX(11) = B2
              AUX(12) = B_xl*f1o4pa
              AUX(13) = B_yl*f1o4pa
              AUX(14) = B_zl*f1o4pa
              AUX(15) = gupxx(i,j,k)/Psi4
              AUX(16) = gupxy(i,j,k)/Psi4
              AUX(17) = gupxz(i,j,k)/Psi4
              AUX(18) = gupyy(i,j,k)/Psi4
              AUX(19) = gupyz(i,j,k)/Psi4
              AUX(20) = gupzz(i,j,k)/Psi4
              AUX(21) = Psi6
              AUX(22) = Psi2
              P_l = max(P(i,j,k),pfloor_gf(i,j,k))
              h_l = max(h_old(i,j,k),1.d0)
              u_scal  = sqrt( (P_l + B2)/(rho_bl*h_l+B2) )
              bb = Psi4*( gxx(i,j,k)*shiftx(i,j,k)**2 + &
                   2.d0*gxy(i,j,k)*shiftx(i,j,k)*shifty(i,j,k) + &
                   2.d0*gxz(i,j,k)*shiftx(i,j,k)*shiftz(i,j,k) + &
                   gyy(i,j,k)*shifty(i,j,k)**2 + &
                   2.d0*gyz(i,j,k)*shifty(i,j,k)*shiftz(i,j,k) + &
                   gzz(i,j,k)*shiftz(i,j,k)**2 )
              u_scal = u_scal + sqrt(bb)
              u_scal = u_scal*Psi2
              sti_scal = max(abs(mhd_st_x(i,j,k)), abs(mhd_st_y(i,j,k)), &
                   abs(mhd_st_z(i,j,k)), &
                   0.001d0*u_scal*Psi2*Psi6*(rho_bl*h_old(i,j,k)+B2) )
              tau_scal = max(abs(tau(i,j,k)), 0.01d0*Psi6*(rho_bl* &
                   eps_scal+0.5d0*B2) )
              UU(1) = u_xll/u_scal
              UU(2) = u_yll/u_scal
              UU(3) = u_zll/u_scal
              UU(4) = max(eps,eps_tiny)/eps_scal
              AUX(23) = u_scal
              AUX(24) = eps_scal
              AUX(25) = sti_scal
              AUX(26) = tau_scal
              AUX(27) = dble(enable_HARM_energyvariable)
              AUX(28) = shiftx(i,j,k)
              AUX(29) = shifty(i,j,k)
              AUX(30) = shiftz(i,j,k)
              AUX(31) = pfloor_gf(i,j,k)
              if(i==52 .and. j==2 .and. k==2) then
                 !              if(i==31 .and. j==31 .and. k==2) then
                 write(*,*) "inside primitives0", X(i,j,k),Y(i,j,k),Z(i,j,k)
                 write(*,*) "inside primitives0", rho_b(i,j,k),mhd_st_x(i,j,k),mhd_st_y(i,j,k),mhd_st_z(i,j,k),h_old(i,j,k),gupxy(i&
  &,j,k)
                 write(*,*) "inside primitives1",u_scal,eps_scal,sti_scal,tau_scal,tau(i,j,k),rho_bl,B2,Psi6,UU
              end if
              call newt2(UU,AUX,funcv_hybrid_disk,fdjac_hybrid_disk,nn,m, &
                   neos,ergo_star,ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab,check)
!              if(i==31 .and. j==31 .and. k==2) then
              if(i==52 .and. j==2 .and. k==2) then
                 write(*,*) "AFTER: AUX= ",u_scal,AUX,gxx(i,j,k),gxy(i,j,k),gxz(i,j,k)
                 write(*,*) "AFTER: UU = ",UU,check
              end if
              recom = .FALSE.
              ! Impose Font fix when the solver fails (check = .TRUE.) or when it
              !  gives negative eps [UU(4) < 0]
              !
              if(check .or. UU(4) .lt. ZERO) then
                 Fontfix_tracker_gf(i,j,k) = 1.0
                 if (Symmetry==AXISYM) then
                    if (j==2) then
                       M_fail_step = M_fail_step + rho_s * X(i,j,k)
                       if (rho_s/rho_max .gt. rho_fail_max_step)  &
                            rho_fail_max_step = rho_s/rho_max
                    end if
                 else
                    if (rho_s/rho_max .gt. rho_fail_max_step)  &
                         rho_fail_max_step = rho_s/rho_max
                    M_fail_step = M_fail_step + rho_s
                 end if
                 UU_font_fix(1) = 1.d0
                 UU_font_fix(2) = 1.d0
                 UU_font_fix(3) = 1.d0
                 if (k .ne. kmax) then
                    count = count + 1
                 end if
                 nn = 3
                 check = .FALSE.
                 call newt2(UU_font_fix, &
                      AUX,funcv_hybrid_font_fix2, &
                      fdjac_hybrid_font_fix2,nn,m, &
                      neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,check)
                 if (force_font_fix_fail == 1) check = .TRUE.
                 if(check) then
                    if (force_font_fix_fail == 0) then
                       write(*,*) 'Secondary solver also fails!'
                       write(*,*) 'Problem at (x,y,z) =', X(i,j,k),Y(i,j,k),Z(i,j,k),i,j,k
                       write(*,*) 'rho_star = ',rho_star(i,j,k)
                       write(*,*) 'tau = ',tau(i,j,k)
                       write(*,*) 'mhd_st_x = ',mhd_st_x(i,j,k)
                       write(*,*) 'mhd_st_y = ',mhd_st_y(i,j,k)
                       write(*,*) 'mhd_st_z = ',mhd_st_z(i,j,k)
                       write(*,*) 'Bx = ',Bx(i,j,k)
                       write(*,*) 'By = ',By(i,j,k)
                       write(*,*) 'Bz = ',Bz(i,j,k)
                       write(*,*) 'B2 = ',AUX(11)
                       write(*,*) 'gamma^xx = ',AUX(15)
                       write(*,*) 'gamma^xy = ',AUX(16)
                       write(*,*) 'gamma^xz = ',AUX(17)
                       write(*,*) 'gamma^yy = ',AUX(18)
                       write(*,*) 'gamma^yz = ',AUX(19)
                       write(*,*) 'gamma^zz = ',AUX(20)
                       write(*,*) 'exp(6 phi) = ',AUX(21)
                    end if
                    failures(i,j,k) = .TRUE.
                    repairs_needed = .TRUE.
                    ! Set everything to zero before calculating 
                    ! these quantities from averages.
                    ! testing
!                    u_xl         = ZERO
!                    u_yl         = ZERO
!                    u_zl         = ZERO
!                    rho_bl       = ZERO
!                    P_l          = ZERO
!                    eps          = ZERO
!                    u0(i,j,k)    = 1.d0/alpn1
                    recom        = .TRUE.
                    rho_bl = rho_b_atm_gf(i,j,k)
                    call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                    neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                    P_l = P_cold
                    if (P_l .lt. pfloor_gf(i,j,k)) then
                       P_l = pfloor_gf(i,j,k)
                    end if
                    eps = eps_cold + (P_l - P_cold)/(gamma_th * rho_bl)
                    u_xl         = ZERO
                    u_yl         = ZERO
                    u_zl         = ZERO
                 else
                    u_xl = UU_font_fix(1)*(UU_font_fix(1)**2 + 1.d0)
                    u_yl = UU_font_fix(2)*(UU_font_fix(2)**2 + 1.d0)
                    u_zl = UU_font_fix(3)*(UU_font_fix(3)**2 + 1.d0)
                    recom = .TRUE.
                    gijuiuj = AUX(15)*u_xl**2 +  &
                         2.d0*AUX(16)*u_xl*u_yl + 2.d0*AUX(17)*u_xl*u_zl + &
                         AUX(18)*u_yl**2 + 2.d0*AUX(19)*u_yl*u_zl + &
                         AUX(20)*u_zl**2
                    au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
                    ! *** Limit velocity 
                    if (au0m1 .gt. max_gamma-1.d0) then
                       fac = sqrt((max_gamma**2 - 1.d0)/((1.d0+au0m1)**2-1.d0))
                       u_xl = fac*u_xl
                       u_yl = fac*u_yl
                       u_zl = fac*u_zl
                       gijuiuj = gijuiuj * fac**2
                       au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
                    end if
                    u0(i,j,k) = (au0m1+1.d0)/alpn1
                    if (rho_s .lt. ZERO) u0(i,j,k) = -u0(i,j,k)
                    rho_bl = rho_s/alpn1/u0(i,j,k)/AUX(21)
                    call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                         neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                    P_l = P_cold
                    eps = eps_cold
                 end if
              else
                 u_xl = UU(1)*u_scal
                 u_yl = UU(2)*u_scal
                 u_zl = UU(3)*u_scal
                 eps  = UU(4)*eps_scal
                 gijuiuj = AUX(15)*u_xl**2 +  &
                      2.d0*AUX(16)*u_xl*u_yl + 2.d0*AUX(17)*u_xl*u_zl + &
                      AUX(18)*u_yl**2 + 2.d0*AUX(19)*u_yl*u_zl + &
                      AUX(20)*u_zl**2
                 au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
                 ! *** Limit velocity
                 if (au0m1 .gt. max_gamma-1.d0) then
                    recom = .TRUE.
                    fac = sqrt((max_gamma**2-1.d0)/((1.d0+au0m1)**2 - 1.d0))
                    u_xl = fac*u_xl
                    u_yl = fac*u_yl
                    u_zl = fac*u_zl
                    gijuiuj = gijuiuj * fac**2
                    au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
                 end if
                 u0(i,j,k) = (au0m1+1.d0)/alpn1
                 if (rho_s .lt. ZERO) u0(i,j,k) = -u0(i,j,k)
                 rho_bl = rho_s/alpn1/u0(i,j,k)/AUX(21)
                 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                      neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                 P_l = P_cold + rho_bl*(gamma_th-1.d0)*(eps-eps_cold)
              end if
              u_xll = u_xl
              u_yll = u_yl
              u_zll = u_zl
              h_l = 1.d0 + P_l/rho_bl + eps
              ux_l = -shiftx(i,j,k)*u0(i,j,k) + AUX(15)*u_xl + &
                   AUX(16)*u_yl + AUX(17)*u_zl
              uy_l = -shifty(i,j,k)*u0(i,j,k) + AUX(16)*u_xl + &
                   AUX(18)*u_yl + AUX(19)*u_zl
              uz_l = -shiftz(i,j,k)*u0(i,j,k) + AUX(17)*u_xl + &
                   AUX(19)*u_yl + AUX(20)*u_zl
              vx(i,j,k) = ux_l/u0(i,j,k)
              vy(i,j,k) = uy_l/u0(i,j,k)
              vz(i,j,k) = uz_l/u0(i,j,k)
              ! Limit P, and impose density floor 
              if (rho_bl .lt. rho_b_atm_gf(i,j,k)) then
                 recom = .TRUE.
                 rho_bl = rho_b_atm_gf(i,j,k)
                 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                      neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                 P_l = pfloor_gf(i,j,k)
                 eps = eps_cold + (Pmin-P_l)/(gamma_th-1.d0)/rho_bl
                 h_l = 1.d0 + P_l/rho_bl + eps
              end if
              call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                   neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
              Pmax = 50.d0*P_cold
              if((P_l .gt. Pmax) .and. (recom .eqv. .FALSE.)) then
                 recom = .TRUE.
                 P_l = Pmax
                 eps = eps_cold + (P_l-P_cold)/(gamma_th-1.d0)/rho_bl
                 h_l = 1.d0 + eps + P_l/rho_bl
              end if
              Pmin = pfloor_gf(i,j,k)
              if((P_l .lt. Pmin) .and. (recom.eqv..FALSE.)) then
                 recom = .TRUE.
                 P_l = Pmin
                 eps = eps_cold + (P_l-P_cold)/(gamma_th-1.d0)/rho_bl
                 h_l = 1.d0 + eps + P_l/rho_bl
              end if
              ! Re-compute rho_star, mhd_st_i and tau if necessary
              if (recom .eqv. .TRUE.) then
                 sb0 = u_xl*AUX(4) + u_yl*AUX(5) + u_zl*AUX(6)
                 sb2 = (AUX(11) + sb0**2)/u0(i,j,k)**2
                 sb_x = (AUX(12) + u_xl*sb0)/u0(i,j,k)
                 sb_y = (AUX(13) + u_yl*sb0)/u0(i,j,k)
                 sb_z = (AUX(14) + u_zl*sb0)/u0(i,j,k)
                 rho_s = alpn1*u0(i,j,k)*rho_bl*AUX(21)
                 rho_star(i,j,k) = rho_s
                 mhd_st_x(i,j,k) = rho_s*h_l*u_xl + &
                      alpn1*AUX(21)*u0(i,j,k)*sb2*u_xl - alpn1*AUX(21)*sb0*sb_x
                 mhd_st_y(i,j,k) = rho_s*h_l*u_yl + &
                      alpn1*AUX(21)*u0(i,j,k)*sb2*u_yl - alpn1*AUX(21)*sb0*sb_y
                 mhd_st_z(i,j,k) = rho_s*h_l*u_zl + &
                      alpn1*AUX(21)*u0(i,j,k)*sb2*u_zl - alpn1*AUX(21)*sb0*sb_z
                 if (enable_HARM_energyvariable.eq.1) then
                    u_0n1 = -alpha(i,j,k) - au0m1*alpn1 + &
                         u_xl*shiftx(i,j,k) + u_yl*shifty(i,j,k) + &
                         u_zl*shiftz(i,j,k)
                    u_0 = u_0n1 - 1.d0
                    u_0hn1 = u_0n1 + u_0*(P_l/rho_bl+eps)
                    u0u0 = u0(i,j,k)*u_0
                    sb_0 = -sb0*alpn1**2 + sb_x*shiftx(i,j,k) + &
                         sb_y*shifty(i,j,k) + sb_z*shiftz(i,j,k)
                    b0b0 = sb0*sb_0
                    tau(i,j,k) = -rho_s*u_0hn1 - alpn1*Psi6* ( &
                         P_l + (u0u0+0.5d0)*sb2 - b0b0 )
                 else
                    tau(i,j,k) = (au0m1+(P_l/rho_bl+eps)*alpn1*u0(i,j,k))*rho_s + &
                         AUX(21)*sb2*(alpn1*u0(i,j,k))**2 &
                         - AUX(21)*(P_l+sb2*0.5d0)-AUX(21)*(alpn1*sb0)**2
                 end if
              end if
              !
              ! hydro sources
              h(i,j,k) = h_l
              w(i,j,k) = alpn1*u0(i,j,k)*rho_s
              st_x_l = rho_s*h_l*u_xl
              st_y_l = rho_s*h_l*u_yl
              st_z_l = rho_s*h_l*u_zl
              st_x(i,j,k) = st_x_l
              st_y(i,j,k) = st_y_l
              st_z(i,j,k) = st_z_l
              fac  = ONE / ( Psi6 * w(i,j,k)* h(i,j,k) )
              rho_b(i,j,k) = rho_bl
              P(i,j,k)   = P_l
              rho(i,j,k) = h(i,j,k) * w(i,j,k) / Psi6 - P_l
              Sx(i,j,k)  = st_x_l / Psi6
              Sy(i,j,k)  = st_y_l / Psi6
              Sz(i,j,k)  = st_z_l / Psi6
              Sxx(i,j,k) = fac * st_x_l*st_x_l + Psi4 * gxx(i,j,k) * P_l
              Sxy(i,j,k) = fac * st_x_l*st_y_l + Psi4 * gxy(i,j,k) * P_l
              Sxz(i,j,k) = fac * st_x_l*st_z_l + Psi4 * gxz(i,j,k) * P_l
              Syy(i,j,k) = fac * st_y_l*st_y_l + Psi4 * gyy(i,j,k) * P_l
              Syz(i,j,k) = fac * st_y_l*st_z_l + Psi4 * gyz(i,j,k) * P_l
              Szz(i,j,k) = fac * st_z_l*st_z_l + Psi4 * gzz(i,j,k) * P_l
           else
              rho_bl = rho_b_atm_gf(i,j,k)
              call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                   neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
              P_l = P_cold
              eps = eps_cold
              if (P_l .lt. pfloor_gf(i,j,k)) then
                 eps = eps + (pfloor_gf(i,j,k)-P_l)/(rho_bl*(gamma_th - 1.d0))
                 P_l = pfloor_gf(i,j,k)
              end if
              h_l = 1.d0 + P_l/rho_bl + eps
              u0(i,j,k) = 1.d0/alpn1
              u_0 = -alpn1
              vx(i,j,k) = -shiftx(i,j,k)
              vy(i,j,k) = -shifty(i,j,k)
              vz(i,j,k) = -shiftz(i,j,k)
              B2 = Psi4*( gxx(i,j,k)*Bx(i,j,k)**2 + &
                   2.d0*gxy(i,j,k)*Bx(i,j,k)*By(i,j,k) + &
                   2.d0*gxz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) + &
                   gyy(i,j,k)*By(i,j,k)**2 + &
                   2.d0*gyz(i,j,k)*By(i,j,k)*Bz(i,j,k) + &
                   gzz(i,j,k)*Bz(i,j,k)**2 )
              sb2 = B2*f1o4p
              rho_star(i,j,k) = Psi6*rho_bl
              if (enable_HARM_energyvariable .eq. 0) then
                 tau(i,j,k)   = Psi6*(rho_bl*eps + sb2*0.5d0)
              else
                 tau(i,j,k)   = -rho_star(i,j,k)*(1.d0 + h_l*u_0) - alpn1*Psi6*P_l
              end if
              st_x(i,j,k)  = ZERO
              st_y(i,j,k)  = ZERO
              st_z(i,j,k)  = ZERO
              mhd_st_x(i,j,k) = ZERO
              mhd_st_y(i,j,k) = ZERO
              mhd_st_z(i,j,k) = ZERO
              u_xll   = ZERO
              u_yll   = ZERO
              u_zll   = ZERO
              h(i,j,k)     = h_l
              w(i,j,k)     = rho_star(i,j,k)
              rho_b(i,j,k) = rho_bl
              P(i,j,k)     = P_l
              rho(i,j,k)   = rho_bl*(1.d0+eps)
              Sx(i,j,k)    = ZERO
              Sy(i,j,k)    = ZERO
              Sz(i,j,k)    = ZERO
              Sxx(i,j,k)   = Psi4 * gxx(i,j,k) * P_l
              Sxy(i,j,k)   = Psi4 * gxy(i,j,k) * P_l
              Sxz(i,j,k)   = Psi4 * gxz(i,j,k) * P_l
              Syy(i,j,k)   = Psi4 * gyy(i,j,k) * P_l
              Syz(i,j,k)   = Psi4 * gyz(i,j,k) * P_l
              Szz(i,j,k)   = Psi4 * gzz(i,j,k) * P_l
           end if
           !
           ! MHD metric sources
           E_xl = Psi6 * ( By(i,j,k)*(vz(i,j,k)+shiftz(i,j,k)) &
                - Bz(i,j,k)*(vy(i,j,k)+shifty(i,j,k)) )/alpn1
           E_yl = Psi6 * ( Bz(i,j,k)*(vx(i,j,k)+shiftx(i,j,k)) &
                - Bx(i,j,k)*(vz(i,j,k)+shiftz(i,j,k)) )/alpn1
           E_zl = Psi6 * ( Bx(i,j,k)*(vy(i,j,k)+shifty(i,j,k)) &
                - By(i,j,k)*(vx(i,j,k)+shiftx(i,j,k)) )/alpn1
           Ex(i,j,k) = (gupxx(i,j,k)*E_xl + gupxy(i,j,k)*E_yl + &
                gupxz(i,j,k)*E_zl)/Psi4
           Ey(i,j,k) = (gupxy(i,j,k)*E_xl + gupyy(i,j,k)*E_yl + &
                gupyz(i,j,k)*E_zl)/Psi4
           Ez(i,j,k) = (gupxz(i,j,k)*E_xl + gupyz(i,j,k)*E_yl + &
                gupzz(i,j,k)*E_zl)/Psi4
           B_xl  = Psi4 * (gxx(i,j,k) * Bx(i,j,k) + gxy(i,j,k) * By(i,j,k) &
                + gxz(i,j,k) * Bz(i,j,k))
           B_yl  = Psi4 * (gxy(i,j,k) * Bx(i,j,k) + gyy(i,j,k) * By(i,j,k) &
                + gyz(i,j,k) * Bz(i,j,k))
           B_zl  = Psi4 * (gxz(i,j,k) * Bx(i,j,k) + gyz(i,j,k) * By(i,j,k) &
                + gzz(i,j,k) * Bz(i,j,k))
           temp = f1o8p*(Ex(i,j,k)*E_xl + Ey(i,j,k)*E_yl + Ez(i,j,k)*E_zl &
                + Bx(i,j,k)*B_xl + By(i,j,k)*B_yl + Bz(i,j,k)*B_zl)
           rho(i,j,k)   = rho(i,j,k) + temp
           Sxx(i,j,k)   = Sxx(i,j,k) + temp*Psi4*gxx(i,j,k) &
                - f1o4p*(E_xl*E_xl + B_xl*B_xl)
           Sxy(i,j,k)   = Sxy(i,j,k) + temp*Psi4*gxy(i,j,k) &
                - f1o4p*(E_xl*E_yl + B_xl*B_yl)
           Sxz(i,j,k)   = Sxz(i,j,k) + temp*Psi4*gxz(i,j,k) &
                - f1o4p*(E_xl*E_zl + B_xl*B_zl)
           Syy(i,j,k)   = Syy(i,j,k) + temp*Psi4*gyy(i,j,k) &
                - f1o4p*(E_yl*E_yl + B_yl*B_yl)
           Syz(i,j,k)   = Syz(i,j,k) + temp*Psi4*gyz(i,j,k) &
                - f1o4p*(E_yl*E_zl + B_yl*B_zl)
           Szz(i,j,k)   = Szz(i,j,k) + temp*Psi4*gzz(i,j,k) &
                - f1o4p*(E_zl*E_zl + B_zl*B_zl)
           Sx(i,j,k)    = Sx(i,j,k) + f1o4p*Psi6*(Ey(i,j,k)*Bz(i,j,k) - Ez(i,j,k)*By(i,j,k))
           Sy(i,j,k)    = Sy(i,j,k) + f1o4p*Psi6*(Ez(i,j,k)*Bx(i,j,k) - Ex(i,j,k)*Bz(i,j,k))
           Sz(i,j,k)    = Sz(i,j,k) + f1o4p*Psi6*(Ex(i,j,k)*By(i,j,k) - Ey(i,j,k)*Bx(i,j,k))
           sbt(i,j,k) = u_xll*Bx(i,j,k) + u_yll*By(i,j,k) + u_zll*Bz(i,j,k)
           sbx(i,j,k) = Bx(i,j,k)/u0(i,j,k) + vx(i,j,k)*sbt(i,j,k)
           sby(i,j,k) = By(i,j,k)/u0(i,j,k) + vy(i,j,k)*sbt(i,j,k)
           sbz(i,j,k) = Bz(i,j,k)/u0(i,j,k) + vz(i,j,k)*sbt(i,j,k)
        end do
     end do
  end do
  if (Symmetry.eq.AXISYM) then
     do i=1,3,2
        sbt(:,i,:) = ZERO
        sbx(:,i,:) = ZERO
        sby(:,i,:) = ZERO
        sbz(:,i,:) = ZERO
         Ex(:,i,:) = ZERO
         Ey(:,i,:) = ZERO
         Ez(:,i,:) = ZERO
      rho_b(:,i,:) = ZERO
          P(:,i,:) = ZERO
         vx(:,i,:) = ZERO
         vy(:,i,:) = ZERO
         vz(:,i,:) = ZERO
          h(:,i,:) = 1.d0
         u0(:,i,:) = 1.d0
   rho_star(:,i,:) = ZERO
       st_x(:,i,:) = ZERO
       st_y(:,i,:) = ZERO
       st_z(:,i,:) = ZERO
   mhd_st_x(:,i,:) = ZERO
   mhd_st_y(:,i,:) = ZERO
   mhd_st_z(:,i,:) = ZERO
        tau(:,i,:) = ZERO
         Bx(:,i,:) = ZERO
         By(:,i,:) = ZERO
         Bz(:,i,:) = ZERO
       Fontfix_tracker_gf(:,i,:) = ZERO
     end do
  end if
  M_fail_step = M_fail_step * dV
  Nfont = count
! Testing
!  if(repairs_needed) then
!     call repair_failures_mhd_alt(ext,X,Z,gamma_th, failures, rho_b, P, & 
!          vx, vy, vz, u0, w, h, rho_star, tau, &
!          st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, & 
!          rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
!          alpha, shiftx, shifty, shiftz, phi, &
!          gxx, gxy, gxz, gyy, gyz, gzz, &
!          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
!          sbt, sbx, sby, sbz, rho_b_atm_gf, & 
!          neos,ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
!          proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
!          glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor_gf,grid_type,enable_HARM_energyvariable)
!  end if
  if(excision_enable.eq.1) then
     !rho_b(:,1,:) = rho_b(:,2,:)
     !rho_b(:,3,:) = rho_b(:,2,:)
     call hydro_ezbc_hybrid(ext,X,Y,Z,rho_star,tau, &
          mhd_st_x,mhd_st_y,mhd_st_z,st_x,st_y,st_z, &
          rho_b,P,h,vx,vy,vz,w,&
          sbt,sbx,sby,sbz,Bx,By,Bz, &
          alpha,shiftx,shifty,shiftz,phi,&
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,&
          Symmetry,excision_zone_gf,gamma_th,neos,ergo_star,ergo_sigma,rho_tab,&
          P_tab,eps_tab,k_tab,gamma_tab);
     call remove_interior2(ext,X,Y,Z,sbt,excision_zone_gf,Symmetry);
     call remove_interior2(ext,X,Y,Z,sbx,excision_zone_gf,Symmetry);
     call remove_interior2(ext,X,Y,Z,sby,excision_zone_gf,Symmetry);
     call remove_interior2(ext,X,Y,Z,sbz,excision_zone_gf,Symmetry);
  end if
end subroutine primitive_vars_alt_disk
subroutine funcv_hybrid_disk(n,x,fvec,m,aux,neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
  implicit none
  integer                    :: n,m, ergo_star
  real*8, dimension(n)       :: x,fvec
  real*8, dimension(m)       :: aux
  real*8 :: rho_s,tau,h,gamma_th, ergo_sigma
  real*8 :: Bx,By,Bz,mhd_stx,mhd_sty,mhd_stz
  real*8 :: B_x,B_y,B_z, B2, betax,betay,betaz
  real*8 :: u_x,u_y,u_z,u0,eps,P,rhob,P_cold,eps_cold
  real*8 :: alp,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,sqrtg
  real*8 :: sb0,sb_x,sb_y,sb_z,sb2,gijuiuj,au0m1,Psi2
  real*8 :: u_scal,eps_scal,sti_scal,tau_scal,u_0n1,u0u0,b0b0,u_0h1
  real*8 :: alpm1, u_0, hm1, sb_0
  integer :: neos,enable_HARM_energyvariable
  logical :: exit_do
  real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab, gamma_tab
!
  enable_HARM_energyvariable = int(aux(27)+0.1d0)
  Psi2 = aux(22)
  u_scal = aux(23)
  eps_scal = aux(24)
  sti_scal = aux(25)
  tau_scal = aux(26)
  u_x    = x(1)*u_scal
  u_y    = x(2)*u_scal
  u_z    = x(3)*u_scal
  eps    = x(4)*eps_scal
  rho_s = aux(1)
  tau    = aux(2)
  gamma_th = aux(3)
  Bx    = aux(4)
  By    = aux(5)
  Bz    = aux(6)
  mhd_stx    = aux(7)
  mhd_sty    = aux(8)
  mhd_stz    = aux(9)
  alpm1 = aux(10)
  alp   = aux(10) + 1.d0
  B2  = aux(11)
  B_x  = aux(12)
  B_y  = aux(13)
  B_z  = aux(14)
  gupxx   = aux(15)
  gupxy   = aux(16)
  gupxz   = aux(17)
  gupyy   = aux(18)
  gupyz   = aux(19)
  gupzz   = aux(20)
  sqrtg = aux(21)
  betax = aux(28)
  betay = aux(29)
  betaz = aux(30)
  gijuiuj = gupxx*u_x**2 + 2.d0*gupxy*u_x*u_y + &
                2.d0*gupxz*u_x*u_z + gupyy*u_y**2 + 2.d0*gupyz*u_y*u_z + &
                gupzz*u_z**2
  au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
  if (rho_s .lt. 0.d0) au0m1 = gijuiuj/( 1.d0-sqrt(1.d0+gijuiuj) )
  u0 = (au0m1+1.d0)/alp
  rhob = abs(rho_s/alp/sqrtg/u0)
  call compute_pcold_epscold(rhob, P_cold, eps_cold, &
                neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
  P = P_cold + (gamma_th - 1.d0)*rhob*(eps - eps_cold)
  if (rhob .ne. 0.d0) then
     hm1   = P/rhob + eps
  else
     hm1   = 0.d0
  end if
  h = 1.d0 + hm1
  sb0 = u_x*Bx + u_y*By + u_z*Bz
  sb2 = (B2 + sb0**2)/u0**2
  sb_x = (B_x + u_x*sb0)/u0
  sb_y = (B_y + u_y*sb0)/u0
  sb_z = (B_z + u_z*sb0)/u0
!
! fvec(1): Eq. for mhd_st_x; fvec(2): Eq. for mhd_st_y; fvec(3): Eq. for mhd_st_z; 
! fvec(4): Eq. for tau.
!
  fvec(1) = (rho_s*h*u_x + alp*sqrtg*u0*sb2*u_x - alp*sqrtg*sb0*sb_x - mhd_stx)/sti_scal
  fvec(2) = (rho_s*h*u_y + alp*sqrtg*u0*sb2*u_y - alp*sqrtg*sb0*sb_y - mhd_sty)/sti_scal
  fvec(3) = (rho_s*h*u_z + alp*sqrtg*u0*sb2*u_z - alp*sqrtg*sb0*sb_z - mhd_stz)/sti_scal
  if (enable_HARM_energyvariable.eq.1) then
     u_0n1 = -alpm1 - alp*au0m1 + u_x*betax + u_y*betay + u_z*betaz
     u_0  = u_0n1 - 1.d0
     u_0h1 = u_0n1 + u_0*hm1
     u0u0 = u_0*u0
     sb_0 = -sb0*alp**2 + sb_x*betax + sb_y*betay + sb_z*betaz
     b0b0 = sb0*sb_0
     fvec(4) = -rho_s*u_0h1 - alp*sqrtg*(P+(u0u0+0.5d0)*sb2-b0b0) - tau
  else
    if (rhob .ne. 0.d0) then
       fvec(4) = (au0m1+alp*u0*(P/rhob+eps))*rho_s + sqrtg*sb2*(alp*u0)**2 &
        - sqrtg*(P+sb2*0.5d0)-sqrtg*(alp*sb0)**2 - tau
    else
       fvec(4) = sqrtg*sb2*(alp*u0)**2 &
        - sqrtg*(P+sb2*0.5d0)-sqrtg*(alp*sb0)**2 - tau
    end if
  end if
  fvec(4) = fvec(4)/tau_scal
end subroutine funcv_hybrid_disk
! Analytic Jacobian definition:
subroutine fdjac_hybrid_disk(n,x,m,aux,neos, ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                        fvec,NP,fjac)
  implicit none
  integer                :: n,m,NP,neos,i, enable_HARM_energyvariable, ergo_star
  real*8, dimension(n)   :: x,fvec
  real*8, dimension(m)   :: aux
  real*8, dimension(NP,NP) :: fjac
  real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab,gamma_tab
  real*8 :: rho_s,h,gamma_th,ergo_sigma
  real*8 :: Bx,By,Bz,mhd_stx,mhd_sty,mhd_stz
  real*8 :: B_x,B_y,B_z,B2,u_x,u_y,u_z,eps
  real*8 :: alp,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,sqrtg
  real*8 :: u0,sb0,sb2,sb_x,sb_y,sb_z,tau,P,rhob
  real*8 :: dPcold_drho, depscold_drho,P_cold,eps_cold
  real*8 :: du0dux,du0duy,du0duz, db2dux,db2duy,db2duz
  real*8 :: dbxdux,dbxduy,dbxduz, dbydux,dbyduy,dbyduz, dbzdux,dbzduy,dbzduz
  real*8 :: dpdux,dpduy,dpduz,gijuiuj,au0m1, c, Psi2
  real*8 :: drho_dux, drho_duy, drho_duz, dP_drho, dhdux,dhduy,dhduz
  real*8 :: du_0dux,du_0duy,du_0duz,db_0dux,db_0duy,db_0duz,tmp
  real*8 :: u_scal,eps_scal,sti_scal,tau_scal, betax,betay,betaz,alpm1
  real*8 :: u_0, hm1, sb_0
  logical :: exit_do
!
  enable_HARM_energyvariable = int(AUX(27)+0.1d0)
  Psi2   = aux(22)
  u_scal = aux(23)
  eps_scal = aux(24)
  sti_scal = aux(25)
  tau_scal = aux(26)
  u_x    = x(1)*u_scal
  u_y    = x(2)*u_scal
  u_z    = x(3)*u_scal
  eps    = x(4)*eps_scal
  rho_s  = aux(1)
  tau    = aux(2)
  gamma_th = aux(3)
  Bx    = aux(4)
  By    = aux(5)
  Bz    = aux(6)
  mhd_stx    = aux(7)
  mhd_sty    = aux(8)
  mhd_stz    = aux(9)
  alp  = aux(10) + 1.d0
  B2   = aux(11)
  B_x  = aux(12)
  B_y  = aux(13)
  B_z  = aux(14)
  gupxx   = aux(15)
  gupxy   = aux(16)
  gupxz   = aux(17)
  gupyy   = aux(18)
  gupyz   = aux(19)
  gupzz   = aux(20)
  sqrtg   = aux(21)
  betax   = aux(28)
  betay   = aux(29)
  betaz   = aux(30)
  gijuiuj = gupxx*u_x**2 + 2.d0*gupxy*u_x*u_y + &
                2.d0*gupxz*u_x*u_z + gupyy*u_y**2 + 2.d0*gupyz*u_y*u_z + &
                gupzz*u_z**2
  au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
  if (rho_s .lt. 0.d0) au0m1 = gijuiuj/( 1.d0-sqrt(1.d0+gijuiuj) )
  u0 = (au0m1+1.d0)/alp
  rhob = abs(rho_s/alp/sqrtg/u0)
  u_0 = u_x*betax + u_y*betay + u_z*betaz - u0*alp**2
  i = 1
  exit_do = .FALSE.
  do
     if (rhob .lt. rho_tab(i)) then
        exit_do = .TRUE.
        P_cold = k_tab(i)*rhob**gamma_tab(i)
        dPcold_drho = gamma_tab(i)*P_cold/rhob
        depscold_drho = P_cold/rhob**2
        if (i.eq.1) then
           eps_cold = P_cold/rhob/(gamma_tab(i)-1.d0)
        else
           eps_cold = eps_tab(i-1) +  &
                (P_cold/rhob - P_tab(i-1)/rho_tab(i-1))/(gamma_tab(i)-1.d0)
        end if
     end if
     if (i.eq.neos .or. exit_do) exit
     i = i + 1
  end do
  if (rhob .gt. rho_tab(neos)) then
     if (ergo_star .eq. 0) then
        P_cold = k_tab(neos+1)*rhob**gamma_tab(neos+1)
        eps_cold = eps_tab(neos) + (P_cold/rhob - P_tab(neos)/rho_tab(neos))/(gamma_tab(neos+1)-1.d0)
        dPcold_drho = gamma_tab(neos+1)*P_cold/rhob
        depscold_drho = P_cold/rhob**2
     else
        P_cold = ((ergo_sigma* (1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma) * rhob**(ergo_sigma+1) + P_ta&
  &b(neos) - ergo_sigma*(1+eps_tab(neos))*rho_tab(neos))/(ergo_sigma+1)
        eps_cold = (((1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma) * rhob**(ergo_sigma+1) - P_tab(neos) + &
  &ergo_sigma*((1+eps_tab(neos))*rho_tab(neos)))/((ergo_sigma+1)*rhob)-1
        dPcold_drho = ergo_sigma*(1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma*rhob
        depscold_drho = (ergo_sigma*(1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma*rhob**(ergo_sigma-1) + (P&
  &_tab(neos) - ergo_sigma*(1+eps_tab(neos))*rho_tab(neos))/rhob**2)/(ergo_sigma+1)
     end if
  end if
  P = P_cold + (gamma_th - 1.d0)*(eps - eps_cold)*rhob
  dP_drho = dPcold_drho + (gamma_th-1.d0)*(eps-eps_cold-rhob*depscold_drho)
  if (rhob .ne. 0.d0) then
     hm1   = P/rhob + eps
  else
     hm1   = 0.d0
  end if
  h = 1.d0 + hm1
  sb0 = u_x*Bx + u_y*By + u_z*Bz
  sb2 = (B2 + sb0**2)/u0**2
  sb_x = (B_x + u_x*sb0)/u0
  sb_y = (B_y + u_y*sb0)/u0
  sb_z = (B_z + u_z*sb0)/u0
  sb_0 = betax*sb_x + betay*sb_y + betaz*sb_z - sb0*alp**2
  du0dux = (gupxx*u_x + gupxy*u_y + gupxz*u_z)/u0/alp**2
  du0duy = (gupxy*u_x + gupyy*u_y + gupyz*u_z)/u0/alp**2
  du0duz = (gupxz*u_x + gupyz*u_y + gupzz*u_z)/u0/alp**2
  db2dux = 2.d0*sb0*Bx/u0**2 - 2.d0*sb2/u0*du0dux
  db2duy = 2.d0*sb0*By/u0**2 - 2.d0*sb2/u0*du0duy
  db2duz = 2.d0*sb0*Bz/u0**2 - 2.d0*sb2/u0*du0duz
  drho_dux = -rhob/u0*du0dux
  drho_duy = -rhob/u0*du0duy
  drho_duz = -rhob/u0*du0duz
  dpdux = dP_drho*drho_dux
  dpduy = dP_drho*drho_duy
  dpduz = dP_drho*drho_duz
  dhdux = dpdux/rhob - drho_dux*P/rhob**2
  dhduy = dpduy/rhob - drho_duy*P/rhob**2
  dhduz = dpduz/rhob - drho_duz*P/rhob**2
  dbxdux = (sb0+u_x*Bx)/u0 - sb_x/u0*du0dux
  dbxduy = u_x*By/u0 - sb_x/u0*du0duy
  dbxduz = u_x*Bz/u0 - sb_x/u0*du0duz
  dbydux = u_y*Bx/u0 - sb_y/u0*du0dux
  dbyduy = (sb0+u_y*By)/u0 - sb_y/u0*du0duy
  dbyduz = u_y*Bz/u0 - sb_y/u0*du0duz
  dbzdux = u_z*Bx/u0 - sb_z/u0*du0dux
  dbzduy = u_z*By/u0 - sb_z/u0*du0duy
  dbzduz = (sb0+u_z*Bz)/u0 - sb_z/u0*du0duz
! f(1) = mhd_st_x; f(2) = mhd_st_y; f(3) = mhd_st_z; f(4) = tau;
! x(1) = u_x/Psi^2; x(2) = u_y/Psi^2; x(3) = u_z/Psi^2; x(4) = eps;
! fjac(i,j) = partial f(i) / partial x(j) 
!
  c = rho_s*h + alp*sqrtg*u0*sb2
  fjac(1,1) = ( c + alp*sqrtg*( u_x*(u0*db2dux + sb2*du0dux) &
         - (Bx*sb_x + sb0*dbxdux) ) + rho_s*u_x*dhdux ) * u_scal/sti_scal
  fjac(1,2) = ( alp*sqrtg*( u_x*(u0*db2duy+sb2*du0duy)  &
         - (By*sb_x + sb0*dbxduy) ) + rho_s*u_x*dhduy ) * u_scal/sti_scal
  fjac(1,3) = ( alp*sqrtg*( u_x*(u0*db2duz+sb2*du0duz)  &
         - (Bz*sb_x + sb0*dbxduz) ) + rho_s*u_x*dhduz ) * u_scal/sti_scal
  fjac(1,4) = gamma_th*rho_s*u_x * eps_scal/sti_scal
  fjac(2,1) = ( alp*sqrtg*( u_y*(u0*db2dux+sb2*du0dux)  &
         - (Bx*sb_y + sb0*dbydux) ) + rho_s*u_y*dhdux ) * u_scal/sti_scal
  fjac(2,2) = ( c + alp*sqrtg*( u_y*(u0*db2duy+sb2*du0duy)  &
         - (By*sb_y + sb0*dbyduy) ) + rho_s*u_y*dhduy ) * u_scal/sti_scal
  fjac(2,3) = ( alp*sqrtg*( u_y*(u0*db2duz+sb2*du0duz)  &
         - (Bz*sb_y + sb0*dbyduz) ) + rho_s*u_y*dhduz ) * u_scal/sti_scal
  fjac(2,4) = gamma_th*rho_s*u_y * eps_scal/sti_scal
  fjac(3,1) = ( alp*sqrtg*( u_z*(u0*db2dux+sb2*du0dux)  &
         - (Bx*sb_z + sb0*dbzdux) ) + rho_s*u_z*dhdux ) * u_scal/sti_scal
  fjac(3,2) = ( alp*sqrtg*( u_z*(u0*db2duy+sb2*du0duy)  &
         - (By*sb_z + sb0*dbzduy) ) + rho_s*u_z*dhduy ) * u_scal/sti_scal
  fjac(3,3) = ( c + alp*sqrtg*( u_z*(u0*db2duz+sb2*du0duz)  &
         - (Bz*sb_z + sb0*dbzduz) ) + rho_s*u_z*dhduz ) * u_scal/sti_scal
  fjac(3,4) = gamma_th*rho_s*u_z * eps_scal/sti_scal
if (enable_HARM_energyvariable.eq.1) then
  du_0dux = betax - du0dux*alp**2
  du_0duy = betay - du0duy*alp**2
  du_0duz = betaz - du0duz*alp**2
  db_0dux = betax*dbxdux + betay*dbydux + betaz*dbzdux - Bx*alp**2
  db_0duy = betax*dbxduy + betay*dbyduy + betaz*dbzduy - By*alp**2
  db_0duz = betax*dbxduz + betay*dbyduz + betaz*dbzduz - Bz*alp**2
  fjac(4,1) = -( rho_s*(h*du_0dux + u_0*dhdux) + alp*sqrtg*( &
                (u0*du_0dux+u_0*du0dux)*sb2 + (u0*u_0+0.5d0)*db2dux &
                - sb0*db_0dux - sb_0*Bx + dpdux ) ) * u_scal/tau_scal
  fjac(4,2) = -( rho_s*(h*du_0duy + u_0*dhduy) + alp*sqrtg*( &
                (u0*du_0duy+u_0*du0duy)*sb2 + (u0*u_0+0.5d0)*db2duy &
                - sb0*db_0duy - sb_0*By + dpduy ) ) * u_scal/tau_scal
  fjac(4,3) = -( rho_s*(h*du_0duz + u_0*dhduz) + alp*sqrtg*( &
                (u0*du_0duz+u_0*du0duz)*sb2 + (u0*u_0+0.5d0)*db2duz &
                - sb0*db_0duz - sb_0*Bz + dpduz ) ) * u_scal/tau_scal
  fjac(4,4) = -rho_s*(u_0*gamma_th + (gamma_th-1.d0)/u0) * eps_scal/tau_scal
else
  fjac(4,1) = ( alp*rho_s*h*du0dux + alp*alp*sqrtg*(2.d0*u0*sb2*du0dux + &
                db2dux*u0**2 - 2.d0*sb0*Bx) - sqrtg*(dpdux + 0.5d0*db2dux) + &
                alp*rho_s*u0*dhdux ) * u_scal/tau_scal
  fjac(4,2) = ( alp*rho_s*h*du0duy + alp*alp*sqrtg*(2.d0*u0*sb2*du0duy + &
                db2duy*u0**2 - 2.d0*sb0*By) - sqrtg*(dpduy + 0.5d0*db2duy) + &
                alp*rho_s*u0*dhduy ) * u_scal/tau_scal
  fjac(4,3) = ( alp*rho_s*h*du0duz + alp*alp*sqrtg*(2.d0*u0*sb2*du0duz + &
                db2duz*u0**2 - 2.d0*sb0*Bz) - sqrtg*(dpduz + 0.5d0*db2duz) + &
                alp*rho_s*u0*dhduz ) * u_scal/tau_scal
  fjac(4,4) = ( alp*u0*gamma_th*rho_s - (gamma_th-1.d0)*rho_s/(alp*u0) ) * &
                eps_scal/tau_scal
end if
end subroutine fdjac_hybrid_disk
! This function is modified to take into account varying pressure
! and density floors.
subroutine funcv_hybrid_font_fix2(n,x,fvec,m,aux,neos,ergo_star, ergo_sigma, rho_tab, &
                        P_tab,eps_tab,k_tab,gamma_tab)
  implicit none
  integer                    :: n,m, ergo_star
  real*8, dimension(n)       :: x,fvec
  real*8, dimension(m)       :: aux
  real*8 :: rho_s,tau,h,gamma_th,ergo_sigma
  real*8 :: Bx,By,Bz,mhd_stx,mhd_sty,mhd_stz
  real*8 :: B_x,B_y,B_z, B2
  real*8 :: u_x,u_y,u_z,u0,eps,P,rhob,pfloor
  real*8 :: alp,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,sqrtg
  real*8 :: sb0,sb_x,sb_y,sb_z,sb2,gijuiuj,au0m1,Psi2, sti_scal
  integer :: neos
  logical exit_do
  real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab, gamma_tab
!
  Psi2   = aux(22)
  sti_scal = aux(25)
  u_x    = x(1)*(x(1)**2 + 1.d0)*Psi2
  u_y    = x(2)*(x(2)**2 + 1.d0)*Psi2
  u_z    = x(3)*(x(3)**2 + 1.d0)*Psi2
  rho_s = aux(1)
  tau    = aux(2)
  gamma_th = aux(3)
  Bx    = aux(4)
  By    = aux(5)
  Bz    = aux(6)
  mhd_stx    = aux(7)
  mhd_sty    = aux(8)
  mhd_stz    = aux(9)
  alp   = aux(10)+1.d0
  B2  = aux(11)
  B_x  = aux(12)
  B_y  = aux(13)
  B_z  = aux(14)
  gupxx   = aux(15)
  gupxy   = aux(16)
  gupxz   = aux(17)
  gupyy   = aux(18)
  gupyz   = aux(19)
  gupzz   = aux(20)
  sqrtg = aux(21)
  pfloor = aux(31)
  gijuiuj = gupxx*u_x**2 + 2.d0*gupxy*u_x*u_y + &
                2.d0*gupxz*u_x*u_z + gupyy*u_y**2 + 2.d0*gupyz*u_y*u_z + &
                gupzz*u_z**2
  au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
  if (rho_s .lt. 0.d0) au0m1 = gijuiuj/( 1.d0-sqrt(1.d0+gijuiuj) )
  u0 = (au0m1+1.d0)/alp
  rhob = rho_s/alp/sqrtg/u0
  call compute_pcold_epscold(rhob, P, eps, &
                neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
  if (P .lt. pfloor) then
     eps = eps + (pfloor-P)/(rhob*(gamma_th-1.d0))
     P = pfloor
  end if
  h     = 1.d0 + P/rhob + eps
  sb0 = u_x*Bx + u_y*By + u_z*Bz
  sb2 = (B2 + sb0**2)/u0**2
  sb_x = (B_x + u_x*sb0)/u0
  sb_y = (B_y + u_y*sb0)/u0
  sb_z = (B_z + u_z*sb0)/u0
!
! fvec(1): Eq. for mhd_st_x; fvec(2): Eq. for mhd_st_y; fvec(3): Eq. for mhd_st_z; 
!
  fvec(1) = (rho_s*h*u_x + alp*sqrtg*u0*sb2*u_x - alp*sqrtg*sb0*sb_x - mhd_stx)/sti_scal
  fvec(2) = (rho_s*h*u_y + alp*sqrtg*u0*sb2*u_y - alp*sqrtg*sb0*sb_y - mhd_sty)/sti_scal
  fvec(3) = (rho_s*h*u_z + alp*sqrtg*u0*sb2*u_z - alp*sqrtg*sb0*sb_z - mhd_stz)/sti_scal
end subroutine funcv_hybrid_font_fix2
! This function is modified to take into account varying pressure
! and density floors.
subroutine fdjac_hybrid_font_fix2(n,x,m,aux,neos, ergo_star, ergo_sigma,rho_tab,P_tab, &
                eps_tab,k_tab,gamma_tab,fvec,NP,fjac)
  implicit none
  integer                :: n,m,NP,neos,i,ergo_star
  real*8, dimension(n)   :: x,fvec
  real*8, dimension(m)   :: aux
  real*8, dimension(NP,NP) :: fjac
  real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab,gamma_tab
  real*8, dimension(3) :: fac
  real*8 :: rho_s,h,gamma_th,ergo_sigma
  real*8 :: Bx,By,Bz,mhd_stx,mhd_sty,mhd_stz
  real*8 :: B_x,B_y,B_z,B2,u_x,u_y,u_z,eps
  real*8 :: alp,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,sqrtg
  real*8 :: u0,sb0,sb2,sb_x,sb_y,sb_z,tau,P,rhob
  real*8 :: du0dux,du0duy,du0duz, db2dux,db2duy,db2duz
  real*8 :: dbxdux,dbxduy,dbxduz, dbydux,dbyduy,dbyduz, dbzdux,dbzduy,dbzduz
  real*8 :: dpdux,dpduy,dpduz,gijuiuj,au0m1, c, pfloor
  real*8 :: dhdux, dhduy, dhduz, temp, Psi2, sti_scal
  real*8 :: dP_drho, deps_drho, drho_dux, drho_duy, drho_duz
  logical :: exit_do
!
  Psi2   = aux(22)
  sti_scal = aux(25)
  u_x    = x(1)*(x(1)**2 + 1.d0)*Psi2
  u_y    = x(2)*(x(2)**2 + 1.d0)*Psi2
  u_z    = x(3)*(x(3)**2 + 1.d0)*Psi2
  fac(1) = (3.d0*x(1)**2 + 1.d0)*Psi2/sti_scal
  fac(2) = (3.d0*x(2)**2 + 1.d0)*Psi2/sti_scal
  fac(3) = (3.d0*x(3)**2 + 1.d0)*Psi2/sti_scal
  rho_s  = aux(1)
  tau    = aux(2)
  gamma_th = aux(3)
  Bx    = aux(4)
  By    = aux(5)
  Bz    = aux(6)
  mhd_stx    = aux(7)
  mhd_sty    = aux(8)
  mhd_stz    = aux(9)
  alp  = aux(10)+1.d0
  B2   = aux(11)
  B_x  = aux(12)
  B_y  = aux(13)
  B_z  = aux(14)
  gupxx   = aux(15)
  gupxy   = aux(16)
  gupxz   = aux(17)
  gupyy   = aux(18)
  gupyz   = aux(19)
  gupzz   = aux(20)
  sqrtg   = aux(21)
  pfloor  = aux(31)
  gijuiuj = gupxx*u_x**2 + 2.d0*gupxy*u_x*u_y + &
                2.d0*gupxz*u_x*u_z + gupyy*u_y**2 + 2.d0*gupyz*u_y*u_z + &
                gupzz*u_z**2
  au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
  if (rho_s .lt. 0.d0) au0m1 = gijuiuj/( 1.d0-sqrt(1.d0+gijuiuj) )
  u0 = (au0m1+1.d0)/alp
  rhob = rho_s/alp/sqrtg/u0
  i = 1
  exit_do = .FALSE.
  do
     if (rhob .lt. rho_tab(i)) then
        exit_do = .TRUE.
        P = k_tab(i)*rhob**gamma_tab(i)
        dP_drho = gamma_tab(i)*P/rhob
        deps_drho = P/rhob**2
        if (i.eq.1) then
           eps = P/rhob/(gamma_tab(i)-1.d0)
        else
           eps = eps_tab(i) +  &
                (P/rhob - P_tab(i-1)/rho_tab(i-1))/(gamma_tab(i)-1.d0)
        end if
     end if
     if (i.eq.neos .or. exit_do) exit
     i = i + 1
  end do
  if (rhob .gt. rho_tab(neos)) then
     P = k_tab(neos+1)*rhob**gamma_tab(neos+1)
     eps = eps_tab(neos) + (P/rhob - P_tab(neos)/rho_tab(neos))
     dP_drho = gamma_tab(neos+1)*P/rhob
     deps_drho = P/rhob**2
  end if
  !before going on, you must gaurantee that P is above the floor value.
  ! this will not work for the hybrid eos
  if (P .lt. pfloor) then
     eps = eps + (pfloor-P)/(rhob*(gamma_th-1.d0))
     P = pfloor
!     dP_drho = gamma_tab(neos+1)*P/rhob
!     deps_drho = P/rhob**2
  end if
  h     = 1.d0 + P/rhob + eps
  sb0 = u_x*Bx + u_y*By + u_z*Bz
  sb2 = (B2 + sb0**2)/u0**2
  sb_x = (B_x + u_x*sb0)/u0
  sb_y = (B_y + u_y*sb0)/u0
  sb_z = (B_z + u_z*sb0)/u0
  du0dux = (gupxx*u_x + gupxy*u_y + gupxz*u_z)/u0/alp**2
  du0duy = (gupxy*u_x + gupyy*u_y + gupyz*u_z)/u0/alp**2
  du0duz = (gupxz*u_x + gupyz*u_y + gupzz*u_z)/u0/alp**2
  db2dux = 2.d0*sb0*Bx/u0**2 - 2.d0*sb2/u0*du0dux
  db2duy = 2.d0*sb0*By/u0**2 - 2.d0*sb2/u0*du0duy
  db2duz = 2.d0*sb0*Bz/u0**2 - 2.d0*sb2/u0*du0duz
  temp = -rhob/u0*dP_drho
  dpdux = temp*du0dux
  dpduy = temp*du0duy
  dpduz = temp*du0duz
  temp = P/rhob/u0 - deps_drho*rhob/u0
  dhdux = temp*du0dux + dpdux/rhob
  dhduy = temp*du0duy + dpduy/rhob
  dhduz = temp*du0duz + dpduz/rhob
  dbxdux = (sb0+u_x*Bx)/u0 - sb_x/u0*du0dux
  dbxduy = u_x*By/u0 - sb_x/u0*du0duy
  dbxduz = u_x*Bz/u0 - sb_x/u0*du0duz
  dbydux = u_y*Bx/u0 - sb_y/u0*du0dux
  dbyduy = (sb0+u_y*By)/u0 - sb_y/u0*du0duy
  dbyduz = u_y*Bz/u0 - sb_y/u0*du0duz
  dbzdux = u_z*Bx/u0 - sb_z/u0*du0dux
  dbzduy = u_z*By/u0 - sb_z/u0*du0duy
  dbzduz = (sb0+u_z*Bz)/u0 - sb_z/u0*du0duz
!
! f(1) = mhd_st_x; f(2) = mhd_st_y; f(3) = mhd_st_z; 
! u_i= Psi^2 [ x(i)^3 + x(i) ]
! fjac(i,j) = partial f(i) / partial x(j) 
!
  c = rho_s*h + alp*sqrtg*u0*sb2
  fjac(1,1) = ( c + alp*sqrtg*( u_x*(u0*db2dux + sb2*du0dux) &
         - (Bx*sb_x + sb0*dbxdux) ) + rho_s*u_x*dhdux ) * fac(1)
  fjac(1,2) = ( alp*sqrtg*( u_x*(u0*db2duy+sb2*du0duy)  &
         - (By*sb_x + sb0*dbxduy) ) + rho_s*u_x*dhduy ) * fac(2)
  fjac(1,3) = ( alp*sqrtg*( u_x*(u0*db2duz+sb2*du0duz)  &
         - (Bz*sb_x + sb0*dbxduz) ) + rho_s*u_x*dhduz ) * fac(3)
  fjac(2,1) = ( alp*sqrtg*( u_y*(u0*db2dux+sb2*du0dux)  &
         - (Bx*sb_y + sb0*dbydux) ) + rho_s*u_y*dhdux ) * fac(1)
  fjac(2,2) = ( c + alp*sqrtg*( u_y*(u0*db2duy+sb2*du0duy)  &
         - (By*sb_y + sb0*dbyduy) ) + rho_s*u_y*dhduy ) * fac(2)
  fjac(2,3) = ( alp*sqrtg*( u_y*(u0*db2duz+sb2*du0duz)  &
         - (Bz*sb_y + sb0*dbyduz) ) + rho_s*u_y*dhduz ) * fac(3)
  fjac(3,1) = ( alp*sqrtg*( u_z*(u0*db2dux+sb2*du0dux)  &
         - (Bx*sb_z + sb0*dbzdux) ) + rho_s*u_z*dhdux ) * fac(1)
  fjac(3,2) = ( alp*sqrtg*( u_z*(u0*db2duy+sb2*du0duy)  &
         - (By*sb_z + sb0*dbzduy) ) + rho_s*u_z*dhduy ) * fac(2)
  fjac(3,3) = ( c + alp*sqrtg*( u_z*(u0*db2duz+sb2*du0duz)  &
         - (Bz*sb_z + sb0*dbzduz) ) + rho_s*u_z*dhduz ) * fac(3)
end subroutine fdjac_hybrid_font_fix2
subroutine repair_failures_mhd_alt(ext,X,Z, gamma_th, failures, rho_b, P, &
     vx, vy, vz, u0, w, h, rho_star, tau, &
     st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
     rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     alpha, betax, betay, betaz, phi, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
     sbt, sbx, sby, sbz, rho_b_atm_gf, &
     neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
     proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
     glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor_gf,energy_type)
  implicit none
!
! Input
!
  integer, dimension(3),intent(in)                :: ext
  real*8                                          :: n, gamma_th, ergo_sigma
  logical, dimension(ext(1),ext(2),ext(3))        :: failures
  real*8, dimension(ext(1),ext(2),ext(3))          :: rho_b, P, vx, vy, vz, X, Z
  real*8, dimension(ext(1),ext(2),ext(3))          :: u0, w, h, rho_star, tau
  real*8, dimension(ext(1),ext(2),ext(3))         :: st_x, st_y, st_z, pfloor_gf
  real*8, dimension(ext(1),ext(2),ext(3))         :: mhd_st_x, mhd_st_y, mhd_st_z
  real*8, dimension(ext(1),ext(2),ext(3))         :: Bx, By, Bz, rho_b_atm_gf
  real*8, dimension(ext(1),ext(2),ext(3))         :: sbt, sbx, sby, sbz
  integer                                           :: neos, ergo_star
  real*8, dimension(neos)                           :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1)                           :: k_tab, gamma_tab
  integer                                         :: proc_imin,proc_jmin,proc_kmin
  integer                                         :: proc_imax,proc_jmax,proc_kmax
  integer                                         :: glob_imax,glob_jmax,glob_kmax
  integer                                         :: energy_type
!
! Sources
!
  real*8, dimension(ext(1),ext(2),ext(3))            :: rho, Sx, Sy, Sz
  real*8, dimension(ext(1),ext(2),ext(3))            :: Sxx, Sxy, Sxz, Syy, Syz, Szz
!
! Metric
!
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: alpha, betax, betay, betaz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: gxx, gxy, gxz, phi
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: gyy, gyz, gzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: gupxx, gupxy, gupxz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: gupyy, gupyz, gupzz
!
! Other variables:
! 
  real*8, dimension(ext(1),ext(2),ext(3))  :: temp, temp2
  integer                               :: i,j,k
  real*8                                :: P_l, Psi6, Psi4
  real*8                                :: rho_s, st_x_l, st_y_l, st_z_l, z_l
  real*8                                :: w_l, tau_l, fac, Pmax, Pmin
  real*8                                :: alpha_l, temp3
  real*8                                 :: rho_bl, h_l, eps, P_cold, eps_cold
  real*8                                :: vx_l, vy_l, vz_l
  real*8                                :: u_xl, u_yl, u_zl
  real*8                                :: B_xl, B_yl, B_zl, B2
  real*8                                :: E_x, E_y, E_z, Ex, Ey, Ez
  real*8                                :: sb0, sb2, sb_x, sb_y, sb_z
  real*8, parameter                     :: TWO  = 2.D0
  real*8, parameter                     :: ZERO = 0.D0
  real*8, parameter                     :: ONE  = 1.D0
  real*8, parameter                     :: FOUR = 4.D0
  real*8, parameter                     :: SIX  = 6.D0
  integer                               :: point_count,ii
  real*8                                :: rho_b_fix_max
  integer                               :: imin, imax, jmin, jmax, kmin, kmax
  real*8                                :: alpn1, f1o4pa, f1o4p, f1o8p, PI
  real*8                                :: er, el, au0r1
  real*8, parameter                     :: fac2 = 0.99D0
  logical                                :: exit_do
  real*8                                :: u_0n1, u_0, u_0hn1, u0u0, sb_0, b0b0
  integer :: Symmetry
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  logical :: x_ogb_p,x_ogb_n,y_ogb_p,y_ogb_n,z_ogb_p,z_ogb_n, outside_excision
  logical :: have_inner_boundary
  integer :: i_lower
  real*8  :: dX
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!
  PI = acos(-ONE)
  f1o4p = ONE/(FOUR*PI)
  f1o8p = f1o4p/TWO
!
! Input translation
!
  imin = lbound(gxx,1)
  jmin = lbound(gxx,2)
  kmin = lbound(gxx,3)
  imax = ubound(gxx,1)
  jmax = ubound(gxx,2)
  kmax = ubound(gxx,3)
  dX = X(imin+1,jmin,kmin) - X(imin,jmin,kmin)
  outside_excision = .true.
  have_inner_boundary = .false.
  ! The following lines will prevent the ghost zones
  ! from being used in the averaging.
  if (proc_imin .eq. 0) imin = imin + 1
  if (proc_kmin .eq. 0) kmin = kmin + 1
!
! Go to each gridpoint
!
  i_lower = imin
  do k=kmin, kmax
     do j=jmin, jmax
        do i=i_lower, imax
           alpn1 = alpha(i,j,k) + ONE
           f1o4pa = sqrt(f1o4p)/alpn1
           outer: if (failures(i,j,k)) then
              rho_bl = 0.d0
              P_l = 0.d0
              vx_l = 0.d0
              vy_l = 0.d0
              vz_l = 0.d0
              point_count = 0
              ! ogb = outer grid boundary, this can be at positive (p) or negative (n) x
              x_ogb_p = (i+1 .eq. imax) .and. (proc_imax .eq. glob_imax)
                 x_ogb_n = (i-1 .eq. imin) .and. (proc_imin .eq. 0) .and. &
                      ((SYMMETRY .eq. NO_SYMM) .or. &
                      (SYMMETRY .eq. EQUATORIAL) .or. (SYMMETRY .eq. PI_SYMM))
              y_ogb_p = (j+1 .eq. jmax) .and. (proc_jmax .eq. glob_jmax)
              y_ogb_n = (j-1 .eq. jmin) .and. (proc_jmin .eq. 0) .and. &
                        ((SYMMETRY .eq. NO_SYMM) .or. (SYMMETRY .eq. EQUATORIAL))
              z_ogb_p = (k+1 .eq. kmax) .and. (proc_kmax .eq. glob_kmax)
              z_ogb_n = (k-1 .eq. kmin) .and. (proc_kmin .eq. 0) .and. &
                   ((SYMMETRY .eq. NO_SYMM) .or. (Z(imin,jmin,kmin) .lt. 0.d0))
              xmax: if ((i .ne. imax) .and. (i .ne. imin)) then
                 p1:  if ( (failures(i+1,j,k) .eqv. .false.) .and. &
                           (failures(i-1,j,k) .eqv. .false.) .and. &
                           (.not. x_ogb_p) .and. (.not. x_ogb_n) ) then
                    rho_bl = rho_bl + rho_b(i+1,j,k) + rho_b(i-1,j,k)
                    P_l = P_l + P(i+1,j,k) + P(i-1,j,k)
                    vx_l = vx_l + vx(i+1,j,k) + vx(i-1,j,k)
                    vy_l = vy_l + vy(i+1,j,k) + vy(i-1,j,k)
                    vz_l = vz_l + vz(i+1,j,k) + vz(i-1,j,k)
                    point_count = point_count + 2
                 end if p1
                 if ((j .ne. jmax) .and. (j .ne. jmin)) then
                    p2: if ((failures(i+1,j+1,k) .eqv. .false.) .and. &
                            (failures(i-1,j-1,k) .eqv. .false.) .and. &
                            (.not. x_ogb_p) .and. (.not. x_ogb_n) .and. &
                            (.not. y_ogb_p) .and. (.not. y_ogb_n) ) then
                       rho_bl = rho_bl + rho_b(i+1,j+1,k) + rho_b(i-1,j-1,k)
                       P_l = P_l + P(i+1,j+1,k) + P(i-1,j-1,k)
                       vx_l = vx_l + vx(i+1,j+1,k) + vx(i-1,j-1,k)
                       vy_l = vy_l + vy(i+1,j+1,k) + vy(i-1,j-1,k)
                       vz_l = vz_l + vz(i+1,j+1,k) + vz(i-1,j-1,k)
                       point_count = point_count + 2
                    end if p2
                    p3: if ((failures(i-1,j+1,k) .eqv. .false.) .and. &
                            (failures(i+1,j-1,k) .eqv. .false.) .and. &
                            (.not. x_ogb_p) .and. (.not. x_ogb_n) .and. &
                            (.not. y_ogb_p) .and. (.not. y_ogb_n) ) then
                       rho_bl = rho_bl + rho_b(i-1,j+1,k) + rho_b(i+1,j-1,k)
                       P_l = P_l + P(i-1,j+1,k) + P(i+1,j-1,k)
                       vx_l = vx_l + vx(i-1,j+1,k) + vx(i+1,j-1,k)
                       vy_l = vy_l + vy(i-1,j+1,k) + vy(i+1,j-1,k)
                       vz_l = vz_l + vz(i-1,j+1,k) + vz(i+1,j-1,k)
                       point_count = point_count + 2
                    end if p3
                    if ((k .ne. kmax) .and. (k .ne. kmin)) then
                       p4: if ((failures(i-1,j-1,k-1) .eqv. .false.) .and. &
                               (failures(i+1,j+1,k+1) .eqv. .false.) .and. &
                               (.not. x_ogb_p) .and. (.not. x_ogb_n) .and. &
                               (.not. y_ogb_p) .and. (.not. y_ogb_n) .and. &
                               (.not. z_ogb_p) .and. (.not. z_ogb_n) ) then
                          rho_bl = rho_bl + rho_b(i-1,j-1,k-1) + rho_b(i+1,j+1,k+1)
                          P_l = P_l + P(i-1,j-1,k-1) + P(i+1,j+1,k+1)
                          vx_l = vx_l + vx(i-1,j-1,k-1) + vx(i+1,j+1,k+1)
                          vy_l = vy_l + vy(i-1,j-1,k-1) + vy(i+1,j+1,k+1)
                          vz_l = vz_l + vz(i-1,j-1,k-1) + vz(i+1,j+1,k+1)
                          point_count = point_count + 2
                       end if p4
                       p5: if ((failures(i+1,j-1,k-1) .eqv. .false.) .and. &
                               (failures(i-1,j+1,k+1) .eqv. .false.) .and. &
                               (.not. x_ogb_p) .and. (.not. x_ogb_n) .and. &
                               (.not. y_ogb_p) .and. (.not. y_ogb_n) .and. &
                               (.not. z_ogb_p) .and. (.not. z_ogb_n) ) then
                          rho_bl = rho_bl + rho_b(i+1,j-1,k-1) + rho_b(i-1,j+1,k+1)
                          P_l = P_l + P(i+1,j-1,k-1) + P(i-1,j+1,k+1)
                          vx_l = vx_l + vx(i+1,j-1,k-1) + vx(i-1,j+1,k+1)
                          vy_l = vy_l + vy(i+1,j-1,k-1) + vy(i-1,j+1,k+1)
                          vz_l = vz_l + vz(i+1,j-1,k-1) + vz(i-1,j+1,k+1)
                          point_count = point_count + 2
                       end if p5
                       p6: if ((failures(i-1,j+1,k-1) .eqv. .false.) .and. &
                               (failures(i+1,j-1,k+1) .eqv. .false.).and. &
                               (.not. x_ogb_p) .and. (.not. x_ogb_n) .and. &
                               (.not. y_ogb_p) .and. (.not. y_ogb_n) .and. &
                               (.not. z_ogb_p) .and. (.not. z_ogb_n) ) then
                          rho_bl = rho_bl + rho_b(i-1,j+1,k-1) + rho_b(i+1,j-1,k+1)
                          P_l = P_l + P(i-1,j+1,k-1) + P(i+1,j-1,k+1)
                          vx_l = vx_l + vx(i-1,j+1,k-1) + vx(i+1,j-1,k+1)
                          vy_l = vy_l + vy(i-1,j+1,k-1) + vy(i+1,j-1,k+1)
                          vz_l = vz_l + vz(i-1,j+1,k-1) + vz(i+1,j-1,k+1)
                          point_count = point_count + 2
                       end if p6
                       p7: if ((failures(i-1,j-1,k+1) .eqv. .false.) .and. &
                               (failures(i+1,j+1,k-1) .eqv. .false.) .and. &
                               (.not. x_ogb_p) .and. (.not. x_ogb_n) .and. &
                               (.not. y_ogb_p) .and. (.not. y_ogb_n) .and. &
                               (.not. z_ogb_p) .and. (.not. z_ogb_n) ) then
                          rho_bl = rho_bl + rho_b(i-1,j-1,k+1) + rho_b(i+1,j+1,k-1)
                          P_l = P_l + P(i-1,j-1,k+1) + P(i+1,j+1,k-1)
                          vx_l = vx_l + vx(i-1,j-1,k+1) + vx(i+1,j+1,k-1)
                          vy_l = vy_l + vy(i-1,j-1,k+1) + vy(i+1,j+1,k-1)
                          vz_l = vz_l + vz(i-1,j-1,k+1) + vz(i+1,j+1,k-1)
                          point_count = point_count + 2
                       end if p7
                    end if
                 end if
                 if ((k .ne. kmax) .and. (k .ne. kmin)) then
                    p8: if ((failures(i-1,j,k-1) .eqv. .false.) .and. &
                            (failures(i+1,j,k+1) .eqv. .false.) .and. &
                            (.not. x_ogb_p) .and. (.not. x_ogb_n) .and. &
                            (.not. z_ogb_p) .and. (.not. z_ogb_n) ) then
                       rho_bl = rho_bl + rho_b(i-1,j,k-1) + rho_b(i+1,j,k+1)
                       P_l = P_l + P(i-1,j,k-1) + P(i+1,j,k+1)
                       vx_l = vx_l + vx(i-1,j,k-1) + vx(i+1,j,k+1)
                       vy_l = vy_l + vy(i-1,j,k-1) + vy(i+1,j,k+1)
                       vz_l = vz_l + vz(i-1,j,k-1) + vz(i+1,j,k+1)
                       point_count = point_count + 2
                    end if p8
                    p9: if ((failures(i+1,j,k-1) .eqv. .false.) .and. &
                            (failures(i-1,j,k+1) .eqv. .false.) .and. &
                            (.not. x_ogb_p) .and. (.not. x_ogb_n) .and. &
                            (.not. z_ogb_p) .and. (.not. z_ogb_n) ) then
                       rho_bl = rho_bl + rho_b(i+1,j,k-1) + rho_b(i-1,j,k+1)
                       P_l = P_l + P(i+1,j,k-1) + P(i-1,j,k+1)
                       vx_l = vx_l + vx(i+1,j,k-1) + vx(i-1,j,k+1)
                       vy_l = vy_l + vy(i+1,j,k-1) + vy(i-1,j,k+1)
                       vz_l = vz_l + vz(i+1,j,k-1) + vz(i-1,j,k+1)
                       point_count = point_count + 2
                    end if p9
                 end if
              end if xmax
              ymax: if ((j .ne. jmax) .and. (j .ne. jmin)) then
                 p10: if ((failures(i,j-1,k) .eqv. .false.) .and. &
                          (failures(i,j+1,k) .eqv. .false.) .and. &
                          (.not. y_ogb_p) .and. (.not. y_ogb_n)) then
                    rho_bl = rho_bl + rho_b(i,j-1,k) + rho_b(i,j+1,k)
                    P_l = P_l + P(i,j-1,k) + P(i,j+1,k)
                    vx_l = vx_l + vx(i,j-1,k) + vx(i,j+1,k)
                    vy_l = vy_l + vy(i,j-1,k) + vy(i,j+1,k)
                    vz_l = vz_l + vz(i,j-1,k) + vz(i,j+1,k)
                    point_count = point_count + 2
                 end if p10
                 if ((k .ne. kmax) .and. (k .ne. kmin)) then
                    p11: if ((failures(i,j-1,k-1) .eqv. .false.) .and. &
                             (failures(i,j+1,k+1) .eqv. .false.) .and. &
                             (.not. y_ogb_p) .and. (.not. y_ogb_n) .and. &
                             (.not. z_ogb_p) .and. (.not. z_ogb_n) ) then
                       rho_bl = rho_bl + rho_b(i,j-1,k-1) + rho_b(i,j+1,k+1)
                       P_l = P_l + P(i,j-1,k-1) + P(i,j+1,k+1)
                       vx_l = vx_l + vx(i,j-1,k-1) + vx(i,j+1,k+1)
                       vy_l = vy_l + vy(i,j-1,k-1) + vy(i,j+1,k+1)
                       vz_l = vz_l + vz(i,j-1,k-1) + vz(i,j+1,k+1)
                       point_count = point_count + 2
                    end if p11
                    p12: if ((failures(i,j+1,k-1) .eqv. .false.) .and. &
                             (failures(i,j-1,k+1) .eqv. .false.) .and. &
                             (.not. y_ogb_p) .and. (.not. y_ogb_n) .and. &
                             (.not. z_ogb_p) .and. (.not. z_ogb_n) ) then
                       rho_bl = rho_bl + rho_b(i,j+1,k-1) + rho_b(i,j-1,k+1)
                       P_l = P_l + P(i,j+1,k-1) + P(i,j-1,k+1)
                       vx_l = vx_l + vx(i,j+1,k-1) + vx(i,j-1,k+1)
                       vy_l = vy_l + vy(i,j+1,k-1) + vy(i,j-1,k+1)
                       vz_l = vz_l + vz(i,j+1,k-1) + vz(i,j-1,k+1)
                       point_count = point_count + 2
                    end if p12
                 end if
              end if ymax
              zmax: if ((k .ne. kmax) .and. (k .ne. kmin)) then
                 p13: if ((failures(i,j,k-1) .eqv. .false.) .and. &
                          (failures(i,j,k+1) .eqv. .false.) .and. &
                          (.not. z_ogb_p) .and. (.not. z_ogb_n) ) then
                    rho_bl = rho_bl + rho_b(i,j,k-1) + rho_b(i,j,k+1)
                    P_l = P_l + P(i,j,k-1) + P(i,j,k+1)
                    vx_l = vx_l + vx(i,j,k-1) + vx(i,j,k+1)
                    vy_l = vy_l + vy(i,j,k-1) + vy(i,j,k+1)
                    vz_l = vz_l + vz(i,j,k-1) + vz(i,j,k+1)
                    point_count = point_count + 2
                 end if p13
              end if zmax
              !Deal specially with the eight corners of the grid
              c1: if (i.eq.imax .and. j.eq.jmax .and. k.eq.kmax) then
                 point_count = 0
              end if c1
              c2: if (i.eq.imax .and. j.eq.jmax .and. k.eq.kmin) then
                 point_count = 0
              end if c2
              c3: if (i.eq.imax .and. j.eq.jmin .and. k.eq.kmax) then
                 point_count = 0
              end if c3
              c4: if (i.eq.imax .and. j.eq.jmin .and. k.eq.kmin) then
                 point_count = 0
              end if c4
              c5: if (i.eq.imin .and. j.eq.jmax .and. k.eq.kmax) then
                 point_count = 0
              end if c5
              c6: if (i.eq.imin .and. j.eq.jmax .and. k.eq.kmin) then
                 point_count = 0
              end if c6
              c7: if (i.eq.imin .and. j.eq.jmin .and. k.eq.kmax) then
                 point_count = 0
              end if c7
              c8: if (i.eq.imin .and. j.eq.jmin .and. k.eq.kmin) then
                 point_count = 0
              end if c8
              if (point_count .eq. 0) then
                 rho_bl = rho_b_atm_gf(i,j,k)
                 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                        neos, ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                 P_l = P_cold
                 if (P_l .lt. pfloor_gf(i,j,k)) then
                    P_l = pfloor_gf(i,j,k)
                 end if
                 vx_l = -betax(i,j,k)
                 vy_l = -betay(i,j,k)
                 vz_l = -betaz(i,j,k)
              else
                 !Normalize averages
                 rho_bl = rho_bl / (dble(point_count))
                 P_l = P_l / (dble(point_count))
                 vx_l = vx_l / (dble(point_count))
                 vy_l = vy_l / (dble(point_count))
                 vz_l = vz_l / (dble(point_count))
              end if
              if (rho_bl .lt. rho_b_atm_gf(i,j,k)) then
                 rho_bl = rho_b_atm_gf(i,j,k)
                 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                        neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                 P_l = P_cold
                 if (P_l .lt. pfloor_gf(i,j,k)) then
                    P_l = pfloor_gf(i,j,k)
                 end if
                 vx_l = -betax(i,j,k)
                 vy_l = -betay(i,j,k)
                 vz_l = -betaz(i,j,k)
              end if
              call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                   neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
              Pmax = 50.d0*P_cold
              if (P_l .gt. Pmax) then
                 P_l = Pmax
              end if
              Pmin = pfloor_gf(i,j,k)
              if (P_l .lt. Pmin) then
                 P_l = Pmin
              end if
              ! Set values
              rho_b(i,j,k) = rho_bl
              P(i,j,k) = P_l
              vx(i,j,k) = vx_l
              vy(i,j,k) = vy_l
              vz(i,j,k) = vz_l
              ! Now recalculate conserved variables and other relevant quantitie
              Psi6 = exp(6.d0*phi(i,j,k))
              Psi4 = exp(4.d0*phi(i,j,k))
              call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                        neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
              eps = eps_cold + (P_l-P_cold)/(rho_bl*(gamma_th-ONE))
              h_l = ONE + P_l/rho_bl + eps
              h(i,j,k) = h_l
              ! Compute al*u0-1
              er = Psi4*(gxx(i,j,k)*(vx(i,j,k) + betax(i,j,k))**2 + &
                   2.d0*gxy(i,j,k)*(vx(i,j,k) + betax(i,j,k))*(vy(i,j,k) + betay(i,j,k)) +  &
                   2.d0*gxz(i,j,k)*(vx(i,j,k) + betax(i,j,k))*(vz(i,j,k) + betaz(i,j,k)) +  &
                   gyy(i,j,k)*(vy(i,j,k) + betay(i,j,k))**2 +                               &
                   2.d0*gyz(i,j,k)*(vy(i,j,k) + betay(i,j,k))*(vz(i,j,k) + betaz(i,j,k)) +  &
                   gzz(i,j,k)*(vz(i,j,k) + betaz(i,j,k))**2 )/alpn1**2
              ! *** Check for superluminal velocity ***
              if (er .gt. 1.d0) then
                 vx(i,j,k) = (vx(i,j,k) + betax(i,j,k))*sqrt(fac2/er)-betax(i,j,k)
                 vy(i,j,k) = (vy(i,j,k) + betay(i,j,k))*sqrt(fac2/er)-betay(i,j,k)
                 vz(i,j,k) = (vz(i,j,k) + betaz(i,j,k))*sqrt(fac2/er)-betaz(i,j,k)
                 er = fac2
              end if
              ! ***************************************
              el = sqrt(1.d0-er)
              au0r1 = er/el/(1.d0+el)
              u0(i,j,k)  = (au0r1+1.d0)/alpn1
              rho_s = alpn1 * Psi6 * rho_bl * u0(i,j,k)
              rho_star(i,j,k) = rho_s
              w_l = alpn1 * u0(i,j,k) * rho_s
              w(i,j,k) = w_l
              u_xl = u0(i,j,k)*Psi4*(gxx(i,j,k)*(vx(i,j,k) + betax(i,j,k)) + &
                   gxy(i,j,k)*(vy(i,j,k) + betay(i,j,k)) + &
                   gxz(i,j,k)*(vz(i,j,k) + betaz(i,j,k)))
              u_yl = u0(i,j,k)*Psi4*(gxy(i,j,k)*(vx(i,j,k) + betax(i,j,k)) + &
                   gyy(i,j,k)*(vy(i,j,k) + betay(i,j,k)) + &
                   gyz(i,j,k)*(vz(i,j,k) + betaz(i,j,k)))
              u_zl = u0(i,j,k)*Psi4*(gxz(i,j,k)*(vx(i,j,k) + betax(i,j,k)) + &
                   gyz(i,j,k)*(vy(i,j,k) + betay(i,j,k)) + &
                   gzz(i,j,k)*(vz(i,j,k) + betaz(i,j,k)))
              st_x_l = rho_s*h_l*u_xl
              st_y_l = rho_s*h_l*u_yl
              st_z_l = rho_s*h_l*u_zl
              st_x(i,j,k) = st_x_l
              st_y(i,j,k) = st_y_l
              st_z(i,j,k) = st_z_l
              B_xl  = Psi4 * (gxx(i,j,k)*Bx(i,j,k) + gxy(i,j,k)*By(i,j,k) + gxz(i,j,k)*Bz(i,j,k))
              B_yl  = Psi4 * (gxy(i,j,k)*Bx(i,j,k) + gyy(i,j,k)*By(i,j,k) + gyz(i,j,k)*Bz(i,j,k))
              B_zl  = Psi4 * (gxz(i,j,k)*Bx(i,j,k) + gyz(i,j,k)*By(i,j,k) + gzz(i,j,k)*Bz(i,j,k))
              B2 = Psi4*(gxx(i,j,k)*Bx(i,j,k)*Bx(i,j,k) + TWO*gxy(i,j,k)*Bx(i,j,k)*By(i,j,k) + &
                   TWO*gxz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) + gyy(i,j,k)*By(i,j,k)*By(i,j,k) + &
                   TWO*gyz(i,j,k)*By(i,j,k)*Bz(i,j,k) + gzz(i,j,k)*Bz(i,j,k)*Bz(i,j,k))
              sb0 = (u_xl*Bx(i,j,k) + u_yl*By(i,j,k) + u_zl*Bz(i,j,k))*f1o4pa
              sb2 = ( B2*f1o4pa**2 + sb0**2)/u0(i,j,k)**2
              sb_x = (B_xl*f1o4pa + u_xl*sb0)/u0(i,j,k)
              sb_y = (B_yl*f1o4pa + u_yl*sb0)/u0(i,j,k)
              sb_z = (B_zl*f1o4pa + u_zl*sb0)/u0(i,j,k)
              mhd_st_x(i,j,k) = rho_s*h_l*u_xl + &
                   alpn1*Psi6*u0(i,j,k)*sb2*u_xl - alpn1*Psi6*sb0*sb_x
              mhd_st_y(i,j,k) = rho_s*h_l*u_yl + &
                   alpn1*Psi6*u0(i,j,k)*sb2*u_yl - alpn1*Psi6*sb0*sb_y
              mhd_st_z(i,j,k) = rho_s*h_l*u_zl + &
                   alpn1*Psi6*u0(i,j,k)*sb2*u_zl - alpn1*Psi6*sb0*sb_z
              if (energy_type.eq.1) then
                 u_0n1 = -alpha(i,j,k) - au0r1*alpn1 + &
                      u_xl*betax(i,j,k) + u_yl*betay(i,j,k) + &
                      u_zl*betaz(i,j,k)
                 u_0 = u_0n1 - 1.d0
                 u_0hn1 = u_0n1 + u_0*(P_l/rho_bl+eps)
                 u0u0 = u0(i,j,k)*u_0
                 sb_0 = -sb0*alpn1**2 + sb_x*betax(i,j,k) + &
                      sb_y*betay(i,j,k) + sb_z*betaz(i,j,k)
                 b0b0 = sb0*sb_0
                 tau(i,j,k) = -rho_s*u_0hn1 - alpn1*Psi6* ( &
                      P_l + (u0u0+0.5d0)*sb2 - b0b0 )
              else
                 tau(i,j,k) = (au0r1+(P_l/rho_bl+eps)*alpn1*u0(i,j,k))*rho_s + &
                      Psi6*sb2*(alpn1*u0(i,j,k))**2 &
                      - Psi6*(P_l+sb2*0.5d0)-Psi6*(alpn1*sb0)**2
              end if
              ! Now calculate the sources....
              fac  = ONE / ( Psi6 * w(i,j,k)* h(i,j,k) )
              rho(i,j,k) = h(i,j,k) * w(i,j,k)/Psi6 - P_l
              Sx(i,j,k)  = st_x_l/Psi6
              Sy(i,j,k)  = st_y_l/Psi6
              Sz(i,j,k)  = st_z_l/Psi6
              Sxx(i,j,k) = fac * st_x_l*st_x_l + Psi4 * gxx(i,j,k) * P_l
              Sxy(i,j,k) = fac * st_x_l*st_y_l + Psi4 * gxy(i,j,k) * P_l
              Sxz(i,j,k) = fac * st_x_l*st_z_l + Psi4 * gxz(i,j,k) * P_l
              Syy(i,j,k) = fac * st_y_l*st_y_l + Psi4 * gyy(i,j,k) * P_l
              Syz(i,j,k) = fac * st_y_l*st_z_l + Psi4 * gyz(i,j,k) * P_l
              Szz(i,j,k) = fac * st_z_l*st_z_l + Psi4 * gzz(i,j,k) * P_l
              !
              ! MHD metric sources
              E_x = Psi6/alpn1 * ( By(i,j,k)*(vz(i,j,k)+betaz(i,j,k)) - Bz(i,j,k)*(vy(i,j,k)+betay(i,j,k)) )
              E_y = Psi6/alpn1 * ( Bz(i,j,k)*(vx(i,j,k)+betax(i,j,k)) - Bx(i,j,k)*(vz(i,j,k)+betaz(i,j,k)) )
              E_z = Psi6/alpn1 * ( Bx(i,j,k)*(vy(i,j,k)+betay(i,j,k)) - By(i,j,k)*(vx(i,j,k)+betax(i,j,k)) )
              Ex = (gupxx(i,j,k)*E_x + gupxy(i,j,k)*E_y + gupxz(i,j,k)*E_z)/Psi4
              Ey = (gupxy(i,j,k)*E_x + gupyy(i,j,k)*E_y + gupyz(i,j,k)*E_z)/Psi4
              Ez = (gupxz(i,j,k)*E_x + gupyz(i,j,k)*E_y + gupzz(i,j,k)*E_z)/Psi4
              temp3 = f1o8p*(Ex*E_x + Ey*E_y + Ez*E_z + Bx(i,j,k)*B_xl + By(i,j,k)*B_yl + Bz(i,j,k)*B_zl)
              rho(i,j,k) = rho(i,j,k) + temp3
              Sxx(i,j,k) = Sxx(i,j,k) + temp3*Psi4*gxx(i,j,k) - f1o4p*(E_x*E_x + B_xl*B_xl)
              Sxy(i,j,k) = Sxy(i,j,k) + temp3*Psi4*gxy(i,j,k) - f1o4p*(E_x*E_y + B_xl*B_yl)
              Sxz(i,j,k) = Sxz(i,j,k) + temp3*Psi4*gxz(i,j,k) - f1o4p*(E_x*E_z + B_xl*B_zl)
              Syy(i,j,k) = Syy(i,j,k) + temp3*Psi4*gyy(i,j,k) - f1o4p*(E_y*E_y + B_yl*B_yl)
              Syz(i,j,k) = Syz(i,j,k) + temp3*Psi4*gyz(i,j,k) - f1o4p*(E_y*E_z + B_yl*B_zl)
              Szz(i,j,k) = Szz(i,j,k) + temp3*Psi4*gzz(i,j,k) - f1o4p*(E_z*E_z + B_zl*B_zl)
              Sx(i,j,k)  = Sx(i,j,k)  + f1o4p*Psi6*(Ey*Bz(i,j,k) - Ez*By(i,j,k))
              Sy(i,j,k)  = Sy(i,j,k)  + f1o4p*Psi6*(Ez*Bx(i,j,k) - Ex*Bz(i,j,k))
              Sz(i,j,k)  = Sz(i,j,k)  + f1o4p*Psi6*(Ex*By(i,j,k) - Ey*Bx(i,j,k))
              sbt(i,j,k) = u_xl*Bx(i,j,k) + u_yl*By(i,j,k) + u_zl*Bz(i,j,k)
              sbx(i,j,k) = Bx(i,j,k)/u0(i,j,k) + vx(i,j,k)*sbt(i,j,k)
              sby(i,j,k) = By(i,j,k)/u0(i,j,k) + vy(i,j,k)*sbt(i,j,k)
              sbz(i,j,k) = Bz(i,j,k)/u0(i,j,k) + vz(i,j,k)*sbt(i,j,k)
!!$              write(*,*) repaired point: i, j, k = , i,j,k
!!$              write(*,*) rho = , rho(i,j,k)
!!$              write(*,*) u0 = , u0(i,j,k)
!!$              write(*,*) vx = , vx(i,j,k)
!!$              write(*,*) vy = , vy(i,j,k)
!!$              write(*,*) vz = , vz(i,j,k)
!!$              write(*,*) Bx = , Bx(i,j,k)
!!$              write(*,*) By = , By(i,j,k)
!!$              write(*,*) Bz = , Bz(i,j,k)
!!$              write(*,*) Ex = , Ex
!!$              write(*,*) Ey = , Ey
!!$              write(*,*) Ez = , Ez
!!$              write(*,*) h = , h(i,j,k)
!!$              write(*,*) w = , w(i,j,k)
!!$              write(*,*) Psi6 = , Psi6
!!$              write(*,*) P_l = , P_l
!!$              write(*,*) rho_s = , rho_s
!!$              write(*,*) rho_bl = , rho_bl
           end if outer
        end do
     end do
  end do
end subroutine repair_failures_mhd_alt
