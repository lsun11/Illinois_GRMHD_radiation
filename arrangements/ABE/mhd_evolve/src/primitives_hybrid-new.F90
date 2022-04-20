!-----------------------------------------------------------------------------
!
! reconstruct primitive variables, compute sources for "hybrid" EOS
! (version for running density floor)
!
!-----------------------------------------------------------------------------
subroutine primitive_vars_hybrid2(ext, cctk_nghostzones,X, Y, Z, rho_star, tau, &
     st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
     neos, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
     w, w_old, rho_b, rho, P, h, Sx, Sy, Sz, &
     Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     phi, alpha, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
     h_old, u0, rho_max, rho_b_atm, rho_fail_max_step, M_fail_step, rhos_max, &
     Bx, By, Bz, Ex, Ey, Ez, vx, vy, vz, sbt, sbx, sby, sbz, &
     proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
     glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor,excision_enable,excision_zone_gf, &
     tau_stildefix_enable)
  implicit none
  ! Hydro Input
  integer, dimension(3),intent(in)                    :: ext,cctk_nghostzones
  integer					      :: neos
  real*8, dimension(neos)			      :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1)			      :: k_tab, gamma_tab
  real*8					      :: gamma_th,rho_b_atm,pfloor
  real*8, dimension(ext(1),ext(2),ext(3)), intent(in) :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3))             :: h_old
  real*8, dimension(ext(1),ext(2),ext(3))             :: rho_star,tau
  real*8, dimension(ext(1),ext(2),ext(3))             :: st_x,st_y,st_z
  real*8, dimension(ext(1),ext(2),ext(3))             :: mhd_st_x,mhd_st_y,mhd_st_z
  real*8                                              :: rho_max,rhos_max
  integer                                             :: proc_imin,proc_jmin,proc_kmin
  integer                                             :: proc_imax,proc_jmax,proc_kmax
  integer                                             :: glob_imax,glob_jmax,glob_kmax
  integer                                             :: excision_enable,tau_stildefix_enable
  integer, dimension(ext(1),ext(2),ext(3))            :: excision_zone_gf


  ! Output
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: w, rho_b, rho
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: P, h
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: w_old
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Sx, Sy, Sz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Sxx, Sxy, Sxz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Syy, Syz, Szz
  real*8, dimension(ext(1),ext(2),ext(3))             :: vx,vy,vz
  ! Metric
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: phi,alpha
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: shiftx,shifty,shiftz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: gupxx, gupxy, gupxz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)  :: gupyy, gupyz, gupzz
  ! MHD
  real*8, dimension(ext(1),ext(2),ext(3))	      :: Bx,By,Bz,Ex,Ey,Ez
  real*8, dimension(ext(1),ext(2),ext(3))	      :: sbt,sbx,sby,sbz
  ! Internal variables:
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  real*8                             :: P_l, eps, eps_star,u2,Psi6,Psi4,Psi2
  real*8                             :: rho_s, st_x_l, st_y_l,st_z_l, Tiny
  real*8                             :: deltE,B_x_l,B_y_l,B_z_l
  real*8                             :: E_x_l,E_y_l,E_z_l,alpn1,Fx_l,Fy_l,Fz_l
  real*8                             :: ux_l,uy_l,uz_l, u_xl,u_yl,u_zl
  real*8                             :: Pmax,Pmin,rho_bl,h_l,P_tiny,BB
  real*8                             :: gijuiuj,au0m1,sb0,sb2,sb_x,sb_y,sb_z
  real*8                             :: fac, B2, dX, dY, dZ, dV
  real*8                             :: rho_fail_max_step, M_fail_step
  real*8                             :: P_cold,eps_cold,eps_tiny
  real*8, parameter                  :: TWO  = 2.D0
  real*8, parameter                  :: ZERO = 0.D0
  real*8, parameter                  :: ONE  = 1.D0
  real*8, parameter                  :: FOUR = 4.D0
  real*8, parameter                  :: SIX  = 6.D0
  real*8                             :: u_scal, eps_scal, sti_scal, tau_scal
  real*8                             :: f, rho_cutoff, factor
  real*8, parameter		     :: max_gamma = 60.d0
  real*8, dimension(ext(1),ext(2),ext(3)) :: E_x, E_y, E_z, B_x, B_y, B_z
  real*8, dimension(ext(1),ext(2),ext(3)) :: u_x, u_y, u_z, u0
  real*8, dimension(ext(1),ext(2),ext(3)) :: temp, psi_6, psi_4
  real*8 :: PI,f1o4p, f1o8p,f1o4pa
  real*8, dimension(4) :: UU
  real*8, dimension(3) :: UU_font_fix
  real*8, dimension(27) :: AUX
  integer :: Symmetry
  integer, parameter :: AXISYM = 4
  integer :: count,nn,m
  logical :: check, recom,exit_do, repairs_needed
  logical, dimension(ext(1),ext(2),ext(3)) :: failures
  external funcv_hybrid,fdjac_hybrid, funcv_hybrid_font_fix,fdjac_hybrid_font_fix
  external funcv_hybrid_font_fix2,fdjac_hybrid_font_fix2
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

  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
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
  rho_fail_max_step = 0.d0
  M_fail_step =0.d0
  !  Font = 0.d0
  !  Nfont = 0
  !-----------------------------------------------------------------------------
  ! Funny parameters...  (See Shibata, Oohara & Nakamura...)
  !-----------------------------------------------------------------------------
  f = 4.D-5
  rho_cutoff = 1D-4
  !  Tiny = 1.D-10*rho_max
  Tiny = 0.d0

  psi_4 = exp(4.d0*phi)
  psi_6 = exp(6.d0*phi)

  u_x = (gxx*(shiftx+vx) + gxy*(shifty+vy) + gxz*(shiftz+vz))*psi_4*u0
  u_y = (gxy*(shiftx+vx) + gyy*(shifty+vy) + gyz*(shiftz+vz))*psi_4*u0
  u_z = (gxz*(shiftx+vx) + gyz*(shifty+vy) + gzz*(shiftz+vz))*psi_4*u0

  B_x  = psi_4 * (gxx * Bx + gxy * By + gxz * Bz)
  B_y  = psi_4 * (gxy * Bx + gyy * By + gyz * Bz)
  B_z  = psi_4 * (gxz * Bx + gyz * By + gzz * Bz)




  if (Symmetry == AXISYM) then 
     u0 = ONE
  end if

  ! Initialize the failures array
  failures = .FALSE.
  repairs_needed = .FALSE.


  ! Set eps_scal
  call compute_pcold_epscold(rho_max, P_cold, eps_scal, &
       neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
  !
  ! Now go to each gridpoint
  !
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax

           ! Set eps_tiny and P_tiny
           call compute_pcold_epscold(rho_b_atm, P_tiny, eps_tiny, &
                neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab) 

   	   Psi4 = psi_4(i,j,k)
   	   Psi6 = psi_6(i,j,k)
           Psi2 = sqrt(Psi4)
	   alpn1 = alpha(i,j,k) + ONE

!!$           rho_s = max(alpn1*rho_star(i,j,k),  & 
!!$		abs(alpn1)*rho_b_atm*Psi6) / alpn1
!!$	   rho_star(i,j,k) = rho_s
           if (rho_b_atm*Psi6 .lt. rho_star(i,j,k)) then
              rho_s = rho_star(i,j,k)
           else
              rho_s = -1.d0
           end if
	   f1o4pa = sqrt(f1o4p)/alpn1

!!$           if (abs(rho_s) .gt. Tiny) then
           if (rho_s .gt. Tiny) then
              if(tau_stildefix_enable == 1) then
                 if(tau(i,j,k).lt.zero) then
                    if(k.gt.1)write(6,*)'fix3:', &
                         i,j,k,X(i,j,k),Y(i,j,k),Z(i,j,k),tau(i,j,k)
                    tau(i,j,k)=1.0e-3
                 endif
                 stxi=mhd_st_x(i,j,k)
                 styi=mhd_st_y(i,j,k)
                 stzi=mhd_st_z(i,j,k)
                 gxxi=gupxx(i,j,k)/Psi4
                 gyyi=gupyy(i,j,k)/Psi4
                 gzzi=gupzz(i,j,k)/Psi4
                 gxyi=gupxy(i,j,k)/Psi4
                 gxzi=gupxz(i,j,k)/Psi4
                 gyzi=gupyz(i,j,k)/Psi4
                 sdots=gxxi*stxi**2+gyyi*styi**2+gzzi*stzi**2+2.0* &
                      (gxyi*stxi*styi+gxzi*stxi*stzi+gyzi*styi*stzi)
                 rhot=tau(i,j,k)*(tau(i,j,k)+2.0*rho_s)
                 if(sdots.gt.0.98*rhot) then
                    if(k.gt.1)write(6,*)'fix2:', &
                         i,j,k,X(i,j,k),Y(i,j,k),Z(i,j,k), &
                         sdots,rhot,tau(i,j,k),rho_s

                    rfact=sdots/(0.98*rhot)
                    mhd_st_x(i,j,k)=mhd_st_x(i,j,k)/rfact
                    mhd_st_y(i,j,k)=mhd_st_y(i,j,k)/rfact
                    mhd_st_z(i,j,k)=mhd_st_z(i,j,k)/rfact
                 end if
              end if
 	      rho_bl = rho_b(i,j,k)
              eps = h_old(i,j,k)-1.d0-P(i,j,k)/rho_bl
              nn = 4
              UU(1) = u_x(i,j,k)/Psi2
	      UU(2) = u_y(i,j,k)/Psi2
	      UU(3) = u_z(i,j,k)/Psi2
	      UU(4) = max(eps,eps_tiny)
              m = 27
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
              AUX(10) = alpn1
              AUX(11) = B2
              AUX(12) = B_x(i,j,k)*f1o4pa
   	      AUX(13) = B_y(i,j,k)*f1o4pa
	      AUX(14) = B_z(i,j,k)*f1o4pa
              AUX(15) = gupxx(i,j,k)/Psi4
              AUX(16) = gupxy(i,j,k)/Psi4
              AUX(17) = gupxz(i,j,k)/Psi4
              AUX(18) = gupyy(i,j,k)/Psi4
              AUX(19) = gupyz(i,j,k)/Psi4
              AUX(20) = gupzz(i,j,k)/Psi4
              AUX(21) = Psi6
	      AUX(22) = Psi2
              !	      P_l = max(P(i,j,k),P_tiny)
	      P_l = max(P(i,j,k),pfloor)
              h_l = max(h_old(i,j,k),1.d0)
              u_scal  = sqrt( (P_l + B2)/(rho_bl*h_l+B2) )
	      bb = Psi4*( gxx(i,j,k)*shiftx(i,j,k)**2 + & 
                   2.d0*gxy(i,j,k)*shiftx(i,j,k)*shifty(i,j,k) + &
                   2.d0*gxz(i,j,k)*shiftx(i,j,k)*shiftz(i,j,k) + &
                   gyy(i,j,k)*shifty(i,j,k)**2 + &
                   2.d0*gyz(i,j,k)*shifty(i,j,k)*shiftz(i,j,k) + &
                   gzz(i,j,k)*shiftz(i,j,k)**2 )
	      u_scal = u_scal + sqrt(bb)
              sti_scal = max(abs(mhd_st_x(i,j,k)), abs(mhd_st_y(i,j,k)), &
                   abs(mhd_st_z(i,j,k)), &
                   0.001d0*u_scal*Psi2*Psi6*(rho_bl*h_old(i,j,k)+B2) )
              tau_scal = max(abs(tau(i,j,k)), 0.01d0*Psi6*(rho_bl* &
                   eps_scal+0.5d0*B2) )
              UU(1) = UU(1)/u_scal
              UU(2) = UU(2)/u_scal
              UU(3) = UU(3)/u_scal
              UU(4) = UU(4)/eps_scal
	      AUX(23) = u_scal
	      AUX(24) = eps_scal
	      AUX(25) = sti_scal
	      AUX(26) = tau_scal
              AUX(27) = pfloor

              call newt2(UU,AUX,funcv_hybrid,fdjac_hybrid,nn,m, & 
                   neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,check)

	      recom = .FALSE.

              ! Impose Font fix when the solver fails (check = .TRUE.) or when it
              !  gives negative eps [UU(4) < 0]
              !
              if(check .or. UU(4) .lt. 0.d0) then
                 !Font(i,j,k) = 1.0
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

                 count = count + 1
                 nn = 3
                 check = .FALSE.
                 call newt2(UU_font_fix, &
                      AUX,funcv_hybrid_font_fix2, & 
                      fdjac_hybrid_font_fix2,nn,m, & 
                      neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,check) 
                 if(check) then 
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

                    failures(i,j,k) = .TRUE.
                    repairs_needed = .TRUE.
                    ! Set everything to zero before calculating 
                    ! these quantities from averages.
                    u_xl         = ZERO
                    u_yl         = ZERO
                    u_zl         = ZERO
                    rho_bl       = ZERO
                    P_l          = ZERO
                    eps          = ZERO
                    u0(i,j,k)    = 1.d0/alpn1
                    recom        = .TRUE.
                 else  !if the Font fix worked, do the following:                 
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
                    if (rho_s .lt. 0.d0) u0(i,j,k) = -u0(i,j,k)
                    rho_bl = rho_s/alpn1/u0(i,j,k)/AUX(21)
		    call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                         neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)

		    P_l = P_cold 
		    eps = eps_cold

                 end if
              else
                 u_xl = UU(1)*Psi2*u_scal
                 u_yl = UU(2)*Psi2*u_scal
                 u_zl = UU(3)*Psi2*u_scal
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
                 if (rho_s .lt. 0.d0) u0(i,j,k) = -u0(i,j,k)
                 rho_bl = rho_s/alpn1/u0(i,j,k)/AUX(21)
		 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                      neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                 P_l = P_cold + rho_bl*(gamma_th-1.d0)*(eps-eps_cold)
              end if
	      u_x(i,j,k) = u_xl
	      u_y(i,j,k) = u_yl
              u_z(i,j,k) = u_zl
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
              !testing *********!!!!!!!!!!!!
              ! the following lines are only for the disk accretion runs, 
              ! which are performed without excision.
!!$	      if (rho_bl .gt. 1.d3) then
!!$		 recom = .TRUE.
!!$		 rho_bl = 1.d3
!!$	       	 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
!!$                        neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
!!$		 P_l = P_cold 
!!$		 eps = eps_cold
!!$                 h_l = 1.d0 + P_l/rho_bl + eps
!!$	      end if
              ! end testing lines

	      if (rho_bl .lt. rho_b_atm) then
		 recom = .TRUE.
		 rho_bl = rho_b_atm
	       	 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                      neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
		 P_l = P_cold 
		 eps = eps_cold
                 h_l = 1.d0 + P_l/rho_bl + eps
	      end if
              !              Pmax = 10.d0*P_cold
              Pmax = 50.d0*P_cold
              if (P_l .gt. Pmax) then
		 recom = .TRUE.
                 P_l = Pmax
                 eps = eps_cold + (Pmax-P_cold)/(gamma_th-1.d0)/rho_bl
                 h_l = 1.d0 + eps + Pmax/rho_bl
	      end if
              !              Pmin = pfloor(i,j,k)
              Pmin = pfloor
              !              Pmin = 0.5d0*P_cold
              if (P_l .lt. Pmin) then
		 recom = .TRUE.
                 P_l = Pmin
                 eps = eps_cold + (Pmin-P_l)/(gamma_th-1.d0)/rho_bl
                 h_l = 1.d0 + eps + Pmin/rho_bl
	      end if
              ! Re-compute rho_star, mhd_st_i and tau if necessary
	      if (recom) then
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
		 tau(i,j,k) = (au0m1+(P_l/rho_bl+eps)*alpn1*u0(i,j,k))*rho_s + &
                      AUX(21)*sb2*(alpn1*u0(i,j,k))**2 &
                      - AUX(21)*(P_l+sb2*0.5d0)-AUX(21)*(alpn1*sb0)**2
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
!!$	      rho_star(i,j,k) = ZERO
!!$              tau(i,j,k)   = ZERO
!!$	      Bx(i,j,k)    = ZERO
!!$              By(i,j,k)    = ZERO
!!$              Bz(i,j,k)    = ZERO
!!$	      st_x(i,j,k)  = ZERO
!!$	      st_y(i,j,k)  = ZERO
!!$	      st_z(i,j,k)  = ZERO
!!$	      vx(i,j,k)    = ZERO
!!$              vy(i,j,k)    = ZERO
!!$              vz(i,j,k)    = ZERO
!!$	      u_x(i,j,k)   = ZERO
!!$              u_y(i,j,k)   = ZERO
!!$              u_z(i,j,k)   = ZERO
!!$	      w(i,j,k)     = ZERO
!!$              h(i,j,k)     = ONE
!!$              rho_b(i,j,k) = ZERO
!!$              P(i,j,k)     = ZERO
!!$              rho(i,j,k)   = ZERO
!!$              Sx(i,j,k)    = ZERO
!!$              Sy(i,j,k)    = ZERO
!!$              Sz(i,j,k)    = ZERO
!!$              Sxx(i,j,k)   = ZERO
!!$              Sxy(i,j,k)   = ZERO
!!$              Sxz(i,j,k)   = ZERO
!!$              Syy(i,j,k)   = ZERO
!!$              Syz(i,j,k)   = ZERO
!!$              Szz(i,j,k)   = ZERO
              ! Set rho_b to the atmosphere density and u_i=0
              rho_bl = rho_b_atm
              call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                   neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
              P_l = P_cold
              eps = eps_cold
              ! testing!
              if (P_l .lt. pfloor) then
                 eps = eps + (pfloor-P_l)/(rho_bl*(gamma_th - 1.d0))
                 P_l = pfloor
              end if

              h_l = 1.d0 + P_l/rho_bl + eps
              u0(i,j,k) = 1.d0/alpn1
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
              tau(i,j,k)   = Psi6*(rho_bl*eps + sb2*0.5d0)
	      st_x(i,j,k)  = ZERO
	      st_y(i,j,k)  = ZERO
	      st_z(i,j,k)  = ZERO
	      mhd_st_x(i,j,k) = ZERO
              mhd_st_y(i,j,k) = ZERO
              mhd_st_z(i,j,k) = ZERO
	      u_x(i,j,k)   = ZERO
              u_y(i,j,k)   = ZERO
              u_z(i,j,k)   = ZERO
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
        end do
     end do
  end do
  !
  ! MHD metric sources
  E_x = psi_6 * ( By*(vz+shiftz) - Bz*(vy+shifty) )/(1.d0+alpha)
  E_y = psi_6 * ( Bz*(vx+shiftx) - Bx*(vz+shiftz) )/(1.d0+alpha)
  E_z = psi_6 * ( Bx*(vy+shifty) - By*(vx+shiftx) )/(1.d0+alpha)
  Ex = (gupxx*E_x + gupxy*E_y + gupxz*E_z)/psi_4
  Ey = (gupxy*E_x + gupyy*E_y + gupyz*E_z)/psi_4
  Ez = (gupxz*E_x + gupyz*E_y + gupzz*E_z)/psi_4
  B_x  = psi_4 * (gxx * Bx + gxy * By + gxz * Bz)
  B_y  = psi_4 * (gxy * Bx + gyy * By + gyz * Bz)
  B_z  = psi_4 * (gxz * Bx + gyz * By + gzz * Bz)
  temp = f1o8p*(Ex*E_x + Ey*E_y + Ez*E_z + Bx*B_x + By*B_y + Bz*B_z)

  rho   = rho + temp
  Sxx   = Sxx + temp*psi_4*gxx - f1o4p*(E_x*E_x + B_x*B_x)
  Sxy   = Sxy + temp*psi_4*gxy - f1o4p*(E_x*E_y + B_x*B_y)
  Sxz   = Sxz + temp*psi_4*gxz - f1o4p*(E_x*E_z + B_x*B_z)
  Syy   = Syy + temp*psi_4*gyy - f1o4p*(E_y*E_y + B_y*B_y)
  Syz   = Syz + temp*psi_4*gyz - f1o4p*(E_y*E_z + B_y*B_z)
  Szz   = Szz + temp*psi_4*gzz - f1o4p*(E_z*E_z + B_z*B_z)
  Sx    = Sx + f1o4p*psi_6*(Ey*Bz - Ez*By)
  Sy    = Sy + f1o4p*psi_6*(Ez*Bx - Ex*Bz)
  Sz    = Sz + f1o4p*psi_6*(Ex*By - Ey*Bx)

  sbt = u_x*Bx + u_y*By + u_z*Bz
  sbx = Bx/u0 + vx*sbt
  sby = By/u0 + vy*sbt
  sbz = Bz/u0 + vz*sbt
  if (Symmetry==AXISYM) then
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
     end do
  end if

  M_fail_step = M_fail_step * dV
  !  Nfont = count
  write(*,*) 'Fixed ',count,' zones'

  if(repairs_needed) then
     call repair_failures_mhd_hybrid2(ext,Z,gamma_th, failures, rho_b, P, & 
          vx, vy, vz, u0, w, h, rho_star, tau, &
          st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, & 
          rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
          alpha, shiftx, shifty, shiftz, phi, &
          gxx, gxy, gxz, gyy, gyz, gzz, &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
          sbt, sbx, sby, sbz, rho_b_atm, & 
          neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
          proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
          glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor)
  end if

  if(excision_enable==1) then
     call hydro_ezbc_hybrid(ext,X,Y,Z,rho_star,tau, &
          mhd_st_x,mhd_st_y,mhd_st_z,st_x,st_y,st_z, &
          rho_b,P,h,vx,vy,vz,w,&
          sbt,sbx,sby,sbz,Bx,By,Bz, &
          alpha,shiftx,shifty,shiftz,phi,& 
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,& 
          Symmetry,excision_zone_gf,gamma_th,neos,rho_tab,&
          P_tab,eps_tab,k_tab,gamma_tab);
     call remove_interior2(ext,X,Y,Z,sbt,excision_zone_gf,Symmetry);
     call remove_interior2(ext,X,Y,Z,sbx,excision_zone_gf,Symmetry);
     call remove_interior2(ext,X,Y,Z,sby,excision_zone_gf,Symmetry);
     call remove_interior2(ext,X,Y,Z,sbz,excision_zone_gf,Symmetry);
  end if


end subroutine primitive_vars_hybrid2


subroutine funcv_hybrid(n,x,fvec,m,aux,neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
  implicit none
  integer                    :: n,m
  real*8, dimension(n)       :: x,fvec
  real*8, dimension(m)       :: aux
  real*8 :: rho_s,tau,h,gamma_th
  real*8 :: Bx,By,Bz,mhd_stx,mhd_sty,mhd_stz
  real*8 :: B_x,B_y,B_z, B2
  real*8 :: u_x,u_y,u_z,u0,eps,P,rhob,P_cold,eps_cold
  real*8 :: alp,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,sqrtg
  real*8 :: sb0,sb_x,sb_y,sb_z,sb2,gijuiuj,au0m1,Psi2
  real*8 :: u_scal2,eps_scal,sti_scal,tau_scal
  integer :: neos
  logical :: exit_do
  real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab, gamma_tab
!
  Psi2 = aux(22)
  u_scal2 = aux(23)*Psi2
  eps_scal = aux(24)
  sti_scal = aux(25)
  tau_scal = aux(26)
  u_x    = x(1)*u_scal2
  u_y    = x(2)*u_scal2
  u_z    = x(3)*u_scal2
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
  alp   = aux(10)
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
  gijuiuj = gupxx*u_x**2 + 2.d0*gupxy*u_x*u_y + &
                2.d0*gupxz*u_x*u_z + gupyy*u_y**2 + 2.d0*gupyz*u_y*u_z + &
                gupzz*u_z**2
  au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
  if (rho_s .lt. 0.d0) au0m1 = gijuiuj/( 1.d0-sqrt(1.d0+gijuiuj) )
  u0 = (au0m1+1.d0)/alp
  rhob = abs(rho_s/alp/sqrtg/u0)
  call compute_pcold_epscold(rhob, P_cold, eps_cold, &
                neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
  P = P_cold + (gamma_th - 1.d0)*rhob*(eps - eps_cold)

  if (rhob .ne. 0.d0) then
     h   = 1.d0 + P/rhob + eps
  else
     h = 1.d0
  end if
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
  if (rhob .ne. 0.d0) then
     fvec(4) = (au0m1+alp*u0*(P/rhob+eps))*rho_s + sqrtg*sb2*(alp*u0)**2 &
	- sqrtg*(P+sb2*0.5d0)-sqrtg*(alp*sb0)**2 - tau
  else
     fvec(4) = sqrtg*sb2*(alp*u0)**2 &
        - sqrtg*(P+sb2*0.5d0)-sqrtg*(alp*sb0)**2 - tau
  end if
  fvec(4) = fvec(4)/tau_scal
end subroutine funcv_hybrid

subroutine fdjac_hybrid(n,x,m,aux,neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
			fvec,NP,fjac)
  implicit none
  integer                :: n,m,NP,neos,i
  real*8, dimension(n)   :: x,fvec
  real*8, dimension(m)   :: aux
  real*8, dimension(NP,NP) :: fjac
  real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab,gamma_tab
  real*8 :: rho_s,h,gamma_th
  real*8 :: Bx,By,Bz,mhd_stx,mhd_sty,mhd_stz
  real*8 :: B_x,B_y,B_z,B2,u_x,u_y,u_z,eps
  real*8 :: alp,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,sqrtg
  real*8 :: u0,sb0,sb2,sb_x,sb_y,sb_z,tau,P,rhob
  real*8 :: dPcold_drho, depscold_drho,P_cold,eps_cold
  real*8 :: du0dux,du0duy,du0duz, db2dux,db2duy,db2duz
  real*8 :: dbxdux,dbxduy,dbxduz, dbydux,dbyduy,dbyduz, dbzdux,dbzduy,dbzduz
  real*8 :: dpdux,dpduy,dpduz,gijuiuj,au0m1, c, Psi2
  real*8 :: drho_dux, drho_duy, drho_duz, dP_drho, dhdux,dhduy,dhduz
  real*8 :: u_scal2,eps_scal,sti_scal,tau_scal
  logical :: exit_do
!
  Psi2   = aux(22)
  u_scal2 = aux(23)*Psi2
  eps_scal = aux(24)
  sti_scal = aux(25)
  tau_scal = aux(26)
  u_x    = x(1)*u_scal2
  u_y    = x(2)*u_scal2
  u_z    = x(3)*u_scal2
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
  alp  = aux(10)
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

  gijuiuj = gupxx*u_x**2 + 2.d0*gupxy*u_x*u_y + &
                2.d0*gupxz*u_x*u_z + gupyy*u_y**2 + 2.d0*gupyz*u_y*u_z + &
                gupzz*u_z**2
  au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
  if (rho_s .lt. 0.d0) au0m1 = gijuiuj/( 1.d0-sqrt(1.d0+gijuiuj) )
  u0 = (au0m1+1.d0)/alp
  rhob = abs(rho_s/alp/sqrtg/u0)

  i = 1
  exit_do = .FALSE.
  do
     if (rhob .lt. rho_tab(i)) then
        exit_do = .TRUE.
        P_cold = k_tab(i)*rhob**gamma_tab(i)
	dPcold_drho = gamma_tab(i)*P_cold/rhob
	depscold_drho = P_cold/rhob**2
        if (i==1) then
           eps_cold = P_cold/rhob/(gamma_tab(i)-1.d0)
        else
           eps_cold = eps_tab(i-1) +  &
                (P_cold/rhob - P_tab(i-1)/rho_tab(i-1))/(gamma_tab(i)-1.d0)
        end if
     end if
     if (i==neos .or. exit_do) exit
     i = i + 1
  end do
  if (rhob .gt. rho_tab(neos)) then
     P_cold = k_tab(neos+1)*rhob**gamma_tab(neos+1)
     eps_cold = eps_tab(neos) + (P_cold/rhob - P_tab(neos)/rho_tab(neos))
     dPcold_drho = gamma_tab(neos+1)*P_cold/rhob
     depscold_drho = P_cold/rhob**2 
  end if
  P = P_cold + (gamma_th - 1.d0)*(eps - eps_cold)*rhob
  dP_drho = dPcold_drho + (gamma_th-1.d0)*(eps-eps_cold-rhob*depscold_drho)

  if (rhob .ne. 0.d0) then
     h   = 1.d0 + P/rhob + eps
  else
     h   = 1.d0
  end if
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
	 - (Bx*sb_x + sb0*dbxdux) ) + rho_s*u_x*dhdux ) * u_scal2/sti_scal
  fjac(1,2) = ( alp*sqrtg*( u_x*(u0*db2duy+sb2*du0duy)  & 
	 - (By*sb_x + sb0*dbxduy) ) + rho_s*u_x*dhduy ) * u_scal2/sti_scal
  fjac(1,3) = ( alp*sqrtg*( u_x*(u0*db2duz+sb2*du0duz)  &
	 - (Bz*sb_x + sb0*dbxduz) ) + rho_s*u_x*dhduz ) * u_scal2/sti_scal
  fjac(1,4) = gamma_th*rho_s*u_x * eps_scal/sti_scal

  fjac(2,1) = ( alp*sqrtg*( u_y*(u0*db2dux+sb2*du0dux)  &
         - (Bx*sb_y + sb0*dbydux) ) + rho_s*u_y*dhdux ) * u_scal2/sti_scal
  fjac(2,2) = ( c + alp*sqrtg*( u_y*(u0*db2duy+sb2*du0duy)  &
         - (By*sb_y + sb0*dbyduy) ) + rho_s*u_y*dhduy ) * u_scal2/sti_scal
  fjac(2,3) = ( alp*sqrtg*( u_y*(u0*db2duz+sb2*du0duz)  &
         - (Bz*sb_y + sb0*dbyduz) ) + rho_s*u_y*dhduz ) * u_scal2/sti_scal
  fjac(2,4) = gamma_th*rho_s*u_y * eps_scal/sti_scal

  fjac(3,1) = ( alp*sqrtg*( u_z*(u0*db2dux+sb2*du0dux)  &
         - (Bx*sb_z + sb0*dbzdux) ) + rho_s*u_z*dhdux ) * u_scal2/sti_scal
  fjac(3,2) = ( alp*sqrtg*( u_z*(u0*db2duy+sb2*du0duy)  &
         - (By*sb_z + sb0*dbzduy) ) + rho_s*u_z*dhduy ) * u_scal2/sti_scal
  fjac(3,3) = ( c + alp*sqrtg*( u_z*(u0*db2duz+sb2*du0duz)  &
         - (Bz*sb_z + sb0*dbzduz) ) + rho_s*u_z*dhduz ) * u_scal2/sti_scal
  fjac(3,4) = gamma_th*rho_s*u_z * eps_scal/sti_scal

  fjac(4,1) = ( alp*rho_s*h*du0dux + alp*alp*sqrtg*(2.d0*u0*sb2*du0dux + &
		db2dux*u0**2 - 2.d0*sb0*Bx) - sqrtg*(dpdux + 0.5d0*db2dux) + &
		alp*rho_s*u0*dhdux ) * u_scal2/tau_scal
  fjac(4,2) = ( alp*rho_s*h*du0duy + alp*alp*sqrtg*(2.d0*u0*sb2*du0duy + &
                db2duy*u0**2 - 2.d0*sb0*By) - sqrtg*(dpduy + 0.5d0*db2duy) + &
		alp*rho_s*u0*dhduy ) * u_scal2/tau_scal
  fjac(4,3) = ( alp*rho_s*h*du0duz + alp*alp*sqrtg*(2.d0*u0*sb2*du0duz + &
                db2duz*u0**2 - 2.d0*sb0*Bz) - sqrtg*(dpduz + 0.5d0*db2duz) + &
		alp*rho_s*u0*dhduz ) * u_scal2/tau_scal
  fjac(4,4) = ( alp*u0*gamma_th*rho_s - (gamma_th-1.d0)*rho_s/(alp*u0) ) * & 
                eps_scal/tau_scal

end subroutine fdjac_hybrid

subroutine funcv_hybrid_font_fix(n,x,fvec,m,aux,neos,rho_tab, & 
			P_tab,eps_tab,k_tab,gamma_tab)
  implicit none
  integer                    :: n,m
  real*8, dimension(n)       :: x,fvec
  real*8, dimension(m)       :: aux
  real*8 :: rho_s,tau,h,gamma_th
  real*8 :: Bx,By,Bz,mhd_stx,mhd_sty,mhd_stz
  real*8 :: B_x,B_y,B_z, B2
  real*8 :: u_x,u_y,u_z,u0,eps,P,rhob
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
  alp   = aux(10)
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
  gijuiuj = gupxx*u_x**2 + 2.d0*gupxy*u_x*u_y + &
                2.d0*gupxz*u_x*u_z + gupyy*u_y**2 + 2.d0*gupyz*u_y*u_z + &
                gupzz*u_z**2
  au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
  if (rho_s .lt. 0.d0) au0m1 = gijuiuj/( 1.d0-sqrt(1.d0+gijuiuj) )
  u0 = (au0m1+1.d0)/alp
  rhob = rho_s/alp/sqrtg/u0
  call compute_pcold_epscold(rhob, P, eps, &
                neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)

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
end subroutine funcv_hybrid_font_fix

subroutine funcv_hybrid_font_fix2(n,x,fvec,m,aux,neos,rho_tab, & 
			P_tab,eps_tab,k_tab,gamma_tab)
  implicit none
  integer                    :: n,m
  real*8, dimension(n)       :: x,fvec
  real*8, dimension(m)       :: aux
  real*8 :: rho_s,tau,h,gamma_th
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
  alp   = aux(10)
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
  pfloor = aux(27)
  gijuiuj = gupxx*u_x**2 + 2.d0*gupxy*u_x*u_y + &
                2.d0*gupxz*u_x*u_z + gupyy*u_y**2 + 2.d0*gupyz*u_y*u_z + &
                gupzz*u_z**2
  au0m1 = gijuiuj/( 1.d0+sqrt(1.d0+gijuiuj) )
  if (rho_s .lt. 0.d0) au0m1 = gijuiuj/( 1.d0-sqrt(1.d0+gijuiuj) )
  u0 = (au0m1+1.d0)/alp
  rhob = rho_s/alp/sqrtg/u0
  call compute_pcold_epscold(rhob, P, eps, &
                neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)

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

subroutine fdjac_hybrid_font_fix(n,x,m,aux,neos,rho_tab,P_tab, & 
		eps_tab,k_tab,gamma_tab,fvec,NP,fjac)
  implicit none
  integer                :: n,m,NP,neos,i
  real*8, dimension(n)   :: x,fvec
  real*8, dimension(m)   :: aux
  real*8, dimension(NP,NP) :: fjac
  real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab,gamma_tab
  real*8, dimension(3) :: fac
  real*8 :: rho_s,h,gamma_th
  real*8 :: Bx,By,Bz,mhd_stx,mhd_sty,mhd_stz
  real*8 :: B_x,B_y,B_z,B2,u_x,u_y,u_z,eps
  real*8 :: alp,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,sqrtg
  real*8 :: u0,sb0,sb2,sb_x,sb_y,sb_z,tau,P,rhob
  real*8 :: du0dux,du0duy,du0duz, db2dux,db2duy,db2duz
  real*8 :: dbxdux,dbxduy,dbxduz, dbydux,dbyduy,dbyduz, dbzdux,dbzduy,dbzduz
  real*8 :: dpdux,dpduy,dpduz,gijuiuj,au0m1, c
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
  alp  = aux(10)
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
        if (i==1) then
           eps = P/rhob/(gamma_tab(i)-1.d0)
        else
           eps = eps_tab(i) +  &
                (P/rhob - P_tab(i-1)/rho_tab(i-1))/(gamma_tab(i)-1.d0)
        end if
     end if
     if (i==neos .or. exit_do) exit
     i = i + 1
  end do
  if (rhob .gt. rho_tab(neos)) then
     P = k_tab(neos+1)*rhob**gamma_tab(neos+1)
     eps = eps_tab(neos) + (P/rhob - P_tab(neos)/rho_tab(neos))
     dP_drho = gamma_tab(neos+1)*P/rhob
     deps_drho = P/rhob**2
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

end subroutine fdjac_hybrid_font_fix

subroutine fdjac_hybrid_font_fix2(n,x,m,aux,neos,rho_tab,P_tab, & 
		eps_tab,k_tab,gamma_tab,fvec,NP,fjac)
  implicit none
  integer                :: n,m,NP,neos,i
  real*8, dimension(n)   :: x,fvec
  real*8, dimension(m)   :: aux
  real*8, dimension(NP,NP) :: fjac
  real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab,gamma_tab
  real*8, dimension(3) :: fac
  real*8 :: rho_s,h,gamma_th
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
  alp  = aux(10)
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
  pfloor  = aux(27)

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
        if (i==1) then
           eps = P/rhob/(gamma_tab(i)-1.d0)
        else
           eps = eps_tab(i) +  &
                (P/rhob - P_tab(i-1)/rho_tab(i-1))/(gamma_tab(i)-1.d0)
        end if
     end if
     if (i==neos .or. exit_do) exit
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

      FUNCTION fmin2(x,aux,fvec,funcv1,n,m,neos,rho_tab,P_tab,eps_tab, & 
			k_tab,gamma_tab)
      implicit none
      INTEGER  :: n,m,NP
      real*8   :: fmin2
      real*8, dimension(n) :: x,fvec
      real*8, dimension(m) :: aux
      integer neos
      real*8 rho_tab(neos), P_tab(neos), eps_tab(neos)
      real*8 k_tab(neos+1), gamma_tab(neos+1)
      external funcv1
      PARAMETER (NP=40)
!      COMMON /newtv/ fvec(NP),n,m
!      SAVE /newtv/
! CU    USES funcv
      INTEGER i
      real*8 sum
      call funcv1(n,x,fvec,m,aux,neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
      sum=0.d0
      do 11 i=1,n
        sum=sum+fvec(i)**2
11    continue
      fmin2=0.5d0*sum
      return
      END

      SUBROUTINE lubksb2(a,n,np,indx,b)
      implicit none
      INTEGER n,np,indx(n)
      real*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      real*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END

      SUBROUTINE ludcmp2(a,n,np,indx,d)
      implicit none
      INTEGER n,np,indx(n),NMAX
      real*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      real*8 aamax,dum,sum,vv(NMAX)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) then
           write(*,*) 'singular matrix in ludcmp, i = ', i
	   stop
        end if
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END

      SUBROUTINE lnsrch2(n,xold,fold,g,p,x,f,stpmax,check,func,funcv1,m,aux,fvec, &
                neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
      implicit none
      INTEGER n,m
      LOGICAL check
      real*8 f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX,aux(m)
      real*8 fvec(n)
      integer neos
      real*8 rho_tab(neos), P_tab(neos), eps_tab(neos)
      real*8 k_tab(neos+1), gamma_tab(neos+1)
      PARAMETER (ALF=1.d-4,TOLX=1.d-15)
      EXTERNAL func, funcv1
! CU    USES func
      INTEGER i
      real*8 a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp, &
           test,tmplam
      check=.false.
      sum=0.d0
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.d0
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.d0
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.d0)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.d0
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=func(x,aux,fvec,funcv1,n,m,neos,rho_tab,P_tab,eps_tab, &
                        k_tab,gamma_tab)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.d0)then
            tmplam=-slope/(2.d0*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.d0)then
              tmplam=-slope/(2.d0*b)
            else
              disc=b*b-3.d0*a*slope
              if(disc.lt.0.d0) write(*,*) 'roundoff problem in lnsrch ', &
		disc,b,a,slope
              tmplam=(-b+sqrt(disc))/(3.d0*a)
            endif
            if(tmplam.gt.0.5d0*alam)tmplam=0.5d0*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,0.1d0*alam)
      goto 1
      END

      SUBROUTINE newt2(x,aux,funcv1,fdjac1,n,m,neos,rho_tab, & 
                P_tab,eps_tab,k_tab,gamma_tab,check)
      implicit none
      INTEGER n,m,nn,mm,NP,MAXITS
      LOGICAL check
      real*8 x(n),aux(m),fvec(n),TOLF,TOLMIN,TOLX,STPMX
      real*8 pos
      PARAMETER (NP=40,MAXITS=1000,TOLF=1.d-13,TOLMIN=1.d-15,TOLX=1.d-13, & 
                 STPMX=0.5d0)
!      COMMON /newtv/ fvec(NP),nn,mm
!      SAVE /newtv/
! CU    USES fdjac,fmin,lnsrch,lubksb,ludcmp
      INTEGER i,its,j,indx(NP)
      real*8 d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP), &
           xold(NP),fmin2
      integer neos
      real*8 rho_tab(neos), P_tab(neos), eps_tab(neos) 
      real*8 k_tab(neos+1), gamma_tab(neos+1)
      EXTERNAL fmin2
      external funcv1,fdjac1
!
      nn=n
      mm=m
      f=fmin2(x,aux,fvec,funcv1,n,m,neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
      test=0.d0
      do 11 i=1,n
        test = max( test, abs(fvec(i)) )
11    continue
      if(test.lt.0.01d0*TOLF)then
        check=.false.
        return
      endif
      sum=0.d0
      do 12 i=1,n
        sum=sum+x(i)**2
12    continue
      stpmax=STPMX*max(sqrt(sum),dble(n))
      do 21 its=1,MAXITS
        call fdjac1(n,x,m,aux,neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, & 
			fvec,NP,fjac)
        do 14 i=1,n
          sum=0.d0
          do 13 j=1,n
            sum=sum+fjac(j,i)*fvec(j)
13        continue
          g(i)=sum
14      continue
        do 15 i=1,n
          xold(i)=x(i)
15      continue
        fold=f
        do 16 i=1,n
          p(i)=-fvec(i)
16      continue
        call ludcmp2(fjac,n,NP,indx,d)
        call lubksb2(fjac,n,NP,indx,p)
        call lnsrch2(n,xold,fold,g,p,x,f,stpmax,check,fmin2,funcv1,m,aux,fvec, & 
		neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
        test=0.d0
        do 17 i=1,n
          test = max( test, abs(fvec(i)) )
17      continue
        if(test.lt.TOLF)then
          check=.false.
          return
        endif
        if(check)then
          test=0.d0
          den=max(f,0.5d0*n)
          do 18 i=1,n
            temp=abs(g(i))*max(abs(x(i)),1.d0)/den
            if(temp.gt.test)test=temp
18        continue
          if(test.lt.TOLMIN)then
            check=.true.
          else
            check=.false.
          endif
          return
        endif
        test=0.d0
        do 19 i=1,n
          temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.d0)
          test = max(test,temp)
19      continue
        if(test.lt.TOLX)return
21    continue
!      write(*,*) 'MAXITS exceeded in newt'
      check=.TRUE.
      END

subroutine repair_failures_mhd_hybrid(ext,Z, gamma_th, failures, rho_b, P, & 
     vx, vy, vz, u0, w, h, rho_star, tau, &
     st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
     rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     alpha, betax, betay, betaz, phi, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
     sbt, sbx, sby, sbz, rho_b_atm, & 
     neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
     proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
     glob_imax,glob_jmax,glob_kmax,Symmetry)
  implicit none
!
! Input
!
  integer, dimension(3),intent(in)                :: ext
  real*8, dimension(ext(3))                       :: Z
  real*8                                          :: n, rho_b_atm, gamma_th
  logical, dimension(ext(1),ext(2),ext(3))        :: failures
  real*8, dimension(ext(1),ext(2),ext(3))	  :: rho_b, P, vx, vy, vz
  real*8, dimension(ext(1),ext(2),ext(3))	  :: u0, w, h, rho_star, tau
  real*8, dimension(ext(1),ext(2),ext(3))         :: st_x, st_y, st_z
  real*8, dimension(ext(1),ext(2),ext(3))         :: mhd_st_x, mhd_st_y, mhd_st_z
  real*8, dimension(ext(1),ext(2),ext(3))         :: Bx, By, Bz
  real*8, dimension(ext(1),ext(2),ext(3))         :: sbt, sbx, sby, sbz
  integer 					  :: neos
  real*8, dimension(neos) 			  :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) 			  :: k_tab, gamma_tab
  integer                                         :: proc_imin,proc_jmin,proc_kmin
  integer                                         :: proc_imax,proc_jmax,proc_kmax
  integer                                         :: glob_imax,glob_jmax,glob_kmax
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
  real*8                                :: w_l, tau_l, fac
  real*8                                :: alpha_l, temp3
  real*8 		                :: rho_bl, h_l, eps, P_cold, eps_cold
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
  logical				:: exit_do
  integer :: Symmetry
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  logical :: x_ogb_p,x_ogb_n,y_ogb_p,y_ogb_n,z_ogb_p,z_ogb_n
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
!
! Go to each gridpoint
!

  do k=kmin, kmax
     do j=jmin, jmax
        do i=imin, imax

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
                        ((SYMMETRY .eq. NO_SYMM) .or. (Z(kmin) .lt. 0.d0))

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

              c1: if (i==imax .and. j==jmax .and. k==kmax) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c1

              c2: if (i==imax .and. j==jmax .and. k==kmin) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c2

              c3: if (i==imax .and. j==jmin .and. k==kmax) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c3

              c4: if (i==imax .and. j==jmin .and. k==kmin) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c4

              c5: if (i==imin .and. j==jmax .and. k==kmax) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c5

              c6: if (i==imin .and. j==jmax .and. k==kmin) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c6

              c7: if (i==imin .and. j==jmin .and. k==kmax) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c7

              c8: if (i==imin .and. j==jmin .and. k==kmin) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c8

              if (point_count == 0) then
                 write(*,*) 'No healthy pairs surrounding point i,j,k = ', i,j,k
                 rho_bl = rho_b_atm
		 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                        neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                 P_l = P_cold
                 vx_l = 0.d0
                 vy_l = 0.d0
                 vz_l = 0.d0
              else
                 !Normalize averages
                 rho_bl = rho_bl / (dble(point_count))
                 P_l = P_l / (dble(point_count))
                 vx_l = vx_l / (dble(point_count))
                 vy_l = vy_l / (dble(point_count))
                 vz_l = vz_l / (dble(point_count))
              end if

              if (rho_bl .lt. rho_b_atm) then
                 rho_bl = rho_b_atm
		 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                        neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                 P_l = P_cold
                 vx_l = 0.d0
                 vy_l = 0.d0
                 vz_l = 0.d0
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
                        neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
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
              tau(i,j,k) = (au0r1 + & 
		   (P_l/rho_bl+eps)*alpn1*u0(i,j,k))*rho_s + &
                   Psi6*sb2*(alpn1*u0(i,j,k))**2 &
                   - Psi6*(P_l+sb2*0.5d0)-Psi6*(alpn1*sb0)**2
              
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

              write(*,*) 'repaired point: i, j, k = ', i,j,k
              write(*,*) 'rho = ', rho(i,j,k)
!!$              write(*,*) 'u0 = ', u0(i,j,k)
!!$              write(*,*) 'vx = ', vx(i,j,k)
!!$              write(*,*) 'vy = ', vy(i,j,k)
!!$              write(*,*) 'vz = ', vz(i,j,k)
!!$              write(*,*) 'Bx = ', Bx(i,j,k)
!!$              write(*,*) 'By = ', By(i,j,k)
!!$              write(*,*) 'Bz = ', Bz(i,j,k)
!!$              write(*,*) 'Ex = ', Ex
!!$              write(*,*) 'Ey = ', Ey
!!$              write(*,*) 'Ez = ', Ez
              write(*,*) 'h = ', h(i,j,k)
              write(*,*) 'w = ', w(i,j,k)
              write(*,*) 'Psi6 = ', Psi6
              write(*,*) 'P_l = ', P_l
              write(*,*) 'rho_s = ', rho_s
              write(*,*) 'rho_bl = ', rho_bl
           end if outer
           
        end do
     end do
  end do

end subroutine repair_failures_mhd_hybrid

subroutine compute_pcold_epscold(rhob, P_cold, eps_cold, &
		neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
  implicit none
  integer :: neos,i
  logical :: exit_do
  real*8  :: rhob, P_cold, eps_cold
  real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab, gamma_tab
!
  i = 1
  exit_do = .FALSE.
  do
     if (rhob .lt. rho_tab(i)) then
        exit_do = .TRUE.
        P_cold = k_tab(i)*rhob**gamma_tab(i)
        if (i==1) then
	   if (rhob .ne. 0.d0) then
              eps_cold = P_cold/rhob/(gamma_tab(i)-1.d0)
	   else
	      eps_cold = 0.d0
	   end if
        else
           eps_cold = eps_tab(i-1) +  &
                (P_cold/rhob - P_tab(i-1)/rho_tab(i-1))/(gamma_tab(i)-1.d0)
        end if
     end if
     if (i==neos .or. exit_do) exit
     i = i + 1
  end do
  if (rhob .gt. rho_tab(neos)) then
     P_cold = k_tab(neos+1)*rhob**gamma_tab(neos+1)
     eps_cold = eps_tab(neos) +  & 
            (P_cold/rhob - P_tab(neos)/rho_tab(neos))/(gamma_tab(neos+1)-1.d0)
  end if
end subroutine compute_pcold_epscold
subroutine repair_failures_mhd_hybrid2(ext,Z, gamma_th, failures, rho_b, P, & 
     vx, vy, vz, u0, w, h, rho_star, tau, &
     st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
     rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     alpha, betax, betay, betaz, phi, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
     sbt, sbx, sby, sbz, rho_b_atm, & 
     neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
     proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
     glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor)
  implicit none
!
! Input
!
  integer, dimension(3),intent(in)                :: ext
  real*8, dimension(ext(3))                       :: Z
  real*8                                          :: n, gamma_th,rho_b_atm, pfloor
  logical, dimension(ext(1),ext(2),ext(3))        :: failures
  real*8, dimension(ext(1),ext(2),ext(3))	  :: rho_b, P, vx, vy, vz
  real*8, dimension(ext(1),ext(2),ext(3))	  :: u0, w, h, rho_star, tau
  real*8, dimension(ext(1),ext(2),ext(3))         :: st_x, st_y, st_z
  real*8, dimension(ext(1),ext(2),ext(3))         :: mhd_st_x, mhd_st_y, mhd_st_z
  real*8, dimension(ext(1),ext(2),ext(3))         :: Bx, By, Bz
  real*8, dimension(ext(1),ext(2),ext(3))         :: sbt, sbx, sby, sbz
  integer 					  :: neos
  real*8, dimension(neos) 			  :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) 			  :: k_tab, gamma_tab
  integer                                         :: proc_imin,proc_jmin,proc_kmin
  integer                                         :: proc_imax,proc_jmax,proc_kmax
  integer                                         :: glob_imax,glob_jmax,glob_kmax
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
  real*8                                :: w_l, tau_l, fac
  real*8                                :: alpha_l, temp3
  real*8 		                :: rho_bl, h_l, eps, P_cold, eps_cold
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
  logical				:: exit_do
  integer :: Symmetry
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  logical :: x_ogb_p,x_ogb_n,y_ogb_p,y_ogb_n,z_ogb_p,z_ogb_n
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
!
! Go to each gridpoint
!

  do k=kmin, kmax
     do j=jmin, jmax
        do i=imin, imax

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
                        ((SYMMETRY .eq. NO_SYMM) .or. (Z(kmin) .lt. 0.d0))

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

              c1: if (i==imax .and. j==jmax .and. k==kmax) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c1

              c2: if (i==imax .and. j==jmax .and. k==kmin) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c2

              c3: if (i==imax .and. j==jmin .and. k==kmax) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c3

              c4: if (i==imax .and. j==jmin .and. k==kmin) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c4

              c5: if (i==imin .and. j==jmax .and. k==kmax) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c5

              c6: if (i==imin .and. j==jmax .and. k==kmin) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c6

              c7: if (i==imin .and. j==jmin .and. k==kmax) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c7

              c8: if (i==imin .and. j==jmin .and. k==kmin) then
                 write(*,*) 'oops: failed in corner ', i,j,k
                 point_count = 0
              end if c8

              if (point_count == 0) then
                 write(*,*) 'No healthy pairs surrounding point i,j,k = ', i,j,k
                 rho_bl = rho_b_atm
		 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                        neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                 P_l = P_cold
                 vx_l = 0.d0
                 vy_l = 0.d0
                 vz_l = 0.d0
              else
                 !Normalize averages
                 rho_bl = rho_bl / (dble(point_count))
                 P_l = P_l / (dble(point_count))
                 vx_l = vx_l / (dble(point_count))
                 vy_l = vy_l / (dble(point_count))
                 vz_l = vz_l / (dble(point_count))
              end if

              if (rho_bl .lt. rho_b_atm) then
                 rho_bl = rho_b_atm
		 call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
                        neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                 P_l = P_cold

                 if (P_l .lt. pfloor) then
                    P_l = pfloor
                 end if
                 vx_l = 0.d0
                 vy_l = 0.d0
                 vz_l = 0.d0
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
                        neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
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
              tau(i,j,k) = (au0r1 + & 
		   (P_l/rho_bl+eps)*alpn1*u0(i,j,k))*rho_s + &
                   Psi6*sb2*(alpn1*u0(i,j,k))**2 &
                   - Psi6*(P_l+sb2*0.5d0)-Psi6*(alpn1*sb0)**2
              
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

              write(*,*) 'repaired point: i, j, k = ', i,j,k
              write(*,*) 'rho = ', rho(i,j,k)
!!$              write(*,*) 'u0 = ', u0(i,j,k)
!!$              write(*,*) 'vx = ', vx(i,j,k)
!!$              write(*,*) 'vy = ', vy(i,j,k)
!!$              write(*,*) 'vz = ', vz(i,j,k)
!!$              write(*,*) 'Bx = ', Bx(i,j,k)
!!$              write(*,*) 'By = ', By(i,j,k)
!!$              write(*,*) 'Bz = ', Bz(i,j,k)
!!$              write(*,*) 'Ex = ', Ex
!!$              write(*,*) 'Ey = ', Ey
!!$              write(*,*) 'Ez = ', Ez
              write(*,*) 'h = ', h(i,j,k)
              write(*,*) 'w = ', w(i,j,k)
              write(*,*) 'Psi6 = ', Psi6
              write(*,*) 'P_l = ', P_l
              write(*,*) 'rho_s = ', rho_s
              write(*,*) 'rho_bl = ', rho_bl
           end if outer
           
        end do
     end do
  end do

end subroutine repair_failures_mhd_hybrid2
