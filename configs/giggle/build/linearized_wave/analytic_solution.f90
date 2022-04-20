! Linearized gravitational wave analytic solution. 
! This is used for setting up initial data and comparing numerical
!     vs. analytic results.
subroutine lin_wave_analytic(ex, X, Y, Z, T, Amp, Width, &
     gxx, gxy, gxz, gyy, gyz, gzz, kxx,kxy,kxz,kyy,kyz,kzz,PhysR, dRdr,mode)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: PhysR, dRdr
  real*8, dimension(ex(1),ex(2),ex(3)),target :: gxx,gxy,gxz
  real*8, dimension(ex(1),ex(2),ex(3)),target :: gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3)),target :: kxx,kxy,kxz
  real*8, dimension(ex(1),ex(2),ex(3)),target :: kyy,kyz,kzz
  real*8                                      :: T, Amp, Width
  integer                                      :: mode
  integer                            :: i,j,k,l,m,n,o
  integer                            :: imin,jmin,kmin,imax,jmax,kmax
  real*8                             :: R,Rho,SinTh,CosTh,SinPh,CosPh,thetaangle,phiangle
  real*8                             :: A,B,C, KK,LL
  real*8                             :: Ad,Bd,Cd, KKd,LLd
  real*8, dimension(3,3)             :: trans,g_pol, k_pol
  real*8                             :: frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp
  real*8                              :: drt,drp,dtt,dtp,dpp
  real*8                             :: ZERO, ONE, TWO, THREE
  parameter ( ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0 )
!
! go to each gridpoint...
!
  do k = 1,ex(3)
     do j = 1,ex(2)
        do i = 1,ex(1)
!
! find polar coordinates from cartesian coordinates
!
           R     = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           Rho   = sqrt(X(i,j,k)**2 + Y(i,j,k)**2)
           SinPh = Y(i,j,k) / Rho
           CosPh = X(i,j,k) / Rho
           SinTh = Rho / R
           CosTh = Z(i,j,k) / R
!!$           !FIXME: UGLY KLUDGE
!!$           if(R==0) R=1.D-10
!!$
!!$           ! Heres a new formulation for the trans() matrix and Sin/Cos Th/Ph, so that
!!$           !   we dont get NaNs from Rho=0 points
!!$           thetaangle=atan2(Rho,z(i,j,k))
!!$           phiangle=atan2(y(i,j,k),x(i,j,k))
!!$
!!$           SinTh = sin(thetaangle)
!!$           CosTh = cos(thetaangle)
!!$           SinPh = sin(phiangle)
!!$           CosPh = cos(phiangle)
!!$           !
!!$           ! find transformation matrix trans_{ij} = dx^{i}/dx^j, where
!!$           ! primes are polar coordinates
!!$           !
!!$           trans(1,1) = SinTh*CosPh
!!$           trans(1,2) = SinTh*SinPh
!!$           trans(1,3) = CosTh
!!$           trans(2,1) = CosTh*CosPh/R
!!$           trans(2,2) = CosTh*SinPh/R
!!$           trans(2,3) = - SinTh/R
!!$           trans(3,1) = - SinPh/(R*SinTh)
!!$           trans(3,2) = CosPh/(R*SinTh)
!!$           trans(3,3) = ZERO
           trans(1,1) = X(i,j,k)/R
           trans(1,2) = Y(i,j,k)/R
           trans(1,3) = Z(i,j,k)/R
           trans(2,1) = X(i,j,k)*Z(1,1,k)/(R**2 * Rho)
           trans(2,2) = Y(i,j,k)*Z(1,1,k)/(R**2 * Rho)
           trans(2,3) = - Rho/R**2
           trans(3,1) = - Y(i,j,k)/Rho**2
           trans(3,2) = X(i,j,k)/Rho**2
           trans(3,3) = ZERO
!
! find coefficients A, B, C (eq. (6) ) and K and L (eq. (9))
!
           call coefficients(PhysR(i,j,k),T,Amp,Width,A,B,C,KK,LL)
            call dt_coefficients(PhysR(i,j,k),T,Amp,Width,Ad,Bd,Cd,KKd,LLd)
!
! find coefficients f and d (eqs. (7)  and (10))
! 
           call fij_dij(frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,drt,drp,dtt,dtp,dpp, &
                        SinTh,CosTh,SinPh,CosPh,mode)
!
! Find metric in polar coordinates (eq. (5) )
!
           g_pol(1,1) = (ONE + A*frr) * dRdr(i,j,k)**2
           g_pol(1,2) = (B*frt+KK*drt)*PhysR(i,j,k) * dRdr(i,j,k)
           g_pol(1,3) = (B*frp+KK*drp)*PhysR(i,j,k)*SinTh * dRdr(i,j,k)
           g_pol(2,2) = (ONE + C*f1tt + A*f2tt + LL*dtt) * PhysR(i,j,k)**2
           g_pol(2,3) = ( (A-2.d0*C)*ftp + LL*dtp)*SinTh*PhysR(i,j,k)**2
           g_pol(3,3) = (ONE + C*f1pp + A*f2pp + LL*dpp)*(PhysR(i,j,k)*SinTh)**2
           g_pol(2,1) = g_pol(1,2)
           g_pol(3,1) = g_pol(1,3)
           g_pol(3,2) = g_pol(2,3)
           k_pol(1,1) = -0.5d0*Ad*frr * dRdr(i,j,k)**2
           k_pol(1,2) = -0.5d0*(Bd*frt+KKd*drt)*PhysR(i,j,k) * dRdr(i,j,k)
           k_pol(1,3) = -0.5d0*(Bd*frp+KKd*drp)*PhysR(i,j,k)*SinTh * dRdr(i,j,k)
           k_pol(2,2) = -0.5d0*(Cd*f1tt + Ad*f2tt + LLd*dtt) * PhysR(i,j,k)**2
           k_pol(2,3) = -0.5d0*((Ad-2.d0*Cd)*ftp+LLd*dtp)*SinTh*PhysR(i,j,k)**2
           k_pol(3,3) = -0.5d0*(Cd*f1pp + Ad*f2pp + LLd*dpp)*(PhysR(i,j,k)*SinTh)**2
           k_pol(2,1) = k_pol(1,2)
           k_pol(3,1) = k_pol(1,3)
           k_pol(3,2) = k_pol(2,3)
!
! Transform metric into cartesian coordinates
!
           gxx(i,j,k) = ZERO
             kxx(i,j,k) = ZERO
           do n = 1,3
              do o = 1,3
                 gxx(i,j,k) = gxx(i,j,k) + &
                      trans(n,1)*trans(o,1)*g_pol(n,o)
                 kxx(i,j,k) = kxx(i,j,k) + &
                      trans(n,1)*trans(o,1)*k_pol(n,o)
              end do
           end do
           gxy(i,j,k) = ZERO
           kxy(i,j,k) = ZERO
           do n = 1,3
              do o = 1,3
                 gxy(i,j,k) = gxy(i,j,k) + &
                      trans(n,1)*trans(o,2)*g_pol(n,o)
                 kxy(i,j,k) = kxy(i,j,k) + &
                      trans(n,1)*trans(o,2)*k_pol(n,o)
              end do
           end do
           gxz(i,j,k) = ZERO
           kxz(i,j,k) = ZERO
           do n = 1,3
              do o = 1,3
                 gxz(i,j,k) = gxz(i,j,k) + &
                      trans(n,1)*trans(o,3)*g_pol(n,o)
                 kxz(i,j,k) = kxz(i,j,k) + &
                      trans(n,1)*trans(o,3)*k_pol(n,o)
              end do
           end do
           gyy(i,j,k) = ZERO
           kyy(i,j,k) = ZERO
           do n = 1,3
              do o = 1,3
                 gyy(i,j,k) = gyy(i,j,k) + &
                      trans(n,2)*trans(o,2)*g_pol(n,o)
                 kyy(i,j,k) = kyy(i,j,k) + &
                      trans(n,2)*trans(o,2)*k_pol(n,o)
              end do
           end do
           gyz(i,j,k) = ZERO
           kyz(i,j,k) = ZERO
           do n = 1,3
              do o = 1,3
                 gyz(i,j,k) = gyz(i,j,k) + &
                      trans(n,2)*trans(o,3)*g_pol(n,o)
                 kyz(i,j,k) = kyz(i,j,k) + &
                      trans(n,2)*trans(o,3)*k_pol(n,o)
              end do
           end do
           gzz(i,j,k) = ZERO
           kzz(i,j,k) = ZERO
           do n = 1,3
              do o = 1,3
                 gzz(i,j,k) = gzz(i,j,k) + &
                      trans(n,3)*trans(o,3)*g_pol(n,o)
                 kzz(i,j,k) = kzz(i,j,k) + &
                      trans(n,3)*trans(o,3)*k_pol(n,o)
            end do
         end do
        end do
     end do
  end do
end subroutine lin_wave_analytic
! Compute the analytic GW. 
!
subroutine gw_anal(T,PhysR,th,ph,hplus,hcross,E_GW,psi4r,psi4i,psi4_22, &
                   psi4_21,psi4_20,Amp,Width,mode)
  implicit none
  integer :: mode
  real*8 :: T,PhysR,th,ph,hplus,hcross
  real*8 :: SinTh,CosTh,SinPh,CosPh,E_GW
  real*8 :: frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp
  real*8 :: drt,drp,dtt,dtp,dpp
  real*8 :: Amp,Width,A,B,C,K,L,psi4r,psi4i,fac,Cdd,Ldd
  real*8 :: psi4_22,psi4_21,psi4_20
  real*8, parameter :: pi = 3.14159265358979323846d0
!
  SinTh = sin(th)
  CosTh = cos(th)
  SinPh = sin(ph)
  CosPh = cos(ph)
  !
  ! find coefficients A, B and C (eq. (6) )
  !
  call coefficients(PhysR,T,Amp,Width,A,B,C,K,L)
  !
  ! find coefficients f or d (eqs. (7) and (10) )
  !
  f1tt = 0.d0
  ftp = 0.d0
  dtt = 0.d0
  dtp = 0.d0
  call fij_dij(frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,drt,drp,dtt,dtp,dpp, &
                        SinTh,CosTh,SinPh,CosPh,mode)
  call compute_E_GW(T,PhysR,E_GW,Amp,Width,mode)
  if (mode==20 .or. mode==22) then
     call compute_Cdd(T,PhysR,Amp,Width,Cdd)
     Ldd = 0.d0
  else
     call compute_Ldd(T,PhysR,Amp,Width,Ldd)
     Cdd = 0.d0
  end if
!  Rlmr = 0.d0
!  Rlmi = 0.d0
!  call dij_21(drt,drp,dtt,dtp,dpp,SinTh,CosTh,SinPh,CosPh)
!  call E_GW_21(T,PhysR,E_GW,Amp,Width)
!  call Q_21(T,PhysR,Rlmr,Rlmi,Amp,Width)
! *** test ***
!! even mode
!  fac = sqrt(2.d0*acos(-1.d0)/15.d0)
!  write(*,*) H_222 = ,2.d0*A*fac
!  write(*,*) h_122 = ,PhysR*B*fac
!  write(*,*) G_22 = ,-0.5d0*(A-2.d0*C)*fac
!  write(*,*) K_22 = ,-A*fac
!  call R_22(T,PhysR,Rlmr,Rlmi,Amp,Width)
!  write(*,*) R_22 = ,Rlmr
!! odd mode
!  fac = sqrt(8.d0*acos(-1.d0)/15.d0)
!  write(*,*) r,t = ,PhysR,T
!  write(*,*) L, K = ,L,K
!  write(*,*) Im(C_21) = ,PhysR*K*fac
!  write(*,*) Im(D_21) = ,0.25d0*L*fac
!  write(*,*) Im(Q_21) = ,Rlmi
!  write(*,*) hplus = ,PhysR*L*dtt
!  write(*,*) hcross = ,PhysR*L*dtp
!  write(*,*) dE/dt = ,E_GW
! ************
  hplus = PhysR*(C*f1tt + L*dtt)
  hcross = PhysR*(L*dtp - 2.d0*C*ftp)
  psi4r = (Cdd*f1tt + Ldd*dtt)
  psi4i = -(Ldd*dtp - 2.d0*Cdd*ftp)
  ! Decompose psi4 into lm modes
  psi4_22 = 0.d0
  psi4_21 = 0.d0
  psi4_20 = 0.d0
  if (mode==22) then
     psi4_22 = 4.d0*Cdd*sqrt(pi/5.d0)
  else if (mode==21) then
     psi4_21 = -2.d0*Ldd*sqrt(pi/5.d0)
  else
     psi4_20 = 4.d0*Cdd*sqrt(6.d0*pi/5.d0)
  end if
!  write(*,*) cdd =,Cdd
  ! for 22 mode ONLY
!  psi4r = Cdd*sqrt(6.D0*3.1415926535D0/5.D0)
!  psi4i = Cdd
!  write(*,*) psi4r = ,psi4r
!  psi4i = 0.D0
end subroutine gw_anal
subroutine coefficients(r,t,amp,lam,A,B,C,K,L)
  implicit none
  real*8,intent(in)           :: r,t,lam,amp
  real*8,intent(out)          :: A,B,C,K,L
! Function f of t+r and t-r and its derivatives
  real*8                      :: fm,fm1,fm2,fm3,fm4,fm5
  real*8                      :: fp,fp1,fp2,fp3,fp4,fp5
  real*8                      :: tmr,tpr,epp, epm, exa
  real*8                      :: F1o4
  parameter ( F1o4 = 0.25D0 )
!
! Construct function f and its first four derivatives
!
  tmr = t - r
  exa = - lam * tmr**2
  if (exa .lt. -691.d0) then
     epm = 0.d0
  else
     epm = exp(exa)
  end if
  fm  = tmr * epm
  fm1 = ( 1.D0 - 2.D0*lam*tmr**2 ) * epm
  fm2 = (- 6.D0*lam*tmr + 4.D0*lam**2*tmr**3 ) * epm
  fm3 = (- 6.D0*lam + 24.D0*lam**2*tmr**2 - 8.D0*lam**3*tmr**4 ) * epm
  fm4 = (60.D0*lam**2*tmr - 80.D0*lam**3*tmr**3 + 16.D0*lam**4*tmr**5 ) * epm
  tpr = t + r
  exa = - lam * tpr**2
  if (exa .lt. -691.d0) then
     epp = 0.d0
  else
     epp = exp(exa)
  end if
  fp  = tpr * epp
  fp1 = ( 1.D0 - 2.D0*lam*tpr**2 ) * epp
  fp2 = (- 6.D0*lam*tpr + 4.D0*lam**2*tpr**3 ) * epp
  fp3 = (- 6.D0*lam + 24.D0*lam**2*tpr**2 - 8.D0*lam**3*tpr**4 ) * epp
  fp4 = (60.D0*lam**2*tpr - 80.D0*lam**3*tpr**3 + 16.D0*lam**4*tpr**5 ) * epp
!
! Construct coefficients A, B, C, K and L
!
  A = 3.D0 * amp * ( (fm2 - fp2)/r**3 + 3.D0*(fm1 + fp1)/r**4 + &
       3.D0*(fm - fp)/r**5 )
  B = - amp*( (fm3 + fp3)/r**2 + 3.D0*(fm2 - fp2)/r**3 + &
       6.D0*(fm1 + fp1)/r**4 + 6.D0*(fm - fp)/r**5 )
  C = F1o4*amp*( (fm4 - fp4)/r + 2.D0*(fm3 + fp3)/r**2 + &
       9.D0*(fm2 - fp2)/r**3 + 21.D0*(fm1 + fp1)/r**4 + 21.D0*(fm - fp)/r**5 )
  K = amp*( (fm2 - fp2)/r**2 + 3.d0*(fm1+fp1)/r**3 + 3*(fm-fp)/r**4 )
  L = amp*( (fm3+fp3)/r + 2.d0*(fm2-fp2)/r**2 + 3.d0*(fm1+fp1)/r**3 + &
                3.d0*(fm-fp)/r**4 )
!
  return
!
end subroutine coefficients
subroutine dt_coefficients(r,t,amp,lam,Ad,Bd,Cd,Kd,Ld)
  implicit none
  real*8,intent(in)           :: r,t,lam,amp
  real*8,intent(out)          :: Ad,Bd,Cd,Kd,Ld
! Function f of t+r and t-r and its derivatives
  real*8                      :: fm1,fm2,fm3,fm4,fm5
  real*8                      :: fp1,fp2,fp3,fp4,fp5
  real*8                      :: tmr,tpr,epp, epm, exa
  real*8                      :: F1o4
  parameter ( F1o4 = 0.25D0 )
!
! Construct function f and its first four derivatives
!
  tmr = t - r
  exa = - lam * tmr**2
  if (exa .lt. -691.d0) then
     epm = 0.d0
  else
     epm = exp(exa)
  end if
  fm1 = ( 1.D0 - 2.D0*lam*tmr**2 ) * epm
  fm2 = (- 6.D0*lam*tmr + 4.D0*lam**2*tmr**3 ) * epm
  fm3 = (- 6.D0*lam + 24.D0*lam**2*tmr**2 - 8.D0*lam**3*tmr**4 ) * epm
  fm4 = (60.D0*lam**2*tmr - 80.D0*lam**3*tmr**3 + 16.D0*lam**4*tmr**5 ) * epm
  fm5 = (60.d0*lam**2 - 360.d0*lam**3*tmr**2 + 240.d0*lam**4*tmr**4  &
              - 32.d0*lam**5*tmr**6) * epm
  tpr = t + r
  exa = - lam * tpr**2
  if (exa .lt. -691.d0) then
     epp = 0.d0
  else
     epp = exp(exa)
  end if
  fp1 = ( 1.D0 - 2.D0*lam*tpr**2 ) * epp
  fp2 = (- 6.D0*lam*tpr + 4.D0*lam**2*tpr**3 ) * epp
  fp3 = (- 6.D0*lam + 24.D0*lam**2*tpr**2 - 8.D0*lam**3*tpr**4 ) * epp
  fp4 = (60.D0*lam**2*tpr - 80.D0*lam**3*tpr**3 + 16.D0*lam**4*tpr**5 ) * epp
  fp5 = (60.d0*lam**2 - 360.d0*lam**3*tpr**2 + 240.d0*lam**4*tpr**4  &
              - 32.d0*lam**5*tpr**6) * epp
!
! Construct coefficients A, B, C, K and L
!
  Ad = 3.D0 * amp * ( (fm3 - fp3)/r**3 + 3.D0*(fm2 + fp2)/r**4 + &
       3.D0*(fm1 - fp1)/r**5 )
  Bd = - amp*( (fm4 + fp4)/r**2 + 3.D0*(fm3 - fp3)/r**3 + &
       6.D0*(fm2 + fp2)/r**4 + 6.D0*(fm1 - fp1)/r**5 )
  Cd = F1o4*amp*( (fm5 - fp5)/r + 2.D0*(fm4 + fp4)/r**2 + &
       9.D0*(fm3 - fp3)/r**3 + 21.D0*(fm2 + fp2)/r**4 + 21.D0*(fm1 - fp1)/r**5 )
  Kd = amp*( (fm3 - fp3)/r**2 + 3.d0*(fm2+fp2)/r**3 + 3*(fm1-fp1)/r**4 )
  Ld = amp*( (fm4+fp4)/r + 2.d0*(fm3-fp3)/r**2 + 3.d0*(fm2+fp2)/r**3 + &
                3.d0*(fm1-fp1)/r**4 )
!
  return
!
end subroutine dt_coefficients
! Compute d^2C/dt^2
!
subroutine compute_Cdd(t,r,amp,lam,Cdd)
  implicit none
  real*8,intent(in)           :: r,t,lam,amp
  real*8,intent(out)          :: Cdd
! Function f of t+r and t-r and its derivatives
  real*8                      :: fm2,fm3,fm4,fm5,fm6
  real*8                      :: fp2,fp3,fp4,fp5,fp6
  real*8                      :: tmr,tpr,epp, epm, exa
  real*8                      :: F1o4
  parameter ( F1o4 = 0.25D0 )
!
  tmr = t - r
  exa = - lam * tmr**2
  if (exa .lt. -691.d0) then
     epm = 0.d0
  else
     epm = exp(exa)
  end if
  fm2 = (- 6.D0*lam*tmr + 4.D0*lam**2*tmr**3 ) * epm
  fm3 = (- 6.D0*lam + 24.D0*lam**2*tmr**2 - 8.D0*lam**3*tmr**4 ) * epm
  fm4 = (60.D0*lam**2*tmr - 80.D0*lam**3*tmr**3 + 16.D0*lam**4*tmr**5 ) * epm
  fm5 = (60.d0*lam**2 - 360.d0*lam**3*tmr**2 + 240.d0*lam**4*tmr**4  &
              - 32.d0*lam**5*tmr**6) * epm
  fm6 = (-840.d0*lam**3*tmr + 1680.d0*lam**4*tmr**3 - 672.d0*lam**5*tmr**5 &
              + 64.d0*lam**6*tmr**7) * epm
  tpr = t + r
  exa = - lam * tpr**2
  if (exa .lt. -691.d0) then
     epp = 0.d0
  else
     epp = exp(exa)
  end if
  fp2 = (- 6.D0*lam*tpr + 4.D0*lam**2*tpr**3 ) * epp
  fp3 = (- 6.D0*lam + 24.D0*lam**2*tpr**2 - 8.D0*lam**3*tpr**4 ) * epp
  fp4 = (60.D0*lam**2*tpr - 80.D0*lam**3*tpr**3 + 16.D0*lam**4*tpr**5 ) * epp
  fp5 = (60.d0*lam**2 - 360.d0*lam**3*tpr**2 + 240.d0*lam**4*tpr**4  &
              - 32.d0*lam**5*tpr**6) * epp
  fp6 = (-840.d0*lam**3*tpr + 1680.d0*lam**4*tpr**3 - 672.d0*lam**5*tpr**5 &
              + 64.d0*lam**6*tpr**7) * epp
!
  Cdd = F1o4*amp*( (fm6 - fp6)/r + 2.D0*(fm5 + fp5)/r**2 + &
       9.D0*(fm4 - fp4)/r**3 + 21.D0*(fm3 + fp3)/r**4 + 21.D0*(fm2 - fp2)/r**5 )
end subroutine compute_Cdd
! Compute d^2L/dt^2
!
subroutine compute_Ldd(t,r,amp,lam,Ldd)
  implicit none
  real*8,intent(in)           :: r,t,lam,amp
  real*8,intent(out)          :: Ldd
! Function f of t+r and t-r and its derivatives
  real*8                      :: fm2,fm3,fm4,fm5,fm6
  real*8                      :: fp2,fp3,fp4,fp5,fp6
  real*8                      :: tmr,tpr,epp, epm, exa
  real*8                      :: F1o4
  parameter ( F1o4 = 0.25D0 )
!
  tmr = t - r
  exa = - lam * tmr**2
  if (exa .lt. -691.d0) then
     epm = 0.d0
  else
     epm = exp(exa)
  end if
  fm2 = (- 6.D0*lam*tmr + 4.D0*lam**2*tmr**3 ) * epm
  fm3 = (- 6.D0*lam + 24.D0*lam**2*tmr**2 - 8.D0*lam**3*tmr**4 ) * epm
  fm4 = (60.D0*lam**2*tmr - 80.D0*lam**3*tmr**3 + 16.D0*lam**4*tmr**5 ) * epm
  fm5 = (60.d0*lam**2 - 360.d0*lam**3*tmr**2 + 240.d0*lam**4*tmr**4  &
              - 32.d0*lam**5*tmr**6) * epm
  tpr = t + r
  exa = - lam * tpr**2
  if (exa .lt. -691.d0) then
     epp = 0.d0
  else
     epp = exp(exa)
  end if
  fp2 = (- 6.D0*lam*tpr + 4.D0*lam**2*tpr**3 ) * epp
  fp3 = (- 6.D0*lam + 24.D0*lam**2*tpr**2 - 8.D0*lam**3*tpr**4 ) * epp
  fp4 = (60.D0*lam**2*tpr - 80.D0*lam**3*tpr**3 + 16.D0*lam**4*tpr**5 ) * epp
  fp5 = (60.d0*lam**2 - 360.d0*lam**3*tpr**2 + 240.d0*lam**4*tpr**4  &
              - 32.d0*lam**5*tpr**6) * epp
  Ldd = amp*( (fm5+fp5)/r + 2.d0*(fm4-fp4)/r**2 + 3.d0*(fm3+fp3)/r**3 + &
                3.d0*(fm2-fp2)/r**4 )
end subroutine compute_Ldd
! Compute r dC/dt 
!
subroutine compute_Cdotxr(Cdotxr,amp,lam,t,r)
   implicit none
   real*8 :: Cdotxr,amp,lam,t,r
   real*8 :: tmr,epm,fm5,tpr,epp,fp5,exa
!
  tmr = t - r
  exa = -lam*tmr**2
  if (exa .lt. -691.d0) then
     epm = 0.d0
  else
     epm = exp(exa)
  end if
  fm5 = (60.D0*lam**2 - 360.D0*lam**3*tmr**2 + 240.D0*lam**4*tmr**4 &
            - 32.d0*lam**5 * tmr**6 ) * epm
! *** TEST ****
!  tpr = t + r 
!  exa = -lam*tpr**2
!  if (exa .lt. -691.d0) then 
!     epp = 0.d0
!  else
!     epp = exp(exa)
!  end if
!  fp5 = (60.D0*lam**2 - 360.D0*lam**3*tpr**2 + 240.D0*lam**4*tpr**4 &
!            - 32.d0*lam**5 * tpr**6 ) * epp
  fp5 = 0.d0
! ************
  Cdotxr = 0.25d0*amp*(fm5-fp5)
end subroutine compute_Cdotxr
! Compute r dL/dt
!
subroutine compute_Ldotxr(Ldotxr,amp,lam,t,r)
   implicit none
   real*8 :: Ldotxr,amp,lam,t,r
   real*8 :: tmr,epm,fm4,exa
!
  tmr = t - r
  exa = -lam*tmr**2
  if (exa .lt. -691.d0) then
     epm = 0.d0
  else
     epm = exp(exa)
  end if
  fm4 = (60.D0*lam**2 - 80.D0*lam**3*tmr**2 + 16.D0*lam**4*tmr**4)*tmr*epm
  Ldotxr = amp*fm4
end subroutine compute_Ldotxr
! Compute R_{22} = sqrt(Pi/5) F^{(4)}
subroutine R_22(t,r,Rlmr,Rlmi,amp,lam)
   implicit none
   real*8 :: t,r,Rlmr,Rlmi,amp,lam
   real*8 :: tmr,epm,fm4,exa
   real*8, parameter :: pi = 3.14159265358979323846d0
! 
  tmr = t - r
  exa = - lam * tmr**2
  if (exa .lt. -691.d0) then
     epm = 0.d0
  else
     epm = exp(exa)
  end if
  fm4 = (60.D0*lam**2*tmr - 80.D0*lam**3*tmr**3 + 16.D0*lam**4*tmr**5 ) * epm
  Rlmr = amp*fm4*sqrt(0.2d0*pi)
  Rlmi = 0.d0
end subroutine R_22
! Compute Q_{21} = (-i/4) sqrt(8 Pi/15) G^{(4)} 
subroutine Q_21(t,r,Qlmr,Qlmi,amp,lam)
   implicit none
   real*8 :: t,r,Qlmr,Qlmi,amp,lam
   real*8 :: tmr,epm,fm4,exa
   real*8, parameter :: pi = 3.14159265358979323846d0
!
  tmr = t - r
  exa = - lam * tmr**2
  if (exa .lt. -691.d0) then
     epm = 0.d0
  else
     epm = exp(exa)
  end if
  fm4 = (60.D0*lam**2*tmr - 80.D0*lam**3*tmr**3 + 16.D0*lam**4*tmr**5 ) * epm
  Qlmr = 0.d0
  Qlmi = -0.25d0*amp*fm4*sqrt(8.d0*pi/15.d0)
end subroutine Q_21
!
subroutine fij_dij(frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,drt,drp,dtt,dtp,dpp, &
                        SinTh,CosTh,SinPh,CosPh,mode)
   implicit none
   integer :: mode
   real*8 :: frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,SinTh,CosTh,SinPh,CosPh
   real*8 :: drt,drp,dtt,dtp,dpp
!
   frr = 0.d0
   frt = 0.d0
   frp = 0.d0
   f1tt = 0.d0
   f2tt = 0.d0
   ftp = 0.d0
   f1pp = 0.d0
   f2pp = 0.d0
   drt = 0.d0
   drp = 0.d0
   dtt = 0.d0
   dtp = 0.d0
   dpp = 0.d0
!
   if (mode==20) call fij_20(frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,SinTh,CosTh,SinPh,CosPh)
   if (mode==21) call dij_21(drt,drp,dtt,dtp,dpp,SinTh,CosTh,SinPh,CosPh)
   if (mode==22) call fij_22(frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,SinTh,CosTh,SinPh,CosPh)
end subroutine fij_dij
! l=2, m=0
subroutine fij_20(frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,SinTh,CosTh, &
                        SinPh,CosPh)
   implicit none
   real*8 :: frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,SinTh,CosTh,SinPh,CosPh
!  
   frr  = 2.d0 - 3.d0*SinTh**2
   frt  = - 3.d0*(SinTh*CosTh)
   frp = 0.d0
   f1tt = 3.d0*SinTh**2
   f2tt = -1.d0
   ftp = 0.d0
   f1pp = -f1tt
   f2pp = 3.d0*SinTh**2 - 1.d0
end subroutine fij_20
subroutine compute_E_GW(t,r,E_GW,amp,lam,mode)
  implicit none
  integer :: mode
  real*8 :: t,r,E_GW,amp,lam
!
  if (mode==20) call E_GW_20(t,r,E_GW,amp,lam)
  if (mode==21) call E_GW_21(t,r,E_GW,amp,lam)
  if (mode==22) call E_GW_22(t,r,E_GW,amp,lam)
end subroutine compute_E_GW
! In general, dE/dt=r^2 dC/dt/(16 pi)* \int (f1tt^2 + 4 ftp^2) d\Omega 
!                  = fac * r^2 dC/dt 
! where fac = (1/16 pi) \int (f1tt^2 + 4 ftp^2) d\Omega = 6/5 for l=2, m=0
!     
subroutine E_GW_20(t,r,E_GW,amp,lam)
  implicit none
  real*8 :: t,r,E_GW,amp,lam,Cdotxr
  real*8, parameter :: fac = 1.2d0
!
  call compute_Cdotxr(Cdotxr,amp,lam,t,r)
  E_GW = fac * Cdotxr*Cdotxr
end subroutine E_GW_20
! In general, dE/dt=r^2 dL/dt/(16 pi)* \int (dtt^2 + dtp^2) d\Omega
!                  = fac * r^2 dL/dt (for odd modes)
! where fac = (1/16 pi) \int (dtt^2 + dtp^2) d\Omega = 1/10 for l=2, m=1
!
subroutine E_GW_21(t,r,E_GW,amp,lam)
  implicit none
  real*8 :: t,r,E_GW,amp,lam,ldotxr
  real*8, parameter :: fac = 0.1d0
!
  call compute_Ldotxr(Ldotxr,amp,lam,t,r)
  E_GW = fac * Ldotxr*Ldotxr
end subroutine E_GW_21
! l=2, m=2 
subroutine fij_22(frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,SinTh,CosTh, &
                        SinPh,CosPh)
   implicit none
   real*8 :: frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,SinTh,CosTh,SinPh,CosPh
!
   frr  = SinTh*SinTh*(CosPh**2 - SinPh**2)
   frt  = SinTh*CosTh*(CosPh**2 - SinPh**2)
   frp = -2.d0*SinTh*SinPh*CosPh
   f1tt = (1.d0+CosTh*CosTh)*(CosPh**2 - SinPh**2)
   f2tt = -(CosPh**2 - SinPh**2)
   ftp = 2.d0*CosTh*SinPh*CosPh
   f1pp = -f1tt
   f2pp = CosTh*CosTh*(CosPh**2 - SinPh**2)
end subroutine fij_22
! l=2, m=1, odd mode
subroutine dij_21(drt,drp,dtt,dtp,dpp,SinTh,CosTh,SinPh,CosPh)
   implicit none
   real*8 :: drt,drp,dtt,dtp,dpp,SinPh,CosPh,SinTh,CosTh
!
   drt = -2.d0*CosTh*CosPh
   drp = 2.d0*(CosTh**2 - SinTh**2)*SinPh
   dtt = -SinTh*CosPh
   dtp = SinTh*CosTh*SinPh
   dpp = SinTh*CosPh
end subroutine dij_21
! In general, dE/dt=r^2 dC/dt/(16 pi)* \int (f1tt^2 + 4 ftp^2) d\Omega
!                  = fac * r^2 dC/dt
! where fac = (1/16 pi) \int (f1tt^2 + 4 ftp^2) d\Omega 
!           = 2/5 for l=2, m=2
!
subroutine E_GW_22(t,r,E_GW,amp,lam)
  implicit none
  real*8 :: t,r,E_GW,amp,lam,Cdotxr
  real*8, parameter :: fac = 0.4d0
!
  call compute_Cdotxr(Cdotxr,amp,lam,t,r)
  E_GW = fac * Cdotxr*Cdotxr
end subroutine E_GW_22
subroutine lin_wave_analytic_gxx(ex, X, Y, Z, T, Amp, Width, &
     gxx,PhysR, dRdr,mode)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: PhysR, dRdr
  real*8, dimension(ex(1),ex(2),ex(3)),target :: gxx
  real*8                                      :: T, Amp, Width
  integer                                      :: mode
  integer                            :: i,j,k,l,m,n,o
  real*8                             :: R,Rho,SinTh,CosTh,SinPh,CosPh,gxxL,gxyL,gxzL,gyyL,gyzL,gzzL,detL,thetaangle,phiangle
  real*8                             :: A,B,C, KK,LL
!  real*8                             :: Ad,Bd,Cd, KKd,LLd
  real*8, dimension(3,3)             :: trans,g_pol, k_pol
  real*8                             :: frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp
  real*8                              :: drt,drp,dtt,dtp,dpp
  real*8                             :: F1o3,ZERO, ONE, TWO, THREE
  parameter ( F1o3 = 1.D0 / 3.D0, ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0 )
  !
  ! go to each gridpoint...
  !
  do k = 1,ex(3)
     do j = 1,ex(2)
        do i = 1,ex(1)
           !
           ! find polar coordinates from cartesian coordinates
           !
           R     = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           Rho   = sqrt(X(i,j,k)**2 + Y(i,j,k)**2)
           SinPh = Y(i,j,k) / Rho
           CosPh = X(i,j,k) / Rho
           SinTh = Rho / R
           CosTh = Z(i,j,k) / R
!!$           !
!!$
!!$           ! Heres a new formulation for the trans() matrix and Sin/Cos Th/Ph, so that
!!$           !   we dont get NaNs from Rho=0 points
!!$           thetaangle=atan2(Rho,z(i,j,k))
!!$           phiangle=atan2(y(i,j,k),x(i,j,k))
!!$
!!$           SinTh = Rho/R
!!$           thetaangle=asin(SinTh)
!!$           CosTh = cos(thetaangle)
!!$           SinPh = sin(phiangle)
!!$           CosPh = cos(phiangle)
!!$           !
!!$           ! find transformation matrix trans_{ij} = dx^{i}/dx^j, where
!!$           ! primes are polar coordinates
!!$           !
!!$           trans(1,1) = SinTh*CosPh
!!$           trans(1,2) = SinTh*SinPh
!!$           trans(1,3) = CosTh
!!$           trans(2,1) = CosTh*CosPh/R
!!$           trans(2,2) = CosTh*SinPh/R
!!$           trans(2,3) = - SinTh/R
!!$           trans(3,1) = - SinPh/(R*SinTh)
!!$           trans(3,2) = CosPh/(R*SinTh)
!!$           trans(3,3) = ZERO
           trans(1,1) = X(i,j,k)/R
           trans(1,2) = Y(i,j,k)/R
           trans(1,3) = Z(i,j,k)/R
           trans(2,1) = X(i,j,k)*Z(1,1,k)/(R**2 * Rho)
           trans(2,2) = Y(i,j,k)*Z(1,1,k)/(R**2 * Rho)
           trans(2,3) = - Rho/R**2
           trans(3,1) = - Y(i,j,k)/Rho**2
           trans(3,2) = X(i,j,k)/Rho**2
           trans(3,3) = ZERO
           ! 
           !
           ! find coefficients A, B, C (eq. (6) ) and K and L (eq. (9))
           !
           call coefficients(PhysR(i,j,k),T,Amp,Width,A,B,C,KK,LL)
           !
           ! find coefficients f and d (eqs. (7)  and (10))
           ! 
           call fij_dij(frr,frt,frp,f1tt,f2tt,ftp,f1pp,f2pp,drt,drp,dtt,dtp,dpp, &
                SinTh,CosTh,SinPh,CosPh,mode)
           !
           ! Find metric in polar coordinates (eq. (5) )
           !
           g_pol(1,1) = (ONE + A*frr) * dRdr(i,j,k)**2
           g_pol(1,2) = (B*frt+KK*drt)*PhysR(i,j,k) * dRdr(i,j,k)
           g_pol(1,3) = (B*frp+KK*drp)*PhysR(i,j,k)*SinTh * dRdr(i,j,k)
           g_pol(2,2) = (ONE + C*f1tt + A*f2tt + LL*dtt) * PhysR(i,j,k)**2
           g_pol(2,3) = ( (A-2.d0*C)*ftp + LL*dtp)*SinTh*PhysR(i,j,k)**2
           g_pol(3,3) = (ONE + C*f1pp + A*f2pp + LL*dpp)*(PhysR(i,j,k)*SinTh)**2
           g_pol(2,1) = g_pol(1,2)
           g_pol(3,1) = g_pol(1,3)
           g_pol(3,2) = g_pol(2,3)
           ! Transform metric into cartesian coordinates
           !
           gxxL = ZERO
           gxyL = ZERO
           gxzL = ZERO
           gyyL = ZERO
           gyzL = ZERO
           gzzL = ZERO
           do n = 1,3
              do o = 1,3
                 gxxL = gxxL + &
                      trans(n,1)*trans(o,1)*g_pol(n,o)
                 gxyL = gxyL + &
                      trans(n,1)*trans(o,2)*g_pol(n,o)
                 gxzL = gxzL + &
                      trans(n,1)*trans(o,3)*g_pol(n,o)
                 gyyL = gyyL + &
                      trans(n,2)*trans(o,2)*g_pol(n,o)
                 gyzL = gyzL + &
                      trans(n,2)*trans(o,3)*g_pol(n,o)
                 gzzL = gzzL + &
                      trans(n,3)*trans(o,3)*g_pol(n,o)
              end do
           end do
           detL = gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL &
                - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL
           gxx(i,j,k) = gxxL / detL**F1o3
        end do
     end do
  end do
end subroutine lin_wave_analytic_gxx
