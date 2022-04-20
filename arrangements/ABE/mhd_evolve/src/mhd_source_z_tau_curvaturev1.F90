!-----------------------------------------------------------------------------
!
! routines for evolving the conducting fluid
! in the presence of a magnetic field
!
!-----------------------------------------------------------------------------
!

subroutine print(ext,X,rho_star,tau,st_x,st_y,st_z)
  implicit none
  integer, dimension(3)				:: ext
  real*8, dimension(ext(1),ext(2),ext(3))	:: X
  integer					:: i,imin,imax
  real*8, dimension(ext(1),ext(2),ext(3))	:: rho_star,tau
  real*8, dimension(ext(1),ext(2),ext(3))        :: st_x,st_y,st_z
  !
  imin = lbound(st_x,1)
  imax = ubound(st_x,1)
  do i=imin,imax
     write(*,*) X(i,1,1),st_x(i,2,2),st_y(i,2,2),st_z(i,2,2),'  st_i'
  end do
  write(*,*) ' '
end subroutine print

!-----------------------------------------------------------------------------
! 
! Add third-order accurate curvature terms to S_i_rhs
!
!-----------------------------------------------------------------------------
subroutine mhd_source_z_tau(ex, dX, dY, dZ, mhd_st_x_rhs, mhd_st_y_rhs, mhd_st_z_rhs, &
     tau_rhs,rho_star, h, u0, vx, vy, vz, P, sbt,sbx,sby,sbz, &
     alpha, betax, betay, betaz, phi, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Symmetry, &
     alpha_f, betax_f, betay_f, betaz_f, gxx_f, gxy_f, gxz_f, &
     gyy_f, gyz_f, gzz_f, phi_f, m,Z, enable_HARM_energyvariable)
  implicit none

  interface
     subroutine ddx_o3(ext,f_f,df,dX,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f_f,df
       real*8                                   :: sym,dX
     end subroutine ddx_o3
     subroutine ddy_o3(ext,f_f,df,dY,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f_f,df
       real*8                                   :: sym,dY
     end subroutine ddy_o3
     subroutine ddz_o3(ext,f_f,df,dZ,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f_f,df
       real*8                                   :: sym,dZ
     end subroutine ddz_o3
  end interface

  integer, dimension(3)                            :: ex
  real*8                                           :: dX, dY, dZ, sfpi
  real*8, dimension(ex(1),ex(2),ex(3))             :: mhd_st_x_rhs, mhd_st_y_rhs, mhd_st_z_rhs
  real*8, dimension(ex(1),ex(2),ex(3))             :: rho_star, h, u0,tau_rhs
  real*8, dimension(ex(1),ex(2),ex(3))             :: vx, vy, vz, P, T_fluid, M_B
  real*8, dimension(ex(1),ex(2),ex(3))             :: sbt,sbx,sby,sbz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: alpha, phi
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gxx, gxy, gxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gyy, gyz, gzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: alpha_f, phi_f
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: betax_f, betay_f, betaz_f
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gxx_f, gxy_f, gxz_f
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gyy_f, gyz_f, gzz_f
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3))		   :: Z
  integer                                     :: Symmetry,m,enable_HARM_energyvariable
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  integer                            :: AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  !
  ! Other variables:
  ! 
  real*8, dimension(ex(1),ex(2),ex(3)) :: g4tt, g4tx, g4ty, g4tz
  real*8, dimension(ex(1),ex(2),ex(3)) :: g4xx, g4xy, g4xz
  real*8, dimension(ex(1),ex(2),ex(3)) :: g4yy, g4yz, g4zz
  real*8, dimension(ex(1),ex(2),ex(3)) :: g4tt_f, g4tx_f, g4ty_f, g4tz_f
  real*8, dimension(ex(1),ex(2),ex(3)) :: g4xx_f, g4xy_f, g4xz_f
  real*8, dimension(ex(1),ex(2),ex(3)) :: g4yy_f, g4yz_f, g4zz_f
  real*8, dimension(ex(1),ex(2),ex(3)) :: g4ttm, g4txm, g4tym, g4tzm
  real*8, dimension(ex(1),ex(2),ex(3)) :: g4xxm, g4xym, g4xzm
  real*8, dimension(ex(1),ex(2),ex(3)) :: g4yym, g4yzm, g4zzm
  real*8, dimension(ex(1),ex(2),ex(3)) :: Psi4, Psi6, al
  real*8, dimension(ex(1),ex(2),ex(3)) :: Psi4_f, al_f, b2
  integer			       :: kmin
  !
  real*8, parameter                    :: ZERO = 0.D0, HALF = 0.5D0, ONE  = 1.D0
  real*8, parameter                    :: TWO  = 2.D0, SYM  = 1.D0, ANTI = -1.D0
  !
  sfpi = sqrt(4.d0*acos(-1.d0))
  !
  ! Input translation
  !
  Psi6  = exp(phi)
  Psi4 = Psi6 * Psi6 * Psi6 * Psi6
  Psi6  = Psi4 * Psi6 * Psi6
  al = ONE + alpha
  Psi4_f = exp(4.d0*phi_f)
  al_f = ONE + alpha_f

  !
  ! Divide b^{\mu} by alpha sqrt(4 pi) 
  !
  sbt = sbt/al/sfpi
  sbx = sbx/al/sfpi
  sby = sby/al/sfpi
  sbz = sbz/al/sfpi
  !
  ! Compute b^2 
  !
  b2 = -(al*sbt)**2 + Psi4*( gxx*(sbx+betax*sbt)**2 + &
       2.d0*gxy*(sbx+betax*sbt)*(sby+betay*sbt) + &
       2.d0*gxz*(sbx+betax*sbt)*(sbz+betaz*sbt) + &
       gyy*(sby+betay*sbt)**2 + &
       2.d0*gyz*(sby+betay*sbt)*(sbz+betaz*sbt) + &
       gzz*(sbz+betaz*sbt)**2 ) 

  !----------------------------------------------------------------------------
  ! Compute the 4-metric
  !----------------------------------------------------------------------------
  g4tt = -al**2 + Psi4*(gxx*betax**2 + TWO*gxy*betax*betay + TWO*gxz*betax*betaz + &
       gyy*betay**2 + TWO*gyz*betay*betaz + gzz*betaz**2)
  g4tx = Psi4*(gxx*betax + gxy*betay + gxz*betaz)
  g4ty = Psi4*(gxy*betax + gyy*betay + gyz*betaz)
  g4tz = Psi4*(gxz*betax + gyz*betay + gzz*betaz)
  g4xx = Psi4*gxx
  g4xy = Psi4*gxy
  g4xz = Psi4*gxz 
  g4yy = Psi4*gyy
  g4yz = Psi4*gyz
  g4zz = Psi4*gzz

  g4tt_f = -al_f**2 + Psi4_f*(gxx_f*betax_f**2 + TWO*gxy_f*betax_f*betay_f + TWO*gxz_f*betax_f*betaz_f + &
       gyy_f*betay_f**2 + TWO*gyz_f*betay_f*betaz_f + gzz_f*betaz_f**2)
  g4tx_f = Psi4_f*(gxx_f*betax_f + gxy_f*betay_f + gxz_f*betaz_f)
  g4ty_f = Psi4_f*(gxy_f*betax_f + gyy_f*betay_f + gyz_f*betaz_f)
  g4tz_f = Psi4_f*(gxz_f*betax_f + gyz_f*betay_f + gzz_f*betaz_f)
  g4xx_f = Psi4_f*gxx_f
  g4xy_f = Psi4_f*gxy_f
  g4xz_f = Psi4_f*gxz_f 
  g4yy_f = Psi4_f*gyy_f
  g4yz_f = Psi4_f*gyz_f
  g4zz_f = Psi4_f*gzz_f

  if(m==1) then
     !-----------------------------------------------------------------------------
     ! Compute the addition to mhd_st_x_rhs and tau_rhs
     !-----------------------------------------------------------------------------

     ! take derivatives of the four-metric:

     if (Symmetry == OCTANT) then
        call ddx_o3(ex,g4tx_f,g4txm,dX,ANTI)
        call ddx_o3(ex,g4xy_f,g4xym,dX,ANTI)
        call ddx_o3(ex,g4xz_f,g4xzm,dX,ANTI)
        call ddx_o3(ex,g4ty_f,g4tym,dX,SYM)
        call ddx_o3(ex,g4yz_f,g4yzm,dX,SYM)
     elseif (Symmetry==AXISYM) then
        call ddx_o3(ex,g4tx_f,g4txm,dX,ANTI)
        call ddx_o3(ex,g4xy_f,g4xym,dX,SYM)
        call ddx_o3(ex,g4xz_f,g4xzm,dX,ANTI)
        call ddx_o3(ex,g4ty_f,g4tym,dX,ANTI)
        call ddx_o3(ex,g4yz_f,g4yzm,dX,ANTI)
     else
        call ddx_o3(ex,g4tx_f,g4txm,dX,SYM)
        call ddx_o3(ex,g4xy_f,g4xym,dX,SYM)
        call ddx_o3(ex,g4xz_f,g4xzm,dX,SYM)
        call ddx_o3(ex,g4ty_f,g4tym,dX,SYM)
        call ddx_o3(ex,g4yz_f,g4yzm,dX,SYM)
     end if

     call ddx_o3(ex,g4tt_f,g4ttm,dX,SYM)
     call ddx_o3(ex,g4tz_f,g4tzm,dX,SYM)
     call ddx_o3(ex,g4xx_f,g4xxm,dX,SYM)
     call ddx_o3(ex,g4yy_f,g4yym,dX,SYM)
     call ddx_o3(ex,g4zz_f,g4zzm,dX,SYM)

     ! add \sqrt{-g}/2 * T^{00}\partial_x g_{00}
     mhd_st_x_rhs = mhd_st_x_rhs + HALF*( rho_star*h*u0 + al*Psi6*b2*u0**2 &
          - Psi6/al*(P+HALF*b2) -al*Psi6*sbt**2)*g4ttm

     ! add \sqrt{-g} * T^{0i}\partial_x g_{0i} 
     mhd_st_x_rhs = mhd_st_x_rhs + ( rho_star*h*u0*vx + al*Psi6*b2*vx*u0**2 + & 
          Psi6*betax/al*(P+HALF*b2) - al*Psi6*sbt*sbx)*g4txm + &
          (rho_star*h*u0*vy + al*Psi6*b2*vy*u0**2 + &
          Psi6*betay/al*(P+HALF*b2) - al*Psi6*sbt*sby)*g4tym + &
          (rho_star*h*u0*vz + al*Psi6*b2*vz*u0**2 + &
          Psi6*betaz/al*(P+HALF*b2) - al*Psi6*sbt*sbz)*g4tzm

     ! add \sqrt{-g}/2 * T^{ij}\partial_x g_{ij}
     mhd_st_x_rhs = mhd_st_x_rhs + HALF*(rho_star*h*u0*vx*vx + al*Psi6*b2*(vx*u0)**2 + &
          (P+HALF*b2)*(al*Psi6*gupxx/Psi4 - betax**2*Psi6/al) &
          - al*Psi6*sbx**2)*g4xxm

     mhd_st_x_rhs = mhd_st_x_rhs + (rho_star*h*u0*vx*vy + al*Psi6*b2*vx*vy*u0**2 + &
          (P+HALF*b2)*(al*Psi6*gupxy/Psi4 - betax*betay*Psi6/al) &
          - al*Psi6*sbx*sby)*g4xym

     mhd_st_x_rhs = mhd_st_x_rhs + (rho_star*h*u0*vx*vz + al*Psi6*b2*vx*vz*u0**2 + &
          (P+HALF*b2)*(al*Psi6*gupxz/Psi4 - betax*betaz*Psi6/al) &
          - al*Psi6*sbx*sbz)*g4xzm

     mhd_st_x_rhs = mhd_st_x_rhs + HALF*(rho_star*h*u0*vy*vy + al*Psi6*b2*(vy*u0)**2 + &
          (P+HALF*b2)*(al*Psi6*gupyy/Psi4 - betay**2*Psi6/al) &
          - al*Psi6*sby**2)*g4yym

     mhd_st_x_rhs = mhd_st_x_rhs + (rho_star*h*u0*vy*vz + al*Psi6*b2*vy*vz*u0**2 + &
          (P+HALF*b2)*(al*Psi6*gupyz/Psi4 - betay*betaz*Psi6/al) & 
          - al*Psi6*sby*sbz)*g4yzm

     mhd_st_x_rhs = mhd_st_x_rhs + HALF*(rho_star*h*u0*vz*vz + al*Psi6*b2*(vz*u0)**2 + &
          (P+HALF*b2)*(al*Psi6*gupzz/Psi4 - betaz**2*Psi6/al) & 
          - al*Psi6*sbz**2)*g4zzm


     ! add -\sqrt{-g} * (T^{00} \betax + T^{0x})\partial_x alpha to tau_rhs
     if(enable_HARM_energyvariable==0) then
        call ddx_o3(ex,alpha_f,g4ttm,dX,SYM)
        tau_rhs = tau_rhs - ( (rho_star*h*u0 + al*Psi6*b2*u0**2 &
             - al*Psi6*sbt**2)*betax + &
             rho_star*h*u0*vx + al*Psi6*b2*vx*u0**2 - al*Psi6*sbt*sbx )*g4ttm
     end if

  else if(m==2) then
     !-----------------------------------------------------------------------------
     ! Compute the addition to mhd_st_y_rhs and tau_rhs
     !-----------------------------------------------------------------------------
     if (Symmetry .ne. AXISYM) then
        if (Symmetry == OCTANT) then
           call ddy_o3(ex,g4ty_f,g4tym,dY,ANTI)
           call ddy_o3(ex,g4xy_f,g4xym,dY,ANTI)
           call ddy_o3(ex,g4yz_f,g4yzm,dY,ANTI)
        else
           call ddy_o3(ex,g4ty_f,g4tym,dY,SYM)
           call ddy_o3(ex,g4xy_f,g4xym,dY,SYM)
           call ddy_o3(ex,g4yz_f,g4yzm,dY,SYM)
        end if

        call ddy_o3(ex,g4tt_f,g4ttm,dY,SYM)
        call ddy_o3(ex,g4tx_f,g4txm,dY,SYM)
        call ddy_o3(ex,g4tz_f,g4tzm,dY,SYM)
        call ddy_o3(ex,g4xx_f,g4xxm,dY,SYM)
        call ddy_o3(ex,g4xz_f,g4xzm,dY,SYM)
        call ddy_o3(ex,g4yy_f,g4yym,dY,SYM)
        call ddy_o3(ex,g4zz_f,g4zzm,dY,SYM)

        ! add \sqrt{-g}/2 * T^{00}\partial_y g_{00}
        mhd_st_y_rhs = mhd_st_y_rhs + HALF*( rho_star*h*u0 + al*Psi6*b2*u0**2 &
             - Psi6/al*(P+HALF*b2) -al*Psi6*sbt**2)*g4ttm

        ! add \sqrt{-g} * T^{0i}\partial_y g_{0i} 
        mhd_st_y_rhs = mhd_st_y_rhs + ( rho_star*h*u0*vx + al*Psi6*b2*vx*u0**2 + &
             Psi6*betax/al*(P+HALF*b2) - al*Psi6*sbt*sbx)*g4txm + &
             (rho_star*h*u0*vy + al*Psi6*b2*vy*u0**2 + &
             Psi6*betay/al*(P+HALF*b2) - al*Psi6*sbt*sby)*g4tym + &
             (rho_star*h*u0*vz + al*Psi6*b2*vz*u0**2 + &
             Psi6*betaz/al*(P+HALF*b2) - al*Psi6*sbt*sbz)*g4tzm

        ! add \sqrt{-g}/2 * T^{ij}\partial_y g_{ij}
        mhd_st_y_rhs = mhd_st_y_rhs + HALF*(rho_star*h*u0*vx*vx + al*Psi6*b2*(vx*u0)**2 + &
             (P+HALF*b2)*(al*Psi6*gupxx/Psi4 - betax**2*Psi6/al) &
             - al*Psi6*sbx**2)*g4xxm
        mhd_st_y_rhs = mhd_st_y_rhs + (rho_star*h*u0*vx*vy + al*Psi6*b2*vx*vy*u0**2 + &
             (P+HALF*b2)*(al*Psi6*gupxy/Psi4 - betax*betay*Psi6/al) &
             - al*Psi6*sbx*sby)*g4xym
        mhd_st_y_rhs = mhd_st_y_rhs + (rho_star*h*u0*vx*vz + al*Psi6*b2*vx*vz*u0**2 + &
             (P+HALF*b2)*(al*Psi6*gupxz/Psi4 - betax*betaz*Psi6/al) &
             - al*Psi6*sbx*sbz)*g4xzm
        mhd_st_y_rhs = mhd_st_y_rhs + HALF*(rho_star*h*u0*vy*vy + al*Psi6*b2*(vy*u0)**2 + &
             (P+HALF*b2)*(al*Psi6*gupyy/Psi4 - betay**2*Psi6/al) &
             - al*Psi6*sby**2)*g4yym
        mhd_st_y_rhs = mhd_st_y_rhs + (rho_star*h*u0*vy*vz + al*Psi6*b2*vy*vz*u0**2 + &
             (P+HALF*b2)*(al*Psi6*gupyz/Psi4 - betay*betaz*Psi6/al) &
             - al*Psi6*sby*sbz)*g4yzm
        mhd_st_y_rhs = mhd_st_y_rhs + HALF*(rho_star*h*u0*vz*vz + al*Psi6*b2*(vz*u0)**2 + &
             (P+HALF*b2)*(al*Psi6*gupzz/Psi4 - betaz**2*Psi6/al) &
             - al*Psi6*sbz**2)*g4zzm

        ! add -\sqrt{-g} * (T^{00} \betay + T^{0y})\partial_y alpha to tau_rhs
        if(enable_HARM_energyvariable==0) then
           call ddy_o3(ex,alpha_f,g4ttm,dY,SYM)
           tau_rhs = tau_rhs - ( (rho_star*h*u0 + al*Psi6*b2*u0**2 &
                - al*Psi6*sbt**2)*betay + &
                rho_star*h*u0*vy + al*Psi6*b2*vy*u0**2 - al*Psi6*sbt*sby )*g4ttm
        end if
     end if

  else if(m==3) then
     !-----------------------------------------------------------------------------
     ! Compute the addition to st_z_rhs 
     !-----------------------------------------------------------------------------
     kmin = lbound(phi,3)
     if (Symmetry .ne. NO_SYMM .and. Z(1,1,kmin+1) .gt. 0.d0) then
        call ddz_o3(ex,g4tz_f,g4tzm,dZ,ANTI)
        call ddz_o3(ex,g4xz_f,g4xzm,dZ,ANTI)
        call ddz_o3(ex,g4yz_f,g4yzm,dZ,ANTI)
     else
        call ddz_o3(ex,g4tz_f,g4tzm,dZ,SYM)
        call ddz_o3(ex,g4xz_f,g4xzm,dZ,SYM)
        call ddz_o3(ex,g4yz_f,g4yzm,dZ,SYM)
     end if

     call ddz_o3(ex,g4tt_f,g4ttm,dZ,SYM)
     call ddz_o3(ex,g4tx_f,g4txm,dZ,SYM)
     call ddz_o3(ex,g4ty_f,g4tym,dZ,SYM)
     call ddz_o3(ex,g4xx_f,g4xxm,dZ,SYM)
     call ddz_o3(ex,g4xy_f,g4xym,dZ,SYM)
     call ddz_o3(ex,g4yy_f,g4yym,dZ,SYM)
     call ddz_o3(ex,g4zz_f,g4zzm,dZ,SYM)

     ! add \sqrt{-g}/2 * T^{00}\partial_z g_{00}
     mhd_st_z_rhs = mhd_st_z_rhs + HALF*( rho_star*h*u0 + al*Psi6*b2*u0**2 &
          - Psi6/al*(P+HALF*b2) -al*Psi6*sbt**2)*g4ttm

     !write(*,*) "g1: ",mhd_st_z_rhs(30,30,2)

     ! add \sqrt{-g} * T^{0i}\partial_z g_{0i} 
     mhd_st_z_rhs = mhd_st_z_rhs + ( rho_star*h*u0*vx + al*Psi6*b2*vx*u0**2 + &
          Psi6*betax/al*(P+HALF*b2) - al*Psi6*sbt*sbx)*g4txm + &
          (rho_star*h*u0*vy + al*Psi6*b2*vy*u0**2 + &
          Psi6*betay/al*(P+HALF*b2) - al*Psi6*sbt*sby)*g4tym + &
          (rho_star*h*u0*vz + al*Psi6*b2*vz*u0**2 + &
          Psi6*betaz/al*(P+HALF*b2) - al*Psi6*sbt*sbz)*g4tzm

     ! add \sqrt{-g}/2 * T^{ij}\partial_z g_{ij}
     mhd_st_z_rhs = mhd_st_z_rhs + HALF*(rho_star*h*u0*vx*vx + al*Psi6*b2*(vx*u0)**2 + &
          (P+HALF*b2)*(al*Psi6*gupxx/Psi4 - betax**2*Psi6/al) &
          - al*Psi6*sbx**2)*g4xxm
     mhd_st_z_rhs = mhd_st_z_rhs + (rho_star*h*u0*vx*vy + al*Psi6*b2*vx*vy*u0**2 + &
          (P+HALF*b2)*(al*Psi6*gupxy/Psi4 - betax*betay*Psi6/al) &
          - al*Psi6*sbx*sby)*g4xym
     mhd_st_z_rhs = mhd_st_z_rhs + (rho_star*h*u0*vx*vz + al*Psi6*b2*vx*vz*u0**2 + &
          (P+HALF*b2)*(al*Psi6*gupxz/Psi4 - betax*betaz*Psi6/al) &
          - al*Psi6*sbx*sbz)*g4xzm
     mhd_st_z_rhs = mhd_st_z_rhs + HALF*(rho_star*h*u0*vy*vy + al*Psi6*b2*(vy*u0)**2 + &
          (P+HALF*b2)*(al*Psi6*gupyy/Psi4 - betay**2*Psi6/al) &
          - al*Psi6*sby**2)*g4yym
     mhd_st_z_rhs = mhd_st_z_rhs + (rho_star*h*u0*vy*vz + al*Psi6*b2*vy*vz*u0**2 + &
          (P+HALF*b2)*(al*Psi6*gupyz/Psi4 - betay*betaz*Psi6/al) &
          - al*Psi6*sby*sbz)*g4yzm
     mhd_st_z_rhs = mhd_st_z_rhs + HALF*(rho_star*h*u0*vz*vz + al*Psi6*b2*(vz*u0)**2 + &
          (P+HALF*b2)*(al*Psi6*gupzz/Psi4 - betaz**2*Psi6/al) &
          - al*Psi6*sbz**2)*g4zzm

     ! add -\sqrt{-g} * (T^{00} \betaz + T^{0z})\partial_z alpha to tau_rhs
     if (enable_HARM_energyvariable==0) then
        call ddz_o3(ex,alpha_f,g4ttm,dZ,SYM)
        tau_rhs = tau_rhs - ( (rho_star*h*u0 + al*Psi6*b2*u0**2 &
             - al*Psi6*sbt**2)*betaz + &
             rho_star*h*u0*vz + al*Psi6*b2*vz*u0**2 - al*Psi6*sbt*sbz )*g4ttm
     end if
  end if

  !
  ! Multiply b^{\mu} by alpha sqrt(4 pi)
  !
  sbt = sbt*al*sfpi
  sbx = sbx*al*sfpi
  sby = sby*al*sfpi
  sbz = sbz*al*sfpi


  ! TEST: compute the fluid temperature here:

  !T_fluid = P * M_B / rho_star

  !write(*,*) "P, M_B, rho_star, T_fluid are", P, M_B, rho_star, T_fluid

end subroutine mhd_source_z_tau

!--------------------------------------------------------------------
!
! Add the extrinsic curvature terms to tau_rhs 
!
!--------------------------------------------------------------------
!
subroutine mhd_tau_curvature(ext,tau_rhs,rho_star,h,P,sbt,sbx,sby,sbz, &
     u0,vx,vy,vz,alpha,betax,betay,betaz, &
     phi,gxx,gxy,gxz,gyy,gyz,gzz, &
     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
     K,Axx,Axy,Axz,Ayy,Ayz,Azz,dX,dY,dZ)
  implicit none
  integer, dimension(3)                            :: ext
  real*8                                           :: dX, dY, dZ, sfpi
  real*8, dimension(ext(1),ext(2),ext(3))          :: tau_rhs
  real*8, dimension(ext(1),ext(2),ext(3))          :: rho_star, h, P
  real*8, dimension(ext(1),ext(2),ext(3))          :: sbt,sbx,sby,sbz,b2
  real*8, dimension(ext(1),ext(2),ext(3))          :: u0,vx,vy,vz
  real*8, dimension(ext(1),ext(2),ext(3))          :: alpha,betax,betay,betaz
  real*8, dimension(ext(1),ext(2),ext(3))          :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))          :: gupxx,gupxy,gupxz
  real*8, dimension(ext(1),ext(2),ext(3))          :: gupyy,gupyz,gupzz
  real*8, dimension(ext(1),ext(2),ext(3))          :: K,Axx,Axy,Axz,Ayy,Ayz,Azz
  !
  real*8, dimension(ext(1),ext(2),ext(3))          :: al,sqrtg,Psi4
  !
  real*8, parameter				   :: THIRD = 1.d0/3.d0
  !

  sfpi = sqrt(4.d0*acos(-1.d0)) 
  al = 1.d0+alpha
  Psi4 = exp(4.d0*phi)
  sqrtg = al*exp(6.d0*phi)
  ! 
  ! Divide b^a by alpha sqrt(4 pi)
  !
  sbt = sbt/al/sfpi
  sbx = sbx/al/sfpi
  sby = sby/al/sfpi
  sbz = sbz/al/sfpi
  !
  ! Compute b^2 
  !
  b2 = -(al*sbt)**2 + Psi4*( gxx*(sbx+betax*sbt)**2 + &
       2.d0*gxy*(sbx+betax*sbt)*(sby+betay*sbt) + &
       2.d0*gxz*(sbx+betax*sbt)*(sbz+betaz*sbt) + &
       gyy*(sby+betay*sbt)**2 + &
       2.d0*gyz*(sby+betay*sbt)*(sbz+betaz*sbt) + &
       gzz*(sbz+betaz*sbt)**2 )
  !
  ! Now add extrinsic curvature terms to tau_rhs
  !
  !  write(*,*) "inside mhd_tau_curv0", tau_rhs(16,2,2),rho_star(16,2,2),h(16,2,2),u0(16,2,2),sqrtg(16,2,2),b2(16,2,2),u0(16,2,2),Psi4(16,2,2),Axx(16,2,2),gxx(16,2,2),K(16,2,2),vx(16,2,2),betax(16,2,2),vx(16,2,2),vy(16,2,2),vz(16,2,2)
  tau_rhs = tau_rhs + (rho_star*h*u0 + sqrtg*b2*u0**2)*Psi4* ( &
       (Axx+THIRD*gxx*K)*(vx+betax)**2 + &
       2.d0*(Axy+THIRD*gxy*K)*(vx+betax)*(vy+betay) + &
       2.d0*(Axz+THIRD*gxz*K)*(vx+betax)*(vz+betaz) + &
       (Ayy+THIRD*gyy*K)*(vy+betay)**2 + &
       2.d0*(Ayz+THIRD*gyz*K)*(vy+betay)*(vz+betaz) + &
       (Azz+THIRD*gzz*K)*(vz+betaz)**2 )
  !  write(*,*) "inside mhd_tau_curv", tau_rhs(16,2,2),rho_star(16,2,2),h(16,2,2),u0(16,2,2),sqrtg(16,2,2),b2(16,2,2),u0(16,2,2),Psi4(16,2,2),Axx(16,2,2),gxx(16,2,2),K(16,2,2),vx(16,2,2),betax(16,2,2),vx(16,2,2),vy(16,2,2),vz(16,2,2)

  tau_rhs = tau_rhs - sqrtg*Psi4* ( &
       (Axx+THIRD*gxx*K)*(sbx+sbt*betax)**2 + &
       2.d0*(Axy+THIRD*gxy*K)*(sbx+sbt*betax)*(sby+sbt*betay) + &
       2.d0*(Axz+THIRD*gxz*K)*(sbx+sbt*betax)*(sbz+sbt*betaz) + &
       (Ayy+THIRD*gyy*K)*(sby+sbt*betay)**2 + &
       2.d0*(Ayz+THIRD*gyz*K)*(sby+sbt*betay)*(sbz+sbt*betaz) + &
       (Azz+THIRD*gzz*K)*(sbz+sbt*betaz)**2 )

  tau_rhs = tau_rhs + sqrtg*(P+0.5d0*b2)*K

  !
  ! Multiply b^a by alpha sqrt(4 pi)
  !
  sbt = sbt*al*sfpi
  sbx = sbx*al*sfpi
  sby = sby*al*sfpi
  sbz = sbz*al*sfpi

end subroutine mhd_tau_curvature
