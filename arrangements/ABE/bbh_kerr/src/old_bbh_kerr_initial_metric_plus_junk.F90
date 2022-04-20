!-----------------------------------------------------------------------------
!
! Non-boosted rotating Kerr-Schild metric initial data + "junk" in the BH interior
!
!-----------------------------------------------------------------------------
  subroutine ks_initial_metric_plus_junk_one_BH(ex, x, y, z,       &
                    gxxm1, gxy, gxz, gyym1, gyz, gzzm1,             &
                    Axx, Axy, Axz, Ayy, Ayz, Azz, alpm1, Bx,By,Bz, &
                    m0, sam, r_junk, lambda, x0,y0,z0)
  implicit none

!
! Input parameters:
!
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z
  real*8,                    intent(in) :: sam,m0
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxxm1,gxy,gxz,gyym1,gyz,gzzm1
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: alpm1,Bx,By,Bz

!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: el(0:3),boo(0:3),phh(0:3),eta(0:3,0:3),pel(0:3,0:3)
  real*8  :: k1,k2,k3,aa,hh,r2,w2,W,wx,wy,wz,tmp,tmq,tmr
  real*8  :: r_junk, lambda, x0,y0,z0
  real*8  :: w0, wh_out, wh_in,a
  real*8  :: xl,yl,zl
!  real*8,  dimension(ex(1),ex(2),ex(3)) :: qxx,qxy,qxz,qyy,qyz,qzz,q
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: f1o3=1.d0/3.d0,inbd=1.d0,hlf=0.5d0,fou=4.d0
!
! BL radius of the outer and inner horizons 
  a = m0*sam
  wh_out = m0 + sqrt(m0*m0 - a*a)
  wh_in  = m0 - sqrt(m0*m0 - a*a)
  ! BL radius inside which junk initial data will be filled 
  if (r_junk .lt. 0.d0) then
    w0 = 0.5d0*(wh_in + wh_out)
    if (abs(sam) .gt. 0.96d0) w0 = 0.6d0*m0
  else
     if (r_junk .gt. wh_out) then
        write(*,*) 'r_junk = ',r_junk
        write(*,*) 'Horizon radius = ',wh_out
        write(*,*) 'r_junk must be smaller than the horizon radius'
        stop
     end if
     w0 = r_junk
  end if
!
! Input translation
!
  imin = lbound(gxxm1,1)
  jmin = lbound(gxxm1,2)
  kmin = lbound(gxxm1,3)
  imax = ubound(gxxm1,1)
  jmax = ubound(gxxm1,2)
  kmax = ubound(gxxm1,3)

!
!-- flat metric
!

  eta = zeo
  eta(0,0) = -one
  eta(1,1) =  one
  eta(2,2) =  one
  eta(3,3) =  one

!
! go to each gridpoint...
!
  do k = kmin, kmax
   do j = jmin, jmax
    do i = imin, imax

!
!-- the general form of the flat-space radius when boosted in xy-plane
!-- r^2=\gamma^2(v_1 x+v_2 y)^2 + x^2 + y^2 +z^2
!
    xl = X(i,j,k)-x0
    yl = Y(i,j,k)-y0
    zl = Z(i,j,k)-z0
    r2 = xl*xl + yl*yl + zl*zl

!-- W^2 = 1/2 *[ r^2 - a^2 + \sqrt(( r^2 - a^2 )^2 + 4 a^2 z^2 )] (sqrure of BL radius)

    w2 = r2 - a * a
    w2 = hlf *( w2 + sqrt( w2 * w2 + fou * a * a * zl * zl )) + m0*1.d-16
     W = sqrt(w2) 

    tmp = two * w2 + a * a - r2
    tmq =       w2 + a * a

     ! \partial_i W:
     wx =   W * xl / tmp
     wy =   W * yl / tmp
     wz = tmq * zl /(tmp * W)

     ! Modify W when W<w0: W -> w0-m0/lambda+m0/lambda*exp[lambda*(W-w0)]
     if (W .lt. w0) then
	tmr = exp(lambda*(W-w0))
        W = w0 - m0/lambda + m0/lambda*tmr
        w2 = W*W
	wx = wx*tmr
	wy = wy*tmr
	wz = wz*tmr
	tmp = two * w2 + a * a - r2
	tmq =       w2 + a * a
     end if

!-- H = M W^3 /( w^4 + a^2 z^2 )

    hh = m0 * W*w2 /( w2*w2 + a * a * zl * zl )

!
!-- the null vector
!
     el(0) = one
     el(1) = ( W * xl + a * yl )/ tmq
     el(2) = ( W * yl - a * xl )/ tmq
     el(3) = zl / W
     if (zl .eq. 0.d0 .and. r2 .lt. a*a) then 
        el(3) = sqrt(1.d0 - r2/(a*a))
     end if

!
!-- \partial_\mu H
!

     tmp = w2 * w2 + a * a * zl * zl
     tmr = thr * a * a * zl * zl - w2 * w2

     phh(0) = zeo
     phh(1) = hh * tmr * wx /( W * tmp )
     phh(2) = hh * tmr * wy /( W * tmp )
     phh(3) = hh *(tmr * wz /  W - two * a * a * zl )/ tmp

!
!-- \partial_\mu L_\nu
!
       pel = zeo

       pel(1,1) = (     W + ( xl - two * W * el(1) )* wx )/ tmq
       pel(2,1) = (   a + ( xl - two * W * el(1) )* wy )/ tmq
       pel(3,1) =           ( xl - two * W * el(1) )* wz  / tmq

       pel(1,2) = ( - a + ( yl - two * W * el(2) )* wx )/ tmq
       pel(2,2) = (     W + ( yl - two * W * el(2) )* wy )/ tmq
       pel(3,2) =           ( yl - two * W * el(2) )* wz  / tmq

       pel(1,3) =      - el(3) * wx / W
       pel(2,3) =      - el(3) * wy / W
       pel(3,3) = (one - el(3) * wz)/ W

!
!-- the lapse and shift
!--  \alpha = 1/sqrt{1+2HL_t^2}
!-- \beta^i = 2 H \alpha^2 L_t L^i
!
! Note that there are more than one BHs, the functions will be added for each BH. 
! These functions must be initialized to 0 before adding to the first BH.
!
               aa = one/sqrt( one + two*hh*el(0)*el(0) )
     alpm1(i,j,k) = alpm1(i,j,k) + (aa-one)

        Bx(i,j,k) = Bx(i,j,k) + two*aa*aa*hh*el(0)*el(1)
        By(i,j,k) = By(i,j,k) + two*aa*aa*hh*el(0)*el(2)
        Bz(i,j,k) = Bz(i,j,k) + two*aa*aa*hh*el(0)*el(3)

!-- the 3-metric - delta_ij
!

     gxxm1(i,j,k) = gxxm1(i,j,k) + two*hh*el(1)*el(1) 
     gyym1(i,j,k) = gyym1(i,j,k) + two*hh*el(2)*el(2)
     gzzm1(i,j,k) = gzzm1(i,j,k) + two*hh*el(3)*el(3) 
     gxy(i,j,k) = gxy(i,j,k) + two*hh*el(1)*el(2) 
     gyz(i,j,k) = gyz(i,j,k) + two*hh*el(2)*el(3)
     gxz(i,j,k) = gxz(i,j,k) + two*hh*el(3)*el(1)

!-- K_{ij} = - ( L_iL_j@_tH + HL_i@_tL_j + HL_j@_tL_i )/\alpha
!--          + \alpha( L_tL_i@_jH + L_tL_j@_iH + HL_i@_jL_t + HL_j@_iL_t
!--                  + HL_t@_iL_j + HL_t@_jL_i)
!--          + 2\alpha HL_tL^k( L_iL_j@_kH + HL_i@_kL_j +HL_j@_kL_i )
!--       store in Aij
!
     k1 = - ( el(1)*el(1)*phh(0) + hh*el(1)*pel(0,1) + hh*el(1)*pel(0,1) )/aa
     k2 = aa*( el(0)*el(1)*phh(1) + el(0)*el(1)*phh(1) + hh*el(1)*pel(1,0) + &
               hh*el(1)*pel(1,0) + hh*el(0)*( pel(1,1) +  pel(1,1) ) )
     k3 = two*aa*hh*el(0)*( &
          el(1)*el(1)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(1)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) + &
          hh*el(1)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) )

     Axx(i,j,k) = Axx(i,j,k) + k1 + k2 + k3

     k1 = - ( el(1)*el(2)*phh(0) + hh*el(1)*pel(0,2) + hh*el(2)*pel(0,1) )/aa
     k2 = aa*( el(0)*el(1)*phh(2) + el(0)*el(2)*phh(1) + hh*el(1)*pel(2,0) + &
               hh*el(2)*pel(1,0) + hh*el(0)*( pel(1,2) +  pel(2,1) ) )
     k3 = two*aa*hh*el(0)*( &
          el(1)*el(2)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(1)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) + &
          hh*el(2)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) )

     Axy(i,j,k) = Axy(i,j,k) + k1 + k2 + k3

     k1 = - ( el(1)*el(3)*phh(0) + hh*el(1)*pel(0,3) + hh*el(3)*pel(0,1) )/aa
     k2 = aa*( el(0)*el(1)*phh(3) + el(0)*el(3)*phh(1) + hh*el(1)*pel(3,0) + &
               hh*el(3)*pel(1,0) + hh*el(0)*( pel(1,3) +  pel(3,1) ) )
     k3 = two*aa*hh*el(0)*( &
          el(1)*el(3)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(1)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) + &
          hh*el(3)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) )

     Axz(i,j,k) = Axz(i,j,k) + k1 + k2 + k3

     k1 = - ( el(2)*el(2)*phh(0) + hh*el(2)*pel(0,2) + hh*el(2)*pel(0,2) )/aa
     k2 = aa*( el(0)*el(2)*phh(2) + el(0)*el(2)*phh(2) + hh*el(2)*pel(2,0) + &
               hh*el(2)*pel(2,0) + hh*el(0)*( pel(2,2) +  pel(2,2) ) )
     k3 = two*aa*hh*el(0)*( &
          el(2)*el(2)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(2)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) + &
          hh*el(2)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) )

     Ayy(i,j,k) = Ayy(i,j,k) + k1 + k2 + k3

     k1 = - ( el(2)*el(3)*phh(0) + hh*el(2)*pel(0,3) + hh*el(3)*pel(0,2) )/aa
     k2 = aa*( el(0)*el(2)*phh(3) + el(0)*el(3)*phh(2) + hh*el(2)*pel(3,0) + &
               hh*el(3)*pel(2,0) + hh*el(0)*( pel(2,3) +  pel(3,2) ) )
     k3 = two*aa*hh*el(0)*( &
          el(2)*el(3)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(2)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) + &
          hh*el(3)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) )

     Ayz(i,j,k) = Ayz(i,j,k) + k1 + k2 + k3

     k1 = - ( el(3)*el(3)*phh(0) + hh*el(3)*pel(0,3) + hh*el(3)*pel(0,3) )/aa
     k2 = aa*( el(0)*el(3)*phh(3) + el(0)*el(3)*phh(3) + hh*el(3)*pel(3,0) + &
               hh*el(3)*pel(3,0) + hh*el(0)*( pel(3,3) +  pel(3,3) ) )
     k3 = two*aa*hh*el(0)*( &
          el(3)*el(3)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(3)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) + &
          hh*el(3)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) )

     Azz(i,j,k) = Azz(i,j,k) + k1 + k2 + k3

    end do
   end do
  end do

  return

end subroutine ks_initial_metric_plus_junk_one_BH

!-----------------------------------------------------------------------------
!
! Find non-boosted rotating Kerr-Schild metric initial data (with junk)
! This is for the case when we're using the Boyer-Lindquist radius
!
!-----------------------------------------------------------------------------
  subroutine ks_initial_metric_bl_junk_one_BH(ex, x, y, z,             &
       gxxm1, gxy, gxz, gyym1, gyz, gzzm1,              &
       Axx, Axy, Axz, Ayy, Ayz, Azz, alpm1, Bx,By,Bz, &
       x0,y0,z0,m0, sam, r_junk, lambda)
  implicit none

!
! Input parameters:
!
  integer, dimension(3), intent(in)                  :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: x,y,z
  real*8,  intent(in)                                :: sam,m0, x0,y0,z0
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxxm1,gxy,gxz,gyym1,gyz,gzzm1
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: alpm1,Bx,By,Bz
!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: k1,k2,k3,aa,hh,r2,w2,W,wx,wy,wz,tmp,tmq,tmr
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: f1o3=1.d0/3.d0,inbd=1.d0,hlf=0.5d0,fou=4.d0

  real*8 :: grr, grph, gthth, gphph, Br, Bph
  real*8 :: Arr, Arth, Arph, Athth, Athph, Aphph
  real*8 :: Sigma, r, sin2th, sinth, cos2th, costh
  real*8 :: Lx_r, Ly_r, Lz_r, Lx_th, Ly_th, Lz_th, Lx_ph, Ly_ph, Lz_ph
  real*8 :: Lr_x, Lth_x, Lph_x, Lr_y, Lth_y, Lph_y, Lr_z, Lth_z, Lph_z
  real*8 :: r_l, x_l, y_l, z_l
  real*8 :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl
  real*8 :: Kxxl,Kxyl,Kxzl,Kyyl,Kyzl,Kzzl
  real*8  :: r0, rh_out, rh_in, r_junk, lambda,a, Sigma0, f1oalpha0
  real*8 :: aaa,bbb,grr0,grph_r,gthth_r2,gphph_r2,dgrr0,dgrph_r,dgthth_r2,dgphph_r2
  real*8 :: Arr0,dArr0,Arth_r,dArth_r,Athth_r2,dAthth_r2,Arph_r,dArph_r
  real*8 :: Athph_r2,dAthph_r2,Aphph_r2,dAphph_r2
  real*8 :: psi40,alp0,dalp0,grr0_pl,alp0_pl,dgrr0_pl,dalp0_pl

!
! BL radius of the outer and inner horizons
  a = m0*sam
  rh_out = m0 + sqrt(m0*m0 - a*a)
  rh_in  = m0 - sqrt(m0*m0 - a*a)
  ! BL radius inside which junk initial data will be filled
  if (r_junk .lt. 0.d0) then
     r0 = 0.5d0*(rh_in + rh_out)
     if (abs(sam) .gt. 0.96d0) r0 = 0.6d0*m0
  else
     if (r_junk .gt. rh_out) then
        write(*,*) 'r_junk = ',r_junk
        write(*,*) 'Horizon radius = ',rh_out
        write(*,*) 'r_junk must be smaller than the horizon radius'
        stop
     end if
     r0 = r_junk
  end if

! Input translation
!
  imin = lbound(gxxm1,1)
  jmin = lbound(gxxm1,2)
  kmin = lbound(gxxm1,3)
  imax = ubound(gxxm1,1)
  jmax = ubound(gxxm1,2)
  kmax = ubound(gxxm1,3)

!
! go to each gridpoint...
!
  do k = kmin, kmax
   do j = jmin, jmax
    do i = imin, imax

    x_l = x(i,j,k) - x0
    y_l = y(i,j,k) - y0
    z_l = z(i,j,k) - z0
    r2 = x_l*x_l + y_l*y_l + z_l*z_l
    r = sqrt(r2)
    costh = z_l/r
    cos2th = costh**2
    sin2th = one - cos2th
    sinth = sqrt(sin2th)
    Sigma = r2 + a*a*cos2th
!
!-- the lapse and shift (nonzero component)
!
     aa = 1.d0/sqrt(one+two*m0*r/Sigma)
      Br = two*m0*r/(Sigma + two*m0*r)
!
!-- the 3-metric (nonzero components)
!
     grr0_pl = 1.d0+2.d0*m0/r0
     dgrr0_pl = -2.d0*m0/r0/r0
     aaa = 0.5d0*dgrr0_pl/r0 
     psi40 = grr0_pl - aaa*r0*r0
     alp0_pl = 1.d0/sqrt(1.d0+2.d0*m0/r0)
     dalp0_pl = m0/r0/r0 * alp0_pl**3
     aaa = 0.5d0*dalp0_pl/r0
     alp0_pl = alp0_pl - aaa*r0*r0
     if (r .ge. r0) then 
        grr = one+two*m0*r/Sigma
        grph = -a*sin2th*grr 
        gthth = Sigma
        gphph = sin2th*(Sigma + a*a*sin2th*grr)
     else
	Sigma0 = r0*r0+a*a*cos2th
	grr0 = one+two*m0*r0/Sigma0
        grph_r = -a*sin2th*grr0/r0
        gthth_r2 = Sigma0/r0/r0
        gphph_r2 = sin2th*(Sigma0 + a*a*sin2th*grr0)/r0/r0
	dgrr0 = -2.d0*m0*(r0*r0-a*a*cos2th)/Sigma0/Sigma0
	dgrph_r = -a*sin2th*dgrr0/r0 + a*sin2th*grr0/r0/r0
	dgthth_r2 = -2.d0*a*a*cos2th/r0**3
	dgphph_r2 = sin2th*(2.d0*r0 + a*a*sin2th*dgrr0)/r0/r0 - 2.d0*gphph_r2/r0
        call compute_ab_bbh(r0,grr0,dgrr0,psi40,aaa,bbb)
	grr = psi40 + aaa*r2 + bbb*r2*r
	call compute_ab_bbh(r0,grph_r,dgrph_r,0.d0,aaa,bbb)
	grph = r*(aaa*r2 + bbb*r2*r)
        call compute_ab_bbh(r0,gthth_r2,dgthth_r2,psi40,aaa,bbb)
	gthth = r2*(psi40 + aaa*r2 + bbb*r2*r)
	call compute_ab_bbh(r0,gphph_r2,dgphph_r2,psi40*sin2th,aaa,bbb)
	gphph = r2*(psi40*sin2th + aaa*r2 + bbb*r2*r)
	alp0 = 1.d0/sqrt(one+two*m0*r0/Sigma0)
	dalp0 = m0*(r0*r0-a*a*cos2th)/Sigma0/Sigma0*alp0**3
        call compute_ab_bbh(r0,alp0,dalp0,alp0_pl,aaa,bbb)
        aa = alp0_pl + aaa*r2 + bbb*r2*r
     end if

     alpm1(i,j,k) = alpm1(i,j,k) + (aa-1.d0)
     if (alpm1(i,j,k)+1.d0 .lt. 0.d0) then
        write(*,*) 'X, Y, Z = ',X(i,j,k),Y(i,j,k),Z(i,j,k)
        write(*,*) 'lapse = ',alpm1(i,j,k)+1.d0, ' < 0!'
        stop
     end if

!
!-- Store K_{ij} in A_{ij} 
!--         
!
     if (r .ge. r0) then 
        Arr = two*m0*(two*Sigma - 4.d0*r2)/(two*aa)
        Arr = Arr/(Sigma**2*(Sigma + two*m0*r)**2)
        Arr = Arr*(3.d0*m0*r*Sigma + two*m0*m0*r2 + Sigma**2)

        Arth = 4.d0*m0*r*a*a*sinth*costh
        Arth = Arth / (Sigma*(Sigma+two*m0*r))/(two*aa)
        
        Arph = -two*m0*a*sin2th*(one - two*r2/Sigma)/Sigma/(two*aa)

        Athth = 4.d0*m0*r2/(Sigma + two*m0*r)/(two*aa)

        Athph = -4.d0*m0*r*a*a*a*sin2th*sinth*costh
        Athph = Athph/(Sigma*(Sigma+two*m0*r))/(two*aa)

        Aphph = 4.d0*m0*r*sin2th/(Sigma + two*m0*r)
        Aphph = Aphph*(r + m0*a*a*sin2th*((Sigma-two*r2)/Sigma**2))/(two*aa)
     else
	Sigma0 = r0*r0 + a*a*cos2th
	f1oalpha0 = sqrt(1.d0+2.d0*m0*r0/Sigma0)

	Arr0 = m0*(two*Sigma0 - 4.d0*r0*r0)*f1oalpha0
        Arr0 = Arr0/(Sigma0**2*(Sigma0 + two*m0*r0)**2)
        Arr0 = Arr0*(3.d0*m0*r0*Sigma0 + two*m0*m0*r0*r0 + Sigma0**2)
	dArr0 = -m0*(r0*r0 - a*a*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*m0*r0/Sigma0) & 
		- 2.d0*r0/(Sigma0 - 2.d0*r0*r0) + & 
		(Sigma0*(4.d0*r0+3.d0*m0)+4.d0*m0*m0*r0+6.d0*m0*r0*r0)/(Sigma0**2 + 3.d0*m0*r0*Sigma0 + 2.d0*m0*m0*r0*r0) & 
		- 4.d0*r0/Sigma0 - 4.d0*(r0+m0)/(Sigma0+2.d0*m0*r0)
	dArr0 = Arr0*dArr0
        call compute_ab_bbh(r0,Arr0,dArr0,0.d0,aaa,bbb)
	Arr = aaa*r2 + bbb*r2*r

        Arth_r = 2.d0*m0*a*a*sinth*costh/(Sigma0*(Sigma0+two*m0*r0))*f1oalpha0
	dArth_r = -m0*(r0*r0 - a*a*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*m0*r0/Sigma0) & 
		  - 2.d0*r0/Sigma0 - 2.d0*(r0+m0)/(Sigma0+2.d0*m0*r0)
	dArth_r = Arth_r*dArth_r
	call compute_ab_bbh(r0,Arth_r,dArth_r,0.d0,aaa,bbb)
	Arth = r*(aaa*r2 + bbb*r2*r)

	Arph_r = -m0*a*sin2th*(one - two*r0*r0/Sigma0)/r0/Sigma0*f1oalpha0
	dArph_r = -m0*(r0*r0 - a*a*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*m0*r0/Sigma0) &
		  - 4.d0*a*a*cos2th*r0/Sigma0/Sigma0/(1.d0-2.d0*r0*r0/Sigma0) & 
		  - 1.d0/r0 - 2.d0*r0/Sigma0
	dArph_r = dArph_r * Arph_r
	call compute_ab_bbh(r0,Arph_r,dArph_r,0.d0,aaa,bbb)
	Arph = r*(aaa*r2 + bbb*r2*r)

        Athth_r2 = 2.d0*m0/(Sigma0 + two*m0*r0)*f1oalpha0
	dAthth_r2 = -(r0*r0 - a*a*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*m0*r0/Sigma0) &
		  	- 2.d0*(r0+m0)/(Sigma0+2.d0*m0*r0)
	dAthth_r2 = dAthth_r2 *Athth_r2
	call compute_ab_bbh(r0,Athth_r2,dAthth_r2,0.d0,aaa,bbb)
	Athth = r2*(aaa*r2 + bbb*r2*r)

        Athph_r2 = -2.d0*m0*a*a*a*sin2th*sinth*costh/r0
        Athph_r2 = Athph_r2/(Sigma0*(Sigma0+two*m0*r0))*f1oalpha0
	dAthph_r2 = -(r0*r0 - a*a*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*m0*r0/Sigma0) &
		    - 1.d0/r0 - 2.d0*r0/Sigma0 - 2.d0*(r0+m0)/(Sigma0+2.d0*m0*r0)
	dAthph_r2 = dAthph_r2 *Athph_r2
	call compute_ab_bbh(r0,Athph_r2,dAthph_r2,0.d0,aaa,bbb)
	Athph = r2*(aaa*r2 + bbb*r2*r)

        Aphph_r2 = 2.d0*m0*sin2th/(Sigma0 + two*m0*r0)/r0
        Aphph_r2 = Aphph_r2*(r0 + m0*a*a*sin2th*((Sigma0-two*r0*r0)/Sigma0**2))*f1oalpha0
	dAphph_r2 = -(r0*r0 - a*a*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*m0*r0/Sigma0) &
		    - 1.d0/r0 - 2.d0*(r0+m0)/(Sigma0+2.d0*m0*r0) + & 
		    (1.d0+4.d0*a*a*sin2th*m0*r0*(r0*r0-a*a*cos2th)/Sigma0**3 & 
		     - 2.d0*a*a*sin2th*r0/Sigma0**2) & 
		    / (r0 + m0*a*a*sin2th*((Sigma0-two*m0*r0*r0)/Sigma0**2))
	dAphph_r2 = dAphph_r2 * Aphph_r2
	call compute_ab_bbh(r0,Aphph_r2,dAphph_r2,0.d0,aaa,bbb)
	Aphph = r2*(aaa*r2 + bbb*r2*r)

     end if


     !#############################################################################
     !Now convert everything to Cartesian
     !#############################################################################

       r_l = r
       r2 = r_l * r_l

       x_l = x(i,j,k)
       y_l = y(i,j,k)
       z_l = z(i,j,k)

       Lx_r = x_l/r_l  ! = sin(th)cos(ph)
       Ly_r = y_l/r_l  ! = sin(th)sin(ph)
       Lz_r = z_l/r_l  ! = cos(ph)

       Lx_th = z_l*x_l/sqrt(r2 - z_l**2)
       Ly_th = z_l*y_l/sqrt(r2 - z_l**2)
       Lz_th = -sqrt(r2 - z_l**2)

       Lx_ph = -y_l
       Ly_ph = x_l
       Lz_ph = 0.d0

       Lr_x = Lx_r
       Lr_y = Ly_r
       Lr_z = Lz_r

       Lth_x = z_l*x_l/sqrt(r2-z_l**2)/(r_l*r_l)
       Lth_y = z_l*y_l/sqrt(r2-z_l**2)/(r_l*r_l)
       Lth_z = -sqrt(r2-z_l**2)/(r_l*r_l)

       Lph_x = -y_l/(x_l**2 + y_l**2)
       Lph_y = x_l/(x_l**2 + y_l**2)
       Lph_z = 0.d0

       Bx(i,j,k) = Lx_r*Br
       By(i,j,k) = Ly_r*Br
       Bz(i,j,k) = Lz_r*Br

       ! Transform the spatial metric
       gxxl = Lr_x**2 * grr + &
                    2.d0*Lr_x*Lph_x * grph + &
                    Lth_x**2 * gthth + &
                    Lph_x**2 * gphph

       gxyl = Lr_x*Lr_y * grr + &
                    (Lr_x*Lph_y + Lr_y*Lph_x) * grph + &
                    Lth_x*Lth_y * gthth + &
                    Lph_x*Lph_y * gphph

       gxzl = Lr_x*Lr_z * grr + &
                    (Lr_x*Lph_z + Lr_z*Lph_x) * grph + &
                    Lth_x*Lth_z * gthth + &
                    Lph_x*Lph_z * gphph

       gyyl = Lr_y**2 * grr + &
                    2.d0*Lr_y*Lph_y * grph + &
                    Lth_y**2 * gthth + &
                    Lph_y**2 * gphph

       gyzl = Lr_y*Lr_z * grr + &
                    (Lr_y*Lph_z + Lr_z*Lph_y) * grph + &
                    Lth_y*Lth_z * gthth + &
                    Lph_y*Lph_z * gphph

       gzzl = Lr_z**2 * grr + &
                    2.d0*Lr_z*Lph_z * grph + &
                    Lth_z**2 * gthth + &
                    Lph_z**2 * gphph

       ! Transform the extrinsic curvature
       Kxxl = Lr_x**2 * Arr + &
		    2.d0*Lr_x*Lth_x * Arth + &
                    2.d0*Lr_x*Lph_x * Arph + &
                    Lth_x**2 * Athth + &
		    2.d0*Lth_x*Lph_x * Athph + &
                    Lph_x**2 * Aphph

       Kxyl = Lr_x*Lr_y * Arr + &
		    (Lr_x*Lth_y + Lr_y*Lth_x) * Arth + &
                    (Lr_x*Lph_y + Lr_y*Lph_x) * Arph + &
                    Lth_x*Lth_y * Athth + &
		    (Lth_x*Lph_y + Lth_y*Lph_x) * Athph + &
                    Lph_x*Lph_y * Aphph

       Kxzl = Lr_x*Lr_z * Arr + &
		    (Lr_x*Lth_z + Lr_z*Lth_x) * Arth + &
                    (Lr_x*Lph_z + Lr_z*Lph_x) * Arph + &
                    Lth_x*Lth_z * Athth + &
		    (Lth_x*Lph_z + Lth_z*Lph_x) * Athph + &
                    Lph_x*Lph_z * Aphph

       Kyyl = Lr_y**2 * Arr + &
		    2.d0*Lr_y*Lth_y * Arth + &
                    2.d0*Lr_y*Lph_y * Arph + &
                    Lth_y**2 * Athth + &
		    2.d0*Lth_y*Lph_y * Athph + &
                    Lph_y**2 * Aphph

       Kyzl = Lr_y*Lr_z * Arr + &
		    (Lr_y*Lth_z + Lr_z*Lth_y) * Arth + &
                    (Lr_y*Lph_z + Lr_z*Lph_y) * Arph + &
                    Lth_y*Lth_z * Athth + &
		    (Lth_y*Lph_z + Lth_z*Lph_y) * Athph + &
                    Lph_y*Lph_z * Aphph

       Kzzl = Lr_z**2 * Arr + &
		    2.d0*Lr_z*Lth_z * Arth + &
                    2.d0*Lr_z*Lph_z * Arph + &
                    Lth_z**2 * Athth + &
		    2.d0*Lth_z*Lph_z * Athph + &
                    Lph_z**2 * Aphph

       gxxm1(i,j,k) = gxxm1(i,j,k) +  (gxxl - 1.d0)
       gxy(i,j,k) = gxy(i,j,k) + gxyl
       gxz(i,j,k) = gxz(i,j,k) + gxzl
       gyym1(i,j,k) = gyym1(i,j,k) + (gyyl - 1.d0)
       gyz(i,j,k) = gyz(i,j,k) + gyzl
       gzzm1(i,j,k) = gzzm1(i,j,k) + (gzzl - 1.d0)

       Axx(i,j,k) = Axx(i,j,k) + Kxxl
       Axy(i,j,k) = Axy(i,j,k) + Kxyl
       Axz(i,j,k) = Axz(i,j,k) + Kxzl
       Ayy(i,j,k) = Ayy(i,j,k) + Kyyl
       Ayz(i,j,k) = Ayz(i,j,k) + Kyzl
       Azz(i,j,k) = Azz(i,j,k) + Kzzl

    end do
   end do
  end do

  return

end subroutine ks_initial_metric_bl_junk_one_BH

!-----------------------------------------------------------------------------
!
! Find non-boosted rotating Kerr-Schild "puncture" initial data 
! Note: values of lapm1, beta^i, K_ij and (g_ij - delta_ij) are added 
!       each time this subroutine is called. They must be initialized 
!       to zero before calling this subroutine for the first time.
!
!-----------------------------------------------------------------------------
  subroutine kerr_initial_metric_puncture_one_BH(ex, x, y, z, x0,y0,z0,   &
       gxxm1, gxy, gxz, gyym1, gyz, gzzm1,             &
       Axx, Axy, Axz, Ayy, Ayz, Azz, alpm1, Bx,By,Bz, &
       sam,m0)
  implicit none

!
! Input parameters:
!
  integer, dimension(3), intent(in)                  :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: x,y,z
  real*8,  intent(in)                                :: sam,m0,x0,y0,z0
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxxm1,gxy,gxz,gyym1,gyz,gzzm1
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: alpm1,Bx,By,Bz
!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: k1,k2,k3,aa,hh,r2,w2,W,wx,wy,wz,tmp,tmq,tmr
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: inbd=1.d0,hlf=0.5d0,fou=4.d0

  real*8 :: grr, grph, gthth, gphph, alp,Bph
  real*8 :: Krr, Krth,Krph,Kthth,Kthph,Kphph
  real*8 :: Sigma, r, sin2th, sinth, cos2th, costh
  real*8 :: Lx_r, Ly_r, Lz_r, Lx_th, Ly_th, Lz_th, Lx_ph, Ly_ph, Lz_ph
  real*8 :: Lr_x, Lth_x, Lph_x, Lr_y, Lth_y, Lph_y, Lr_z, Lth_z, Lph_z
  real*8 :: r_l, x_l, y_l, z_l, Delta, dtdr, dpdr
  real*8 :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl
  real*8 :: Kxxl,Kxyl,Kxzl,Kyyl,Kyzl,Kzzl
  real*8  :: rh_out, rh_in, a, riso, dr_driso, dr_driso_alp

!
! BL radius of the outer and inner horizons
!
  a = m0*sam
  rh_out = m0 + sqrt(m0*m0 - a*a)
  rh_in  = m0 - sqrt(m0*m0 - a*a)

! Input translation
!
  imin = lbound(gxxm1,1)
  jmin = lbound(gxxm1,2)
  kmin = lbound(gxxm1,3)
  imax = ubound(gxxm1,1)
  jmax = ubound(gxxm1,2)
  kmax = ubound(gxxm1,3)

!
! go to each gridpoint...
!
  do k = kmin, kmax
   do j = jmin, jmax
    do i = imin, imax
 
    x_l = x(i,j,k) - x0
    y_l = y(i,j,k) - y0
    z_l = z(i,j,k) - z0
    riso = sqrt(x_l*x_l + y_l*y_l + z_l*z_l)
    r = riso*(1.d0 + 0.25d0*rh_out/riso)**2
    dr_driso = 1.d0 - (0.25d0*rh_out/riso)**2
    costh = z_l/riso
    cos2th = costh**2
    sin2th = one - cos2th
    sinth = sqrt(sin2th)
    r2 = r * r
    Sigma = r2 + a*a*cos2th
    Delta = r2 - 2.d0*m0*r + a*a
    aa = (r2+a*a)**2 - a*a*Delta*sin2th
!
!-- the lapse and shift (nonzero component)
!
     alp = sqrt(Delta*Sigma/aa)
     alpm1(i,j,k) = alpm1(i,j,k) + (alp - 1.d0)
      Bph = -2.d0*m0*a*r/aa

! Compute (dr/driso)/alpha
     !!if (riso .ge. 0.25d0*rh_out) then 
     !!   dr_driso_alp = 0.125d0*(sqrt(r-rh_out)+sqrt(r))*(4.d0*riso+rh_out)/riso/riso &
     !!   		* sqrt(aa/Sigma/(r-rh_in))
     !!else
     !!   dr_driso_alp = -0.125d0*rh_out/(sqrt(r-rh_out)+sqrt(r))* & 
     !!   		(4.d0*riso+rh_out)/riso/riso * sqrt(aa/Sigma/(r-rh_in))
     !!end if
    
     dr_driso_alp = sqrt(aa*riso/Sigma/(r-rh_in))*(riso+0.25d0*rh_out)/(riso*riso)

!
!-- the 3-metric (nonzero components)
!
!     grr = Sigma/Delta
     if (riso .ge. 0.25d0*rh_out) then 
        grr = Sigma/64.d0/(r-rh_in)*( (sqrt(r-rh_out)+sqrt(r))*(4*riso+rh_out)/riso/riso )**2
     else
         grr = Sigma/64.d0/(r-rh_in)*( rh_out/(sqrt(r-rh_out)+sqrt(r))*(4*riso+rh_out)/riso/riso )**2
     end if
     grph = 0.d0
     gthth = Sigma
     gphph = aa*sin2th/Sigma
!
!--  K_{ij} in "isotropic" Boyer-Lindquist coordinates
!--         
!
     Krr = 0.d0
     Krth = 0.d0
     Krph = (m0*a*sin2th*(-a**4 + 2.d0*a*a*r2 + 3.d0*r2*r2 + &
            a*a*(a - r)*(a + r)*sin2th))/Sigma/aa * dr_driso_alp
     Kthth = 0.d0
     Kthph = -2.d0*a*a*a*m0*r*costh*sinth*sin2th/Sigma * sqrt(Delta/Sigma/aa)
     if (riso .lt. 0.25d0*rh_out) Kthph = -Kthph
     Kphph = 0.d0

     !#############################################################################
     !Now convert everything to Cartesian
     !#############################################################################

       r_l = riso
       r2 = r_l * r_l

       Lx_r = x_l/r_l  ! = sin(th)cos(ph)
       Ly_r = y_l/r_l  ! = sin(th)sin(ph)
       Lz_r = z_l/r_l  ! = cos(ph)

       Lx_th = z_l*x_l/sqrt(r2 - z_l**2)
       Ly_th = z_l*y_l/sqrt(r2 - z_l**2)
       Lz_th = -sqrt(r2 - z_l**2)

       Lx_ph = -y_l
       Ly_ph = x_l
       Lz_ph = 0.d0

       Lr_x = Lx_r
       Lr_y = Ly_r
       Lr_z = Lz_r

       Lth_x = z_l*x_l/sqrt(r2-z_l**2)/(r_l*r_l)
       Lth_y = z_l*y_l/sqrt(r2-z_l**2)/(r_l*r_l)
       Lth_z = -sqrt(r2-z_l**2)/(r_l*r_l)

       Lph_x = -y_l/(x_l**2 + y_l**2)
       Lph_y = x_l/(x_l**2 + y_l**2)
       Lph_z = 0.d0

       Bx(i,j,k) = Bx(i,j,k) + Lx_ph*Bph
       By(i,j,k) = By(i,j,k) + Ly_ph*Bph

       ! Transform the spatial metric
       gxxl = Lr_x**2 * grr + &
                    2.d0*Lr_x*Lph_x * grph + &
                    Lth_x**2 * gthth + &
                    Lph_x**2 * gphph

       gxyl = Lr_x*Lr_y * grr + &
                    (Lr_x*Lph_y + Lr_y*Lph_x) * grph + &
                    Lth_x*Lth_y * gthth + &
                    Lph_x*Lph_y * gphph

       gxzl = Lr_x*Lr_z * grr + &
                    (Lr_x*Lph_z + Lr_z*Lph_x) * grph + &
                    Lth_x*Lth_z * gthth + &
                    Lph_x*Lph_z * gphph

       gyyl = Lr_y**2 * grr + &
                    2.d0*Lr_y*Lph_y * grph + &
                    Lth_y**2 * gthth + &
                    Lph_y**2 * gphph

       gyzl = Lr_y*Lr_z * grr + &
                    (Lr_y*Lph_z + Lr_z*Lph_y) * grph + &
                    Lth_y*Lth_z * gthth + &
                    Lph_y*Lph_z * gphph

       gzzl = Lr_z**2 * grr + &
                    2.d0*Lr_z*Lph_z * grph + &
                    Lth_z**2 * gthth + &
                    Lph_z**2 * gphph

       ! Transform the extrinsic curvature
       Kxxl = Lr_x**2 * Krr + &
		    2.d0*Lr_x*Lth_x * Krth + &
                    2.d0*Lr_x*Lph_x * Krph + &
                    Lth_x**2 * Kthth + &
		    2.d0*Lth_x*Lph_x * Kthph + &
                    Lph_x**2 * Kphph

       Kxyl = Lr_x*Lr_y * Krr + &
		    (Lr_x*Lth_y + Lr_y*Lth_x) * Krth + &
                    (Lr_x*Lph_y + Lr_y*Lph_x) * Krph + &
                    Lth_x*Lth_y * Kthth + &
		    (Lth_x*Lph_y + Lth_y*Lph_x) * Kthph + &
                    Lph_x*Lph_y * Kphph

       Kxzl = Lr_x*Lr_z * Krr + &
		    (Lr_x*Lth_z + Lr_z*Lth_x) * Krth + &
                    (Lr_x*Lph_z + Lr_z*Lph_x) * Krph + &
                    Lth_x*Lth_z * Kthth + &
		    (Lth_x*Lph_z + Lth_z*Lph_x) * Kthph + &
                    Lph_x*Lph_z * Kphph

       Kyyl = Lr_y**2 * Krr + &
		    2.d0*Lr_y*Lth_y * Krth + &
                    2.d0*Lr_y*Lph_y * Krph + &
                    Lth_y**2 * Kthth + &
		    2.d0*Lth_y*Lph_y * Kthph + &
                    Lph_y**2 * Kphph

       Kyzl = Lr_y*Lr_z * Krr + &
		    (Lr_y*Lth_z + Lr_z*Lth_y) * Krth + &
                    (Lr_y*Lph_z + Lr_z*Lph_y) * Krph + &
                    Lth_y*Lth_z * Kthth + &
		    (Lth_y*Lph_z + Lth_z*Lph_y) * Kthph + &
                    Lph_y*Lph_z * Kphph

       Kzzl = Lr_z**2 * Krr + &
		    2.d0*Lr_z*Lth_z * Krth + &
                    2.d0*Lr_z*Lph_z * Krph + &
                    Lth_z**2 * Kthth + &
		    2.d0*Lth_z*Lph_z * Kthph + &
                    Lph_z**2 * Kphph

       gxxm1(i,j,k) = gxxm1(i,j,k) + (gxxl-1.d0)
       gxy(i,j,k) = gxy(i,j,k) + gxyl
       gxz(i,j,k) = gxz(i,j,k) + gxzl
       gyym1(i,j,k) = gyym1(i,j,k) + (gyyl-1.d0)
       gyz(i,j,k) = gyz(i,j,k) + gyzl
       gzzm1(i,j,k) = gzzm1(i,j,k) + (gzzl-1.d0)
       Axx(i,j,k) = Axx(i,j,k) + Kxxl
       Axy(i,j,k) = Axy(i,j,k) + Kxyl
       Axz(i,j,k) = Axz(i,j,k) + Kxzl 
       Ayy(i,j,k) = Ayy(i,j,k) + Kyyl
       Ayz(i,j,k) = Ayz(i,j,k) + Kyzl
       Azz(i,j,k) = Azz(i,j,k) + Kzzl

    end do
   end do
  end do

  return

end subroutine kerr_initial_metric_puncture_one_BH

!-----------------------------------------------------------------------------
!
! Given a function f(r), find the values of a and b so that
!   g(r) = f0 + a r^2 + b r^3, with g(r0) = f(r0) and g'(r0)=f'(r0).
! Here f0 is a constant.
!
!-----------------------------------------------------------------------------
subroutine compute_ab_bbh(r0,fr0,dfr0,f0,a,b)
   implicit none
   real*8 :: f0,fr0,dfr0,a,b,r0
!
   b = dfr0/r0**2 + 2.d0*(f0-fr0)/r0**3
   a = (fr0-f0)/r0**2 - b*r0
end subroutine compute_ab_bbh
