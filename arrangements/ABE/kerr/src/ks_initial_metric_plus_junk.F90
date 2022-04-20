!-----------------------------------------------------------------------------
!
! Find non-boosted rotating Kerr-Schild metric initial data
!
!-----------------------------------------------------------------------------
  subroutine ks_initial_metric_nofisheye(ex, x, y, z, PhysR,       &
                    gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
                    Axx, Axy, Axz, Ayy, Ayz, Azz, alpha, Bx,By,Bz, &
                    sam)
  implicit none

!
! Input parameters:
!
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z, PhysR
  real*8,                    intent(in) :: sam
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: trK,alpha,Bx,By,Bz

!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: el(0:3),boo(0:3),phh(0:3),eta(0:3,0:3),pel(0:3,0:3)
  real*8  :: k1,k2,k3,aa,hh,r2,w2,W,wx,wy,wz,tmp,tmq,tmr
!  real*8,  dimension(ex(1),ex(2),ex(3)) :: qxx,qxy,qxz,qyy,qyz,qzz,q
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: f1o3=1.d0/3.d0,m0=1.d0,inbd=1.d0,hlf=0.5d0,fou=4.d0

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
    r2 = PhysR(i,j,k)
    r2 = r2 * r2

!-- W^2 = 1/2 *[ r^2 - a^2 + \sqrt(( r^2 - a^2 )^2 + 4 a^2 z^2 )]

    w2 = r2 - sam * sam
    w2 = hlf *( w2 + sqrt( w2 * w2 + fou * sam * sam * z(i,j,k) * z(i,j,k) )) 
     W = sqrt(w2) 

    tmp = two * w2 + sam * sam - r2
    tmq =       w2 + sam * sam

     wx =   W * x(i,j,k) / tmp
     wy =   W * y(i,j,k) / tmp
     wz = tmq * z(i,j,k) /(tmp * W)

!-- H = M W^3 /( w^4 + a^2 z^2 )

    hh = m0 * w2 * W /( w2 * w2 + sam * sam * z(i,j,k) * z(i,j,k) )

!
!-- the null vector
!
     el(0) = one
     el(1) = ( W * x(i,j,k) + sam * y(i,j,k) )/ tmq
     el(2) = ( W * y(i,j,k) - sam * x(i,j,k) )/ tmq
     el(3) =       z(i,j,k)               / W

!
!-- @_\mu H
!

     tmp = w2 * w2 + sam * sam * z(i,j,k) * z(i,j,k)
     tmr = thr * sam * sam * z(i,j,k) * z(i,j,k) - w2 * w2

     phh(0) = zeo
     phh(1) = hh * tmr * wx /( W * tmp )
     phh(2) = hh * tmr * wy /( W * tmp )
     phh(3) = hh *(tmr * wz /  W - two * sam * sam * z(i,j,k) )/ tmp

!
!-- @_\mu L_\nu
!
       pel = zeo

       pel(1,1) = (     W + ( x(i,j,k) - two * W * el(1) )* wx )/ tmq
       pel(2,1) = (   sam + ( x(i,j,k) - two * W * el(1) )* wy )/ tmq
       pel(3,1) =           ( x(i,j,k) - two * W * el(1) )* wz  / tmq

       pel(1,2) = ( - sam + ( y(i,j,k) - two * W * el(2) )* wx )/ tmq
       pel(2,2) = (     W + ( y(i,j,k) - two * W * el(2) )* wy )/ tmq
       pel(3,2) =           ( y(i,j,k) - two * W * el(2) )* wz  / tmq

       pel(1,3) =      - el(3) * wx / W
       pel(2,3) =      - el(3) * wy / W
       pel(3,3) = (one - el(3) * wz)/ W

!
!-- the lapse and shift
!--  \alpha = 1/sqrt{1+2HL_t^2}
!-- \beta^i = 2 H \alpha^2 L_t L^i
!
     alpha(i,j,k) = one/sqrt( one + two*hh*el(0)*el(0) )
               aa = alpha(i,j,k)

        Bx(i,j,k) = two*aa*aa*hh*el(0)*el(1)
        By(i,j,k) = two*aa*aa*hh*el(0)*el(2)
        Bz(i,j,k) = two*aa*aa*hh*el(0)*el(3)

!
!-- the 3-metric
!

     gxx(i,j,k) = two*hh*el(1)*el(1) + one
     gyy(i,j,k) = two*hh*el(2)*el(2) + one
     gzz(i,j,k) = two*hh*el(3)*el(3) + one
     gxy(i,j,k) = two*hh*el(1)*el(2)
     gyz(i,j,k) = two*hh*el(2)*el(3)
     gxz(i,j,k) = two*hh*el(3)*el(1)

!
!-- the scalar extrinsic curvature
!-- K = \alpha L_t ( 2H @_iL^i -2H @_tL_t -L_t @_tH)
!--    + 2\alpha^3 (1 + HL_t^2)L_tL^i @_i H
!--    + [4\alpha^3(1+HL_t^2) - 2\alpha]HL^i @_iL_t
!
     k1 = aa*el(0)*( two*hh*( pel(1,1) + pel(2,2) + pel(3,3) - pel(0,0) ) - &
          el(0)*phh(0) )
     k2 = two*aa*aa*aa*el(0)*( el(1)*phh(1) + el(2)*phh(2) +el(3)*phh(3) )* &
          ( one + hh*el(0)*el(0) )
     k3 = two*aa*hh*( el(1)*pel(1,0) + el(2)*pel(2,0) + el(3)*pel(3,0) )* &
          ( two*aa*aa*( one + hh*el(0)*el(0) ) - one )

     trK(i,j,k) = k1 + k2 + k3

!
!-- A_{ij} = - ( L_iL_j@_tH + HL_i@_tL_j + HL_j@_tL_i )/\alpha
!--          + \alpha( L_tL_i@_jH + L_tL_j@_iH + HL_i@_jL_t + HL_j@_iL_t
!--                  + HL_t@_iL_j + HL_t@_jL_i)
!--          + 2\alpha HL_tL^k( L_iL_j@_kH + HL_i@_kL_j +HL_j@_kL_i )
!--          - \gamma_{ij}K/3
!
     k1 = - ( el(1)*el(1)*phh(0) + hh*el(1)*pel(0,1) + hh*el(1)*pel(0,1) )/aa
     k2 = aa*( el(0)*el(1)*phh(1) + el(0)*el(1)*phh(1) + hh*el(1)*pel(1,0) + &
               hh*el(1)*pel(1,0) + hh*el(0)*( pel(1,1) +  pel(1,1) ) )
     k3 = two*aa*hh*el(0)*( &
          el(1)*el(1)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(1)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) + &
          hh*el(1)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) )

     Axx(i,j,k) = k1 + k2 + k3 - f1o3*gxx(i,j,k)*trK(i,j,k)

     k1 = - ( el(1)*el(2)*phh(0) + hh*el(1)*pel(0,2) + hh*el(2)*pel(0,1) )/aa
     k2 = aa*( el(0)*el(1)*phh(2) + el(0)*el(2)*phh(1) + hh*el(1)*pel(2,0) + &
               hh*el(2)*pel(1,0) + hh*el(0)*( pel(1,2) +  pel(2,1) ) )
     k3 = two*aa*hh*el(0)*( &
          el(1)*el(2)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(1)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) + &
          hh*el(2)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) )

     Axy(i,j,k) = k1 + k2 + k3 - f1o3*gxy(i,j,k)*trK(i,j,k)

     k1 = - ( el(1)*el(3)*phh(0) + hh*el(1)*pel(0,3) + hh*el(3)*pel(0,1) )/aa
     k2 = aa*( el(0)*el(1)*phh(3) + el(0)*el(3)*phh(1) + hh*el(1)*pel(3,0) + &
               hh*el(3)*pel(1,0) + hh*el(0)*( pel(1,3) +  pel(3,1) ) )
     k3 = two*aa*hh*el(0)*( &
          el(1)*el(3)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(1)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) + &
          hh*el(3)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) )

     Axz(i,j,k) = k1 + k2 + k3 - f1o3*gxz(i,j,k)*trK(i,j,k)

     k1 = - ( el(2)*el(2)*phh(0) + hh*el(2)*pel(0,2) + hh*el(2)*pel(0,2) )/aa
     k2 = aa*( el(0)*el(2)*phh(2) + el(0)*el(2)*phh(2) + hh*el(2)*pel(2,0) + &
               hh*el(2)*pel(2,0) + hh*el(0)*( pel(2,2) +  pel(2,2) ) )
     k3 = two*aa*hh*el(0)*( &
          el(2)*el(2)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(2)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) + &
          hh*el(2)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) )

     Ayy(i,j,k) = k1 + k2 + k3 - f1o3*gyy(i,j,k)*trK(i,j,k)

     k1 = - ( el(2)*el(3)*phh(0) + hh*el(2)*pel(0,3) + hh*el(3)*pel(0,2) )/aa
     k2 = aa*( el(0)*el(2)*phh(3) + el(0)*el(3)*phh(2) + hh*el(2)*pel(3,0) + &
               hh*el(3)*pel(2,0) + hh*el(0)*( pel(2,3) +  pel(3,2) ) )
     k3 = two*aa*hh*el(0)*( &
          el(2)*el(3)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(2)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) + &
          hh*el(3)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) )

     Ayz(i,j,k) = k1 + k2 + k3 - f1o3*gyz(i,j,k)*trK(i,j,k)

     k1 = - ( el(3)*el(3)*phh(0) + hh*el(3)*pel(0,3) + hh*el(3)*pel(0,3) )/aa
     k2 = aa*( el(0)*el(3)*phh(3) + el(0)*el(3)*phh(3) + hh*el(3)*pel(3,0) + &
               hh*el(3)*pel(3,0) + hh*el(0)*( pel(3,3) +  pel(3,3) ) )
     k3 = two*aa*hh*el(0)*( &
          el(3)*el(3)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(3)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) + &
          hh*el(3)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) )

     Azz(i,j,k) = k1 + k2 + k3 - f1o3*gzz(i,j,k)*trK(i,j,k)

!~~~~~> test
!
!     q(i,j,k)=two * aa * hh * ( one + aa * aa * hh ) / r
!
!qxx(i,j,k)=two*aa*hh*(one-(two+hh)*el(1)*el(1))/r - f1o3*gxx(i,j,k)*q(i,j,k)
!qyy(i,j,k)=two*aa*hh*(one-(two+hh)*el(2)*el(2))/r - f1o3*gyy(i,j,k)*q(i,j,k)
!qzz(i,j,k)=two*aa*hh*(one-(two+hh)*el(3)*el(3))/r - f1o3*gzz(i,j,k)*q(i,j,k)
!qxy(i,j,k)=two*aa*hh*(zeo-(two+hh)*el(1)*el(2))/r - f1o3*gxy(i,j,k)*q(i,j,k)
!qxz(i,j,k)=two*aa*hh*(zeo-(two+hh)*el(1)*el(3))/r - f1o3*gxz(i,j,k)*q(i,j,k)
!qyz(i,j,k)=two*aa*hh*(zeo-(two+hh)*el(2)*el(3))/r - f1o3*gyz(i,j,k)*q(i,j,k)
!
!       q(i,j,k) =   q(i,j,k) - trK(i,j,k)
!     qxx(i,j,k) = qxx(i,j,k) - Axx(i,j,k)
!     qxy(i,j,k) = qxy(i,j,k) - Axy(i,j,k)
!     qxz(i,j,k) = qxz(i,j,k) - Axz(i,j,k)
!     qyy(i,j,k) = qyy(i,j,k) - Ayy(i,j,k)
!     qyz(i,j,k) = qyz(i,j,k) - Ayz(i,j,k)
!     qzz(i,j,k) = qzz(i,j,k) - Azz(i,j,k)

    end do
   end do
  end do

  alpha = alpha - one

  return

end subroutine ks_initial_metric_nofisheye

!-----------------------------------------------------------------------------
!
! Non-boosted rotating Kerr-Schild metric initial data + "junk" in the BH interior
!
!-----------------------------------------------------------------------------
  subroutine ks_initial_metric_plus_junk(ex, x, y, z, PhysR,       &
                    gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
                    Axx, Axy, Axz, Ayy, Ayz, Azz, alpha, Bx,By,Bz, &
                    gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                    phi, sam, r_junk, lambda, vx_boost)
  implicit none

!
! Input parameters:
!
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z, PhysR
  real*8,                    intent(in) :: sam, vx_boost
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: trK,alpha,Bx,By,Bz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: phi

!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: el(0:3),boo(0:3),phh(0:3),eta(0:3,0:3),pel(0:3,0:3)
  real*8  :: k1,k2,k3,aa,hh,r2,w2,W,wx,wy,wz,tmp,tmq,tmr
  real*8  :: w0, wh_out, wh_in, psi4, detg, r_junk, lambda
  real*8  :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, lorentz_gamma
  real*8  :: gupxxl,gupxyl,gupxzl,gupyyl,gupyzl,gupzzl
  real*8  :: Kxxl,Kxyl,Kxzl,Kyyl,Kyzl,Kzzl
!  real*8,  dimension(ex(1),ex(2),ex(3)) :: qxx,qxy,qxz,qyy,qyz,qzz,q
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: f1o3=1.d0/3.d0,m0=1.d0,inbd=1.d0,hlf=0.5d0,fou=4.d0
!
! BL radius of the outer and inner horizons 
  wh_out = m0 + sqrt(m0*m0 - sam*sam)
  wh_in  = m0 - sqrt(m0*m0 - sam*sam)
  ! BL radius inside which junk initial data will be filled 
  if (r_junk .lt. 0.d0) then 
     w0 = 0.5d0*(wh_in + wh_out)
     if (abs(sam) .gt. 0.96d0) w0 = 0.6d0
  else
     if (r_junk .gt. wh_out) then 
	write(*,*) 'r_junk = ',r_junk
	write(*,*) 'Horizon radius = ',wh_out
	write(*,*) 'r_junk must be smaller than the horizon radius'
	stop
     end if
     w0 = r_junk
  end if
  
  lorentz_gamma = 1.d0/sqrt(1.d0-vx_boost**2)
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
    !!r2 = PhysR(i,j,k) 
    !!r2 = r2 * r2

    r2 = (lorentz_gamma*X(i,j,k))**2 + Y(i,j,k)**2 + Z(i,j,k)**2

!-- W^2 = 1/2 *[ r^2 - a^2 + \sqrt(( r^2 - a^2 )^2 + 4 a^2 z^2 )] (sqrure of BL radius)

    w2 = r2 - sam * sam
    w2 = hlf *( w2 + sqrt( w2 * w2 + fou * sam * sam * z(i,j,k) * z(i,j,k) )) + m0*1.d-16
     W = sqrt(w2) 

    tmp = two * w2 + sam * sam - r2
    tmq =       w2 + sam * sam

     ! \partial_i W:
     wx =   W * x(i,j,k) / tmp
     wy =   W * y(i,j,k) / tmp
     wz = tmq * z(i,j,k) /(tmp * W)

     ! Modify W when W<w0: W -> w0-m0/lambda+m0/lambda*exp[lambda*(W-w0)]
     if (W .lt. w0) then
	tmr = exp(lambda*(W-w0))
        W = w0 - m0/lambda + m0/lambda*tmr
        w2 = W*W
	wx = wx*tmr
	wy = wy*tmr
	wz = wz*tmr
	tmp = two * w2 + sam * sam - r2
	tmq =       w2 + sam * sam
     end if

!-- H = M W^3 /( w^4 + a^2 z^2 )

    hh = m0 * W*w2 /( w2*w2 + sam * sam * z(i,j,k) * z(i,j,k) )

!
!-- the null vector
!
     el(0) = one
     el(1) = ( W * lorentz_gamma*x(i,j,k) + sam * y(i,j,k) )/ tmq
     el(2) = ( W * y(i,j,k) - sam * lorentz_gamma*x(i,j,k) )/ tmq
     el(3) = z(i,j,k) / W 
     if (z(i,j,k) .eq. 0.d0 .and. r2 .lt. sam*sam) el(3) = sqrt(1.d0 - r2/(sam*sam))

!
!-- \partial_\mu H
!

     tmp = w2 * w2 + sam * sam * z(i,j,k) * z(i,j,k)
     tmr = thr * sam * sam * z(i,j,k) * z(i,j,k) - w2 * w2

     phh(0) = zeo
     phh(1) = hh * tmr * wx /( W * tmp )
     phh(2) = hh * tmr * wy /( W * tmp )
     phh(3) = hh *(tmr * wz /  W - two * sam * sam * z(i,j,k) )/ tmp

!
!-- \partial_\mu L_\nu
!
       pel = zeo

       pel(1,1) = (     W + ( lorentz_gamma*x(i,j,k) - two * W * el(1) )* wx )/ tmq
       pel(2,1) = (   sam + ( lorentz_gamma*x(i,j,k) - two * W * el(1) )* wy )/ tmq
       pel(3,1) =           ( lorentz_gamma*x(i,j,k) - two * W * el(1) )* wz  / tmq

       pel(1,2) = ( - sam + ( y(i,j,k) - two * W * el(2) )* wx )/ tmq
       pel(2,2) = (     W + ( y(i,j,k) - two * W * el(2) )* wy )/ tmq
       pel(3,2) =           ( y(i,j,k) - two * W * el(2) )* wz  / tmq

       pel(1,3) =      - el(3) * wx / W
       pel(2,3) =      - el(3) * wy / W
       pel(3,3) = (one - el(3) * wz)/ W

!
!-- the lapse and shift
!--  \alpha = 1/sqrt{1+2HL_t^2}
!-- \beta^i = 2 H \alpha^2 L_t L^i
!
     alpha(i,j,k) = one/sqrt( one + two*hh*el(0)*el(0) )
               aa = alpha(i,j,k)

        Bx(i,j,k) = two*aa*aa*hh*el(0)*el(1)
        By(i,j,k) = two*aa*aa*hh*el(0)*el(2)
        Bz(i,j,k) = two*aa*aa*hh*el(0)*el(3)

!-- the 3-metric
!

     gxxl = two*hh*el(1)*el(1) + one
     gyyl = two*hh*el(2)*el(2) + one
     gzzl = two*hh*el(3)*el(3) + one
     gxyl = two*hh*el(1)*el(2)
     gyzl = two*hh*el(2)*el(3)
     gxzl = two*hh*el(3)*el(1)

! Convert to the tilde metric
     detg = gxxl * gyyl * gzzl + & 
            gxyl * gyzl * gxzl + & 
            gxzl * gxyl * gyzl &
	    - gxzl * gyyl * gxzl & 
            - gxyl * gxyl * gzzl & 
            - gxxl * gyzl * gyzl

     if (detg .le. 0.d0) then 
        write(*,*) 'Non-positive determinant in Kerr-Schild metric + junk! '
        write(*,*) 'x, y, z, w, detg: ',x(i,j,k),y(i,j,k),z(i,j,k),W,detg
	stop
     end if
     phi(i,j,k) = log(detg)/12.d0
     psi4 = detg**(1.d0/3.d0)
     gxx(i,j,k) = gxxl/psi4
     gxy(i,j,k) = gxyl/psi4
     gxz(i,j,k) = gxzl/psi4
     gyy(i,j,k) = gyyl/psi4
     gyz(i,j,k) = gyzl/psi4
     gzz(i,j,k) = gzzl/psi4

     gupxxl =   ( gyyl * gzzl - gyzl * gyzl )/detg
     gupxyl = - ( gxyl * gzzl - gyzl * gxzl )/detg
     gupxzl =   ( gxyl * gyzl - gyyl * gxzl )/detg
     gupyyl =   ( gxxl * gzzl - gxzl * gxzl )/detg
     gupyzl = - ( gxxl * gyzl - gxyl * gxzl )/detg
     gupzzl =   ( gxxl * gyyl - gxyl * gxyl )/detg

!
!-- the scalar extrinsic curvature
!-- K = \alpha L_t ( 2H @_iL^i -2H @_tL_t -L_t @_tH)
!--    + 2\alpha^3 (1 + HL_t^2)L_tL^i @_i H
!--    + [4\alpha^3(1+HL_t^2) - 2\alpha]HL^i @_iL_t
!
     !!k1 = aa*el(0)*( two*hh*( pel(1,1) + pel(2,2) + pel(3,3) - pel(0,0) ) - &
     !!     el(0)*phh(0) )
     !!k2 = two*aa*aa*aa*el(0)*( el(1)*phh(1) + el(2)*phh(2) +el(3)*phh(3) )* &
     !!     ( one + hh*el(0)*el(0) )
     !!k3 = two*aa*hh*( el(1)*pel(1,0) + el(2)*pel(2,0) + el(3)*pel(3,0) )* &
     !!     ( two*aa*aa*( one + hh*el(0)*el(0) ) - one )

     !!trK(i,j,k) = k1 + k2 + k3

!
!-- A_{ij} = - ( L_iL_j@_tH + HL_i@_tL_j + HL_j@_tL_i )/\alpha
!--          + \alpha( L_tL_i@_jH + L_tL_j@_iH + HL_i@_jL_t + HL_j@_iL_t
!--                  + HL_t@_iL_j + HL_t@_jL_i)
!--          + 2\alpha HL_tL^k( L_iL_j@_kH + HL_i@_kL_j +HL_j@_kL_i )
!--          - \gamma_{ij}K/3
!
     k1 = - ( el(1)*el(1)*phh(0) + hh*el(1)*pel(0,1) + hh*el(1)*pel(0,1) )/aa
     k2 = aa*( el(0)*el(1)*phh(1) + el(0)*el(1)*phh(1) + hh*el(1)*pel(1,0) + &
               hh*el(1)*pel(1,0) + hh*el(0)*( pel(1,1) +  pel(1,1) ) )
     k3 = two*aa*hh*el(0)*( &
          el(1)*el(1)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(1)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) + &
          hh*el(1)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) )

     Kxxl = k1 + k2 + k3
     !!Axx(i,j,k) = (k1 + k2 + k3)/psi4 - f1o3*gxx(i,j,k)*trK(i,j,k)

     k1 = - ( el(1)*el(2)*phh(0) + hh*el(1)*pel(0,2) + hh*el(2)*pel(0,1) )/aa
     k2 = aa*( el(0)*el(1)*phh(2) + el(0)*el(2)*phh(1) + hh*el(1)*pel(2,0) + &
               hh*el(2)*pel(1,0) + hh*el(0)*( pel(1,2) +  pel(2,1) ) )
     k3 = two*aa*hh*el(0)*( &
          el(1)*el(2)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(1)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) + &
          hh*el(2)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) )

     Kxyl = k1 + k2 + k3
     !!Axy(i,j,k) = (k1 + k2 + k3)/psi4 - f1o3*gxy(i,j,k)*trK(i,j,k)

     k1 = - ( el(1)*el(3)*phh(0) + hh*el(1)*pel(0,3) + hh*el(3)*pel(0,1) )/aa
     k2 = aa*( el(0)*el(1)*phh(3) + el(0)*el(3)*phh(1) + hh*el(1)*pel(3,0) + &
               hh*el(3)*pel(1,0) + hh*el(0)*( pel(1,3) +  pel(3,1) ) )
     k3 = two*aa*hh*el(0)*( &
          el(1)*el(3)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(1)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) + &
          hh*el(3)*( el(1)*pel(1,1) + el(2)*pel(2,1) + el(3)*pel(3,1) ) )

     Kxzl = k1 + k2 + k3
     !!Axz(i,j,k) = (k1 + k2 + k3)/psi4 - f1o3*gxz(i,j,k)*trK(i,j,k)

     k1 = - ( el(2)*el(2)*phh(0) + hh*el(2)*pel(0,2) + hh*el(2)*pel(0,2) )/aa
     k2 = aa*( el(0)*el(2)*phh(2) + el(0)*el(2)*phh(2) + hh*el(2)*pel(2,0) + &
               hh*el(2)*pel(2,0) + hh*el(0)*( pel(2,2) +  pel(2,2) ) )
     k3 = two*aa*hh*el(0)*( &
          el(2)*el(2)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(2)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) + &
          hh*el(2)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) )

     Kyyl = k1 + k2 + k3
     !!Ayy(i,j,k) = (k1 + k2 + k3)/psi4 - f1o3*gyy(i,j,k)*trK(i,j,k)

     k1 = - ( el(2)*el(3)*phh(0) + hh*el(2)*pel(0,3) + hh*el(3)*pel(0,2) )/aa
     k2 = aa*( el(0)*el(2)*phh(3) + el(0)*el(3)*phh(2) + hh*el(2)*pel(3,0) + &
               hh*el(3)*pel(2,0) + hh*el(0)*( pel(2,3) +  pel(3,2) ) )
     k3 = two*aa*hh*el(0)*( &
          el(2)*el(3)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(2)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) + &
          hh*el(3)*( el(1)*pel(1,2) + el(2)*pel(2,2) + el(3)*pel(3,2) ) )

     Kyzl = k1 + k2 + k3
     !!Ayz(i,j,k) = (k1 + k2 + k3)/psi4 - f1o3*gyz(i,j,k)*trK(i,j,k)

     k1 = - ( el(3)*el(3)*phh(0) + hh*el(3)*pel(0,3) + hh*el(3)*pel(0,3) )/aa
     k2 = aa*( el(0)*el(3)*phh(3) + el(0)*el(3)*phh(3) + hh*el(3)*pel(3,0) + &
               hh*el(3)*pel(3,0) + hh*el(0)*( pel(3,3) +  pel(3,3) ) )
     k3 = two*aa*hh*el(0)*( &
          el(3)*el(3)*( el(1)*phh(1) + el(2)*phh(2) + el(3)*phh(3) ) + &
          hh*el(3)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) + &
          hh*el(3)*( el(1)*pel(1,3) + el(2)*pel(2,3) + el(3)*pel(3,3) ) )

     Kzzl = k1 + k2 + k3
     !!Azz(i,j,k) = (k1 + k2 + k3)/psi4 - f1o3*gzz(i,j,k)*trK(i,j,k)

     trK(i,j,k) = gupxxl*Kxxl + 2.d0*gupxyl*Kxyl + 2.d0*gupxzl*Kxzl + & 
     		  gupyyl*Kyyl + 2.d0*gupyzl*Kyzl + gupzzl*Kzzl

     Axx(i,j,k) = (Kxxl - f1o3*gxxl*trK(i,j,k))/psi4
     Axy(i,j,k) = (Kxyl - f1o3*gxyl*trK(i,j,k))/psi4
     Axz(i,j,k) = (Kxzl - f1o3*gxzl*trK(i,j,k))/psi4
     Ayy(i,j,k) = (Kyyl - f1o3*gyyl*trK(i,j,k))/psi4
     Ayz(i,j,k) = (Kyzl - f1o3*gyzl*trK(i,j,k))/psi4
     Azz(i,j,k) = (Kzzl - f1o3*gzzl*trK(i,j,k))/psi4

    end do
   end do
  end do

  alpha = alpha - one

  gupxx =   ( gyy * gzz - gyz * gyz )
  gupxy = - ( gxy * gzz - gyz * gxz )
  gupxz =   ( gxy * gyz - gyy * gxz )
  gupyy =   ( gxx * gzz - gxz * gxz )
  gupyz = - ( gxx * gyz - gxy * gxz )
  gupzz =   ( gxx * gyy - gxy * gxy )

  return

end subroutine ks_initial_metric_plus_junk

!-----------------------------------------------------------------------------
!
! Find non-boosted rotating Kerr-Schild metric initial data (with junk)
! This is for the case when we're using the Boyer-Lindquist radius
!
!-----------------------------------------------------------------------------
  subroutine ks_initial_metric_bl_junk(ex, x, y, z, Rp,             &
       gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
       Axx, Axy, Axz, Ayy, Ayz, Azz, alpha, Bx,By,Bz, &
       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
       phi, sam, r_junk, lambda, vx_boost)
  implicit none

!
! Input parameters:
!
  integer, dimension(3), intent(in)                  :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: x,y,z
  real*8,  intent(in)                                :: sam, vx_boost
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: trK,alpha,Bx,By,Bz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: phi
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: Rp
!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: k1,k2,k3,aa,hh,r2,w2,W,wx,wy,wz,tmp,tmq,tmr, lorentz_gamma
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: f1o3=1.d0/3.d0,m0=1.d0,inbd=1.d0,hlf=0.5d0,fou=4.d0

  real*8 :: grr, grph, gthth, gphph, Br, Bph
  real*8 :: guprr, guprph, gupthth, gupphph,psi4
  real*8 :: Arr, Arth, Arph, Athth, Athph, Aphph
  real*8 :: Sigma, r, sin2th, sinth, cos2th, costh
  real*8 :: Lx_r, Ly_r, Lz_r, Lx_th, Ly_th, Lz_th, Lx_ph, Ly_ph, Lz_ph
  real*8 :: Lr_x, Lth_x, Lph_x, Lr_y, Lth_y, Lph_y, Lr_z, Lth_z, Lph_z
  real*8 :: r_l, x_l, y_l, z_l
  real*8 :: grrtemp, grphtemp, gththtemp, gphphtemp, gupxxtemp
  real*8 :: gupxxl, gupxyl, gupxzl, gupyyl, gupyzl, gupzzl, detg_l
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

  lorentz_gamma = 1.d0/sqrt(1.d0-vx_boost**2)

! Input translation
!
  imin = lbound(gxx,1)
  jmin = lbound(gxx,2)
  kmin = lbound(gxx,3)
  imax = ubound(gxx,1)
  jmax = ubound(gxx,2)
  kmax = ubound(gxx,3)

!
! go to each gridpoint...
!
  do k = kmin, kmax
   do j = jmin, jmax
    do i = imin, imax

    !r = Rp(i,j,k)
    r = sqrt( (lorentz_gamma*X(i,j,k))**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
    z_l = z(i,j,k)
    costh = z_l/r
    cos2th = costh**2
    sin2th = one - cos2th
    sinth = sqrt(sin2th)
    !if (r .lt. r0) then
    !   ! Modify r when r<r0: r -> r0-m0/lambda+m0/lambda*exp[lambda*(r-r0)]
    !   r = r0 - m0/lambda + m0/lambda*exp(lambda*(r-r0))
    !end if
    r2 = r * r
    Sigma = r2 + sam*sam*cos2th
!
!-- the lapse and shift (nonzero component)
!
     alpha(i,j,k) = 1.d0/sqrt(one+two*r/Sigma)
      Br = two*r/(Sigma + two*r)
!
!-- the 3-metric (nonzero components)
!
     grr0_pl = 1.d0+2.d0/r0
     dgrr0_pl = -2.d0/r0/r0
     aaa = 0.5d0*dgrr0_pl/r0 
     psi40 = grr0_pl - aaa*r0*r0
     alp0_pl = 1.d0/sqrt(1.d0+2.d0/r0)
     dalp0_pl = 1.d0/r0/r0 * alp0_pl**3
     aaa = 0.5d0*dalp0_pl/r0
     alp0_pl = alp0_pl - aaa*r0*r0
     if (r .ge. r0) then 
        grr = one+two*r/Sigma
        grph = -sam*sin2th*grr 
        gthth = Sigma
        gphph = sin2th*(Sigma + sam**2*sin2th*grr)
     else
	Sigma0 = r0*r0+sam*sam*cos2th
	grr0 = one+two*r0/Sigma0
        grph_r = -sam*sin2th*grr0/r0
        gthth_r2 = Sigma0/r0/r0
        gphph_r2 = sin2th*(Sigma0 + sam**2*sin2th*grr0)/r0/r0
	dgrr0 = -2.d0*(r0*r0-sam*sam*cos2th)/(r0*r0+sam*sam*cos2th)**2
	dgrph_r = -sam*sin2th*dgrr0/r0 + sam*sin2th*grr0/r0/r0
	dgthth_r2 = -2.d0*sam*sam*cos2th/r0**3
	dgphph_r2 = sin2th*(2.d0*r0 + sam*sam*sin2th*dgrr0)/r0/r0 - 2.d0*gphph_r2/r0
        call compute_ab(r0,grr0,dgrr0,psi40,aaa,bbb)
	grr = psi40 + aaa*r2 + bbb*r2*r
	call compute_ab(r0,grph_r,dgrph_r,0.d0,aaa,bbb)
	grph = r*(aaa*r2 + bbb*r2*r)
        call compute_ab(r0,gthth_r2,dgthth_r2,psi40,aaa,bbb)
	gthth = r2*(psi40 + aaa*r2 + bbb*r2*r)
	call compute_ab(r0,gphph_r2,dgphph_r2,psi40*sin2th,aaa,bbb)
	gphph = r2*(psi40*sin2th + aaa*r2 + bbb*r2*r)
	alp0 = 1.d0/sqrt(one+two*r0/Sigma0)
	dalp0 = (r0*r0-sam*sam*cos2th)/Sigma0/Sigma0*alp0**3
        call compute_ab(r0,alp0,dalp0,alp0_pl,aaa,bbb)
	alpha(i,j,k) = alp0_pl + aaa*r2 + bbb*r2*r
	if (alpha(i,j,k) .lt. 0.d0) then 
	   write(*,*) 'X, Y, Z = ',X(i,j,k),Y(i,j,k),Z(i,j,k)
	   write(*,*) 'lapse = ',alpha(i,j,k), ' < 0!'
	   stop
	end if
       
        !!grr = (1.d0 + 0.5d0/r)**4
        !!grph = 0.d0
        !!gthth = grr*r2
	!!gphph = sin2th*gthth
        
     end if

     if (r .ge. r0) then 
        guprr = (Sigma + sam**2*sin2th*grr)/(grr*Sigma)
        guprph = sam/Sigma
        gupthth = 1.d0/Sigma
        gupphph = 1.d0/(sin2th*Sigma)
     else
	detg_l = grr*gthth*gphph - grph*grph*gthth
	guprr = gthth*gphph/detg_l 
	guprph = -gthth*grph/detg_l
	gupthth = (grr*gphph - grph*grph)/detg_l
	gupphph = grr*gthth/detg_l
     end if


!
!-- Store K_{ij} in A_{ij} 
!--         
!
     if (r .ge. r0) then 
        Arr = two*(two*Sigma - 4.d0*r2)/(two*alpha(i,j,k))
        Arr = Arr/(Sigma**2*(Sigma + two*r)**2)
        Arr = Arr*(3.d0*r*Sigma + two*r2 + Sigma**2)

        Arth = 4.d0*r*sam**2*sinth*costh
        Arth = Arth / (Sigma*(Sigma+two*r))/(two*alpha(i,j,k))
        
        Arph = -two*sam*sin2th*(one - two*r2/Sigma)/Sigma/(two*alpha(i,j,k))

        Athth = 4.d0*r2/(Sigma + two*r)/(two*alpha(i,j,k))

        Athph = -4.d0*r*sam**3*sin2th*sinth*costh
        Athph = Athph/(Sigma*(Sigma+two*r))/(two*alpha(i,j,k))

        Aphph = 4.d0*r*sin2th/(Sigma + two*r)
        Aphph = Aphph*(r + sam*sam*sin2th*((Sigma-two*r2)/Sigma**2))/(two*alpha(i,j,k))
     else
	Sigma0 = r0*r0 + sam*sam*cos2th
	f1oalpha0 = sqrt(1.d0+2.d0*r0/Sigma0)

	Arr0 = (two*Sigma0 - 4.d0*r0*r0)*f1oalpha0
        Arr0 = Arr0/(Sigma0**2*(Sigma0 + two*r0)**2)
        Arr0 = Arr0*(3.d0*r0*Sigma0 + two*r0*r0 + Sigma0**2)
	dArr0 = -(r0*r0 - sam*sam*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*r0/Sigma0) & 
		- 2.d0*r0/(Sigma0 - 2.d0*r0*r0) + & 
		(Sigma0*(4.d0*r0+3.d0)+4.d0*r0+6.d0*r0*r0)/(Sigma0**2 + 3.d0*r0*Sigma0 + 2.d0*r0*r0) & 
		- 4.d0*r0/Sigma0 - 4.d0*(r0+1.d0)/(Sigma0+2.d0*r0)
	dArr0 = Arr0*dArr0
        call compute_ab(r0,Arr0,dArr0,0.d0,aaa,bbb)
	Arr = aaa*r2 + bbb*r2*r

        Arth_r = 2.d0*sam*sam*sinth*costh/(Sigma0*(Sigma0+two*r0))*f1oalpha0
	dArth_r = -(r0*r0 - sam*sam*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*r0/Sigma0) & 
		  - 2.d0*r0/Sigma0 - 2.d0*(r0+1.d0)/(Sigma0+2.d0*r0)
	dArth_r = Arth_r*dArth_r
	call compute_ab(r0,Arth_r,dArth_r,0.d0,aaa,bbb)
	Arth = r*(aaa*r2 + bbb*r2*r)

	Arph_r = -sam*sin2th*(one - two*r0*r0/Sigma0)/r0/Sigma0*f1oalpha0
	dArph_r = -(r0*r0 - sam*sam*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*r0/Sigma0) &
		  - 4.d0*sam*sam*cos2th*r0/Sigma0/Sigma0/(1.d0-2.d0*r0*r0/Sigma0) & 
		  - 1.d0/r0 - 2.d0*r0/Sigma0
	dArph_r = dArph_r * Arph_r
	call compute_ab(r0,Arph_r,dArph_r,0.d0,aaa,bbb)
	Arph = r*(aaa*r2 + bbb*r2*r)

        Athth_r2 = 2.d0/(Sigma0 + two*r0)*f1oalpha0
	dAthth_r2 = -(r0*r0 - sam*sam*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*r0/Sigma0) &
		  	- 2.d0*(r0+1.d0)/(Sigma0+2.d0*r0)
	dAthth_r2 = dAthth_r2 *Athth_r2
	call compute_ab(r0,Athth_r2,dAthth_r2,0.d0,aaa,bbb)
	Athth = r2*(aaa*r2 + bbb*r2*r)

        Athph_r2 = -2.d0*sam**3*sin2th*sinth*costh/r0
        Athph_r2 = Athph_r2/(Sigma0*(Sigma0+two*r0))*f1oalpha0
	dAthph_r2 = -(r0*r0 - sam*sam*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*r0/Sigma0) &
		    - 1.d0/r0 - 2.d0*r0/Sigma0 - 2.d0*(r0+1.d0)/(Sigma0+2.d0*r0)
	dAthph_r2 = dAthph_r2 *Athph_r2
	call compute_ab(r0,Athph_r2,dAthph_r2,0.d0,aaa,bbb)
	Athph = r2*(aaa*r2 + bbb*r2*r)

        Aphph_r2 = 2.d0*sin2th/(Sigma0 + two*r0)/r0
        Aphph_r2 = Aphph_r2*(r0 + sam*sam*sin2th*((Sigma0-two*r0*r0)/Sigma0**2))*f1oalpha0
	dAphph_r2 = -(r0*r0 - sam*sam*cos2th)/Sigma0/Sigma0/(1.d0+2.d0*r0/Sigma0) &
		    - 1.d0/r0 - 2.d0*(r0+1.d0)/(Sigma0+2.d0*r0) + & 
		    (1.d0+4.d0*sam*sam*sin2th*r0*(r0*r0-sam*sam*cos2th)/Sigma0**3 & 
		     - 2.d0*sam*sam*sin2th*r0/Sigma0**2) & 
		    / (r0 + sam*sam*sin2th*((Sigma0-two*r0*r0)/Sigma0**2))
	dAphph_r2 = dAphph_r2 * Aphph_r2
	call compute_ab(r0,Aphph_r2,dAphph_r2,0.d0,aaa,bbb)
	Aphph = r2*(aaa*r2 + bbb*r2*r)

	!!Arr = 0.d0
    	!!Arth = 0.d0
	!!Arph = 0.d0
	!!Athth = 0.d0
	!!Athph = 0.d0
	!!Aphph = 0.d0
     end if


!
!-- the scalar extrinsic curvature
!
     !!trK(i,j,k) = 4.d0/(Sigma*(Sigma+two*r)**3)*( &
     !!     (Sigma+sam*sam*sin2th*(one+two*r/Sigma))* &
     !!     (3.d0*r*Sigma + two*r2 + Sigma*Sigma)* &
     !!     (one - two*r2/Sigma))/(two*alpha(i,j,k))
     !!trK(i,j,k) = trK(i,j,k) - 4.d0*sam*sam*sin2th/Sigma**2* &
     !!     (one - two*r2/Sigma)/(two*alpha(i,j,k)) 
     !!trK(i,j,k) = trK(i,j,k) + 4.d0*r/(Sigma*(Sigma+two*r))* &
     !!     (two*r + sam*sam*sin2th*((Sigma-two*r2)/Sigma**2))/(two*alpha(i,j,k))

     !!trK(i,j,K) = guprr*Arr + two*guprph*Arph + &
     !!             gupthth*Athth + gupphph*Aphph

!
!-- A_{ij} = K_{ij} - \gamma_{ij}K/3 
!--         
!
     !!Arr = Arr - f1o3*grr*trK(i,j,k)
     !!Arph = Arph - f1o3*grph*trK(i,j,k)
     !!Athth = Athth - f1o3*gthth*trK(i,j,k)
     !!Aphph = Aphph - f1o3*gphph*trK(i,j,k)


     !#############################################################################
     !Now convert everything to Cartesian
     !#############################################################################

       r_l = sqrt( (lorentz_gamma*x(i,j,k))**2 + y(i,j,k)**2 + z(i,j,k)**2)
       !r_l = r
       r2 = r_l * r_l

       x_l = x(i,j,k)*lorentz_gamma
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

       detg_l = gxxl * gyyl * gzzl + &
                gxyl * gyzl * gxzl + &
                gxzl * gxyl * gyzl &
                - gxzl * gyyl * gxzl &
                - gxyl * gxyl * gzzl &
                - gxxl * gyzl * gyzl

       phi(i,j,k) = log(detg_l)/12.d0
       psi4 = detg_l**(1.d0/3.d0)
       gxx(i,j,k) = gxxl/psi4
       gxy(i,j,k) = gxyl/psi4
       gxz(i,j,k) = gxzl/psi4
       gyy(i,j,k) = gyyl/psi4
       gyz(i,j,k) = gyzl/psi4
       gzz(i,j,k) = gzzl/psi4

       gupxxl =   ( gyyl * gzzl - gyzl * gyzl )/detg_l
       gupxyl = - ( gxyl * gzzl - gyzl * gxzl )/detg_l
       gupxzl =   ( gxyl * gyzl - gyyl * gxzl )/detg_l
       gupyyl =   ( gxxl * gzzl - gxzl * gxzl )/detg_l
       gupyzl = - ( gxxl * gyzl - gxyl * gxzl )/detg_l
       gupzzl =   ( gxxl * gyyl - gxyl * gxyl )/detg_l

       trK(i,j,k) = gupxxl*Kxxl + 2.d0*gupxyl*Kxyl + 2.d0*gupxzl*Kxzl + &
                    gupyyl*Kyyl + 2.d0*gupyzl*Kyzl + gupzzl*Kzzl

       Axx(i,j,k) = (Kxxl - f1o3*gxxl*trK(i,j,k))/psi4
       Axy(i,j,k) = (Kxyl - f1o3*gxyl*trK(i,j,k))/psi4
       Axz(i,j,k) = (Kxzl - f1o3*gxzl*trK(i,j,k))/psi4
       Ayy(i,j,k) = (Kyyl - f1o3*gyyl*trK(i,j,k))/psi4
       Ayz(i,j,k) = (Kyzl - f1o3*gyzl*trK(i,j,k))/psi4
       Azz(i,j,k) = (Kzzl - f1o3*gzzl*trK(i,j,k))/psi4

    end do
   end do
  end do

  alpha = alpha - one

  gupxx =   ( gyy * gzz - gyz * gyz )
  gupxy = - ( gxy * gzz - gyz * gxz )
  gupxz =   ( gxy * gyz - gyy * gxz )
  gupyy =   ( gxx * gzz - gxz * gxz )
  gupyz = - ( gxx * gyz - gxy * gxz )
  gupzz =   ( gxx * gyy - gxy * gxy )

  !!trK = -trK
  !!Axx = -Axx
  !!Axy = -Axy
  !!Axz = -Axz
  !!Ayy = -Ayy
  !!Ayz = -Ayz
  !!Azz = -Azz

  return

end subroutine ks_initial_metric_bl_junk

!-----------------------------------------------------------------------------
!
! Find non-boosted rotating Kerr-Schild "puncture" initial data 
!
!-----------------------------------------------------------------------------
  subroutine kerr_initial_metric_puncture(ex, x, y, z, Rp,             &
       gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
       Axx, Axy, Axz, Ayy, Ayz, Azz, alpha, Bx,By,Bz, &
       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
       phi, sam, vx_boost)
  implicit none

!
! Input parameters:
!
  integer, dimension(3), intent(in)                  :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: x,y,z
  real*8,  intent(in)                                :: sam, vx_boost
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: trK,alpha,Bx,By,Bz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: phi
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: Rp
!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: k1,k2,k3,aa,hh,r2,w2,W,wx,wy,wz,tmp,tmq,tmr
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: f1o3=1.d0/3.d0,m0=1.d0,inbd=1.d0,hlf=0.5d0,fou=4.d0

  real*8 :: grr, grph, gthth, gphph, alp,Bph
  real*8 :: guprr, guprph, gupthth, gupphph,psi4
  real*8 :: Arr, Arth, Arph, Athth, Athph, Aphph
  real*8 :: Krr, Krth,Krph,Kthth,Kthph,Kphph
  real*8 :: Sigma, r, sin2th, sinth, cos2th, costh, lorentz_gamma
  real*8 :: Lx_r, Ly_r, Lz_r, Lx_th, Ly_th, Lz_th, Lx_ph, Ly_ph, Lz_ph
  real*8 :: Lr_x, Lth_x, Lph_x, Lr_y, Lth_y, Lph_y, Lr_z, Lth_z, Lph_z
  real*8 :: r_l, x_l, y_l, z_l, Delta, dtdr, dpdr
  real*8 :: grrtemp, grphtemp, gththtemp, gphphtemp, gupxxtemp
  real*8 :: gupxxl, gupxyl, gupxzl, gupyyl, gupyzl, gupzzl, detg_l
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
  imin = lbound(gxx,1)
  jmin = lbound(gxx,2)
  kmin = lbound(gxx,3)
  imax = ubound(gxx,1)
  jmax = ubound(gxx,2)
  kmax = ubound(gxx,3)

  lorentz_gamma = 1.d0/sqrt(1.d0-vx_boost*vx_boost)

!
! go to each gridpoint...
!
  do k = kmin, kmax
   do j = jmin, jmax
    do i = imin, imax

    riso = sqrt( (lorentz_gamma*X(i,j,k))**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
    z_l = Z(i,j,k)
    r = riso*(1.d0 + 0.25d0*rh_out/riso)**2
    dr_driso = 1.d0 - (0.25d0*rh_out/riso)**2
    costh = z_l/riso
    cos2th = costh**2
    sin2th = one - cos2th
    sinth = sqrt(sin2th)
    r2 = r * r
    Sigma = r2 + sam*sam*cos2th
    Delta = r2 - 2.d0*r + a*a
    aa = (r2+a*a)**2 - a*a*Delta*sin2th
!
!-- the lapse and shift (nonzero component)
!
     alp = sqrt(Delta*Sigma/aa)
     alpha(i,j,k) = alp
      Bph = -2.d0*a*r/aa

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
     Krph = (a*sin2th*(-a**4 + 2.d0*a*a*r2 + 3.d0*r2*r2 + &
            a*a*(a - r)*(a + r)*sin2th))/Sigma/aa * dr_driso_alp
     Kthth = 0.d0
     Kthph = -2.d0*a*a*a*r*costh*sinth*sin2th/Sigma * sqrt(Delta/Sigma/aa)
     if (riso .lt. 0.25d0*rh_out) Kthph = -Kthph
     Kphph = 0.d0
!
!-- the scalar extrinsic curvature
!
     trK(i,j,k) = 0.d0
!
!-- A_{ij} = K_{ij} - \gamma_{ij}K/3 
!--         
!
     Arr = Krr
     Arth = Krth
     Arph = Krph
     Athth = Kthth
     Athph = Kthph
     Aphph = Kphph

     !#############################################################################
     !Now convert everything to Cartesian
     !#############################################################################
       r_l = sqrt( (lorentz_gamma*x(i,j,k))**2 + y(i,j,k)**2 + z(i,j,k)**2)
       r2 = r_l * r_l

       x_l = x(i,j,k)*lorentz_gamma
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

       Bx(i,j,k) = Lx_ph*Bph
       By(i,j,k) = Ly_ph*Bph
       Bz(i,j,k) = 0.d0

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
       Axx(i,j,k) = Lr_x**2 * Arr + &
		    2.d0*Lr_x*Lth_x * Arth + &
                    2.d0*Lr_x*Lph_x * Arph + &
                    Lth_x**2 * Athth + &
		    2.d0*Lth_x*Lph_x * Athph + &
                    Lph_x**2 * Aphph

       Axy(i,j,k) = Lr_x*Lr_y * Arr + &
		    (Lr_x*Lth_y + Lr_y*Lth_x) * Arth + &
                    (Lr_x*Lph_y + Lr_y*Lph_x) * Arph + &
                    Lth_x*Lth_y * Athth + &
		    (Lth_x*Lph_y + Lth_y*Lph_x) * Athph + &
                    Lph_x*Lph_y * Aphph

       Axz(i,j,k) = Lr_x*Lr_z * Arr + &
		    (Lr_x*Lth_z + Lr_z*Lth_x) * Arth + &
                    (Lr_x*Lph_z + Lr_z*Lph_x) * Arph + &
                    Lth_x*Lth_z * Athth + &
		    (Lth_x*Lph_z + Lth_z*Lph_x) * Athph + &
                    Lph_x*Lph_z * Aphph

       Ayy(i,j,k) = Lr_y**2 * Arr + &
		    2.d0*Lr_y*Lth_y * Arth + &
                    2.d0*Lr_y*Lph_y * Arph + &
                    Lth_y**2 * Athth + &
		    2.d0*Lth_y*Lph_y * Athph + &
                    Lph_y**2 * Aphph

       Ayz(i,j,k) = Lr_y*Lr_z * Arr + &
		    (Lr_y*Lth_z + Lr_z*Lth_y) * Arth + &
                    (Lr_y*Lph_z + Lr_z*Lph_y) * Arph + &
                    Lth_y*Lth_z * Athth + &
		    (Lth_y*Lph_z + Lth_z*Lph_y) * Athph + &
                    Lph_y*Lph_z * Aphph

       Azz(i,j,k) = Lr_z**2 * Arr + &
		    2.d0*Lr_z*Lth_z * Arth + &
                    2.d0*Lr_z*Lph_z * Arph + &
                    Lth_z**2 * Athth + &
		    2.d0*Lth_z*Lph_z * Athph + &
                    Lph_z**2 * Aphph

       detg_l = gxxl * gyyl * gzzl + &
                gxyl * gyzl * gxzl + &
                gxzl * gxyl * gyzl &
                - gxzl * gyyl * gxzl &
                - gxyl * gxyl * gzzl &
                - gxxl * gyzl * gyzl

       phi(i,j,k) = log(detg_l)/12.d0
       psi4 = detg_l**(1.d0/3.d0)
       gxx(i,j,k) = gxxl/psi4
       gxy(i,j,k) = gxyl/psi4
       gxz(i,j,k) = gxzl/psi4
       gyy(i,j,k) = gyyl/psi4
       gyz(i,j,k) = gyzl/psi4
       gzz(i,j,k) = gzzl/psi4
       Axx(i,j,k) = Axx(i,j,k)/psi4
       Axy(i,j,k) = Axy(i,j,k)/psi4
       Axz(i,j,k) = Axz(i,j,k)/psi4
       Ayy(i,j,k) = Ayy(i,j,k)/psi4
       Ayz(i,j,k) = Ayz(i,j,k)/psi4
       Azz(i,j,k) = Azz(i,j,k)/psi4

    end do
   end do
  end do

  alpha = alpha - one

  gupxx =   ( gyy * gzz - gyz * gyz )
  gupxy = - ( gxy * gzz - gyz * gxz )
  gupxz =   ( gxy * gyz - gyy * gxz )
  gupyy =   ( gxx * gzz - gxz * gxz )
  gupyz = - ( gxx * gyz - gxy * gxz )
  gupzz =   ( gxx * gyy - gxy * gxy )

  return

end subroutine kerr_initial_metric_puncture


!-----------------------------------------------------------------------------
!
! Given a function f(r), find the values of a and b so that 
!   g(r) = f0 + a r^2 + b r^3, with g(r0) = f(r0) and g'(r0)=f'(r0). 
! Here f0 is a constant.
!
!-----------------------------------------------------------------------------
subroutine compute_ab(r0,fr0,dfr0,f0,a,b)
   implicit none
   real*8 :: f0,fr0,dfr0,a,b,r0 
!
   b = dfr0/r0**2 + 2.d0*(f0-fr0)/r0**3
   a = (fr0-f0)/r0**2 - b*r0
end subroutine compute_ab


!-----------------------------------------------------------------------------
!
! Convert untilded metric to tilded metric with analytical value 
!
!-----------------------------------------------------------------------------
  subroutine convert_a(ex,    x,     y,     z, PhysR,    phi,          &
                       gxx,   gxy,   gxz,   gyy,   gyz,   gzz,   & 
                       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                       gxxx,  gxxy,  gxxz,  gxyx,  gxyy,  gxyz,  &
                       gxzx,  gxzy,  gxzz,  gyyx,  gyyy,  gyyz,  &
                       gyzx,  gyzy,  gyzz,  gzzx,  gzzy,  gzzz,  &
                       gamx,  gamy,  gamz,  phix,  phiy,  phiz,  &
                       lapx,  lapy,  lapz)
  implicit none

!~~~~~> input variables

  integer, intent(in) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z,PhysR
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout):: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: phi
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: gxxx,gxxy,gxxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: gxyx,gxyy,gxyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: gxzx,gxzy,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: gyyx,gyyy,gyyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: gyzx,gyzy,gyzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: gzzx,gzzy,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: gamx,gamy,gamz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: phix,phiy,phiz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)  :: lapx,lapy,lapz

!~~~~~> local variable

  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,m,n
  real*8  :: a2, hh, rr, dtm, idtm1o3, dtm1o3, a1
  real*8  :: el(0:3),boo(0:3),phh(0:3),eta(0:3,0:3),pel(0:3,0:3)
  real*8  :: LA,LH,LX,LY,LZ,PL,WX,WY,WZ
  real*8, parameter :: f1o3=1.d0/3.d0, f1o12=1.d0/12.d0, m0=1.d0, thr=3.d0
  real*8, parameter :: zeo= 0.d0, one=1.d0, two=2.d0, fou=4.d0, f4o3=4.d0/3.d0
  real*8, parameter :: f1o6=1.d0/6.d0
 
!~~~~~> Input translation
 
  imin = lbound(x,1)
  jmin = lbound(y,1)
  kmin = lbound(z,1)
  imax = ubound(x,1)
  jmax = ubound(y,1)
  kmax = ubound(z,1)
 
!
!-- flat metric
!

  eta      =  zeo
  eta(0,0) = -one
  eta(1,1) =  one
  eta(2,2) =  one
  eta(3,3) =  one

  boo(0)   = -one
  boo(1)   =  zeo
  boo(2)   =  zeo
  boo(3)   =  zeo

!~~~~~~> computation starts
 
  do k=kmin,kmax,1
   do j=jmin,jmax,1
    do i=imin,imax,1

!~~~~~~> H = M / r, 

     rr = PhysR(i,j,k)
     hh = m0 / rr

!~~~~~~> the null vector
!   \ell_t = \gamma -\gamma^2 (v_x x + v_y y) / r
!   \ell_x = x / r - v_x \ell_t
!   \ell_y = y / r - v_y \ell_t
!   \ell_z = z / r

     el(0) = one
     el(1) = x(i,j,k) / rr
     el(2) = y(i,j,k) / rr
     el(3) = z(i,j,k) / rr

!~~~~~~> @_\mu H = - (H/r) [L_\mu + B_\mu]
 
     phh = - hh *( el + boo )/ rr

!--- @_\mu L_\nu = (eta_{\mu\nu} - L_\mu B_\nu -L_\nu B_\mu - L_\mu L_\nu )/r
!
     do m = 0,3,1
      do n = 0,3,1
       pel(m,n) = ( eta(m,n) - el(m)*boo(n) - el(n)*boo(m) - el(m)*el(n) )/ rr
      end do
     end do


!~~~~~~>    detm = 1+2H(\ell_t)^2 = e^(12\phi)
!        alpha^2 = 1/[1+2H(\ell_t)^2], 

     dtm = one + two * hh * el(0) * el(0)
      a2 = one / dtm
      a1 = sqrt(a2)

     dtm1o3 = dtm**f1o3
     idtm1o3 = a2**f1o3

!~~~~~~> derivative of phi
!       \phi_{,i} = \alpha^2 \ell_t ( \ell_t H_{,i} + 2H @_i \ell_t ) / 6
 
     phix(i,j,k) = f1o6 * a2 * el(0) *( el(0) * phh(1) + two * hh * pel(1,0) )
     phiy(i,j,k) = f1o6 * a2 * el(0) *( el(0) * phh(2) + two * hh * pel(2,0) )
     phiz(i,j,k) = f1o6 * a2 * el(0) *( el(0) * phh(3) + two * hh * pel(3,0) )
    
!~~~~~~> derivative of lapse
!       \alpha_{,i} = - \alpha^3 \ell_t ( \ell_t H_{,i} + 2H @_i \ell_t )

     lapx(i,j,k) = - a2 * a1 * el(0) *( el(0) * phh(1) + two * hh * pel(1,0) )
     lapy(i,j,k) = - a2 * a1 * el(0) *( el(0) * phh(2) + two * hh * pel(2,0) )
     lapz(i,j,k) = - a2 * a1 * el(0) *( el(0) * phh(3) + two * hh * pel(3,0) )
 
!~~~~~~> derivative of shift
!        B^i_{,j} = 2\alpha ( 2 H \ell_t \ell^i \alpha_{,j}
!                             + \alpha \ell_t \ell^i H_{,j}
!                             + \alpha H \ell^i @_j \ell_t
!                             + \alpha H \ell_t @_j \ell^i )

!sfxx(i,j,k)= two*a1*( two*hh*el(0)*el(1)*lapx(i,j,k) + a1*el(0)*el(1)*phh(1) +&
!                                a1*hh*el(1)*pel(1,0) + a1*hh*el(0)*pel(1,1)   )
!sfxy(i,j,k)= two*a1*( two*hh*el(0)*el(1)*lapy(i,j,k) + a1*el(0)*el(1)*phh(2) +&
!                                a1*hh*el(1)*pel(2,0) + a1*hh*el(0)*pel(2,1)   )
!sfxz(i,j,k)= two*a1*( two*hh*el(0)*el(1)*lapz(i,j,k) + a1*el(0)*el(1)*phh(3) +&
!                                a1*hh*el(1)*pel(3,0) + a1*hh*el(0)*pel(3,1)   )

!sfyx(i,j,k)= two*a1*( two*hh*el(0)*el(2)*lapx(i,j,k) + a1*el(0)*el(2)*phh(1) +&
!                                a1*hh*el(2)*pel(1,0) + a1*hh*el(0)*pel(1,2)   )
!sfyy(i,j,k)= two*a1*( two*hh*el(0)*el(2)*lapy(i,j,k) + a1*el(0)*el(2)*phh(2) +&
!                                a1*hh*el(2)*pel(2,0) + a1*hh*el(0)*pel(2,2)   )
!sfyz(i,j,k)= two*a1*( two*hh*el(0)*el(2)*lapz(i,j,k) + a1*el(0)*el(2)*phh(3) +&
!                                a1*hh*el(2)*pel(3,0) + a1*hh*el(0)*pel(3,2)   )

!sfzx(i,j,k)= two*a1*( two*hh*el(0)*el(3)*lapx(i,j,k) + a1*el(0)*el(3)*phh(1) +&
!                                a1*hh*el(3)*pel(1,0) + a1*hh*el(0)*pel(1,3)   )
!sfzy(i,j,k)= two*a1*( two*hh*el(0)*el(3)*lapy(i,j,k) + a1*el(0)*el(3)*phh(2) +&
!                                a1*hh*el(3)*pel(2,0) + a1*hh*el(0)*pel(2,3)   )
!sfzz(i,j,k)= two*a1*( two*hh*el(0)*el(3)*lapz(i,j,k) + a1*el(0)*el(3)*phh(3) +&
!                                a1*hh*el(3)*pel(3,0) + a1*hh*el(0)*pel(3,3)   )

!~~~~~~> divide metric by determinant
!        \tilde{g_ij} = e^{-4\phi}g_{ij}

     gxx(i,j,k) = idtm1o3 * gxx(i,j,k)
     gxy(i,j,k) = idtm1o3 * gxy(i,j,k)
     gxz(i,j,k) = idtm1o3 * gxz(i,j,k)
     gyy(i,j,k) = idtm1o3 * gyy(i,j,k)
     gyz(i,j,k) = idtm1o3 * gyz(i,j,k)
     gzz(i,j,k) = idtm1o3 * gzz(i,j,k)

!~~~~~~> determine phi
!        \phi = ln(detm)/12 = ln(1+2H)/12

     phi(i,j,k) = f1o12 * log(dtm)

!~~~~~~> invert metric
!        g^{ij} = \eta^{ij} - 2\alpha^2 H \ell^i \ell^j
!        \tilde{g^{ij}} = e^{4\phi}g^{ij}

     gupxx(i,j,k) = dtm1o3 *( eta(1,1) - two * a2 * hh * el(1) * el(1) )
     gupxy(i,j,k) = dtm1o3 *( eta(1,2) - two * a2 * hh * el(1) * el(2) )
     gupxz(i,j,k) = dtm1o3 *( eta(1,3) - two * a2 * hh * el(1) * el(3) )
     gupyy(i,j,k) = dtm1o3 *( eta(2,2) - two * a2 * hh * el(2) * el(2) )
     gupyz(i,j,k) = dtm1o3 *( eta(2,3) - two * a2 * hh * el(2) * el(3) )
     gupzz(i,j,k) = dtm1o3 *( eta(3,3) - two * a2 * hh * el(3) * el(3) )

!~~~~~~> determine gamx.....
!        \Gamma^i = - \tilde{g^{ij}}_{,j}
! \tilde{g^{ij}}_{,j} = 2 \alpha e^(4\phi) ( 2 H \ell^i \ell^j \alpha_{,j}
!                                            + \alpha \ell^i \ell^j H_{,j}
!                                            + \alpha H \ell^j @_j \ell^i
!                                            + \alpha H \ell^i @_j \ell^j )
!                                          - 4 \tilde{g^{ij}} \phi_{,j}

     LA = el(1) * lapx(i,j,k) + el(2) * lapy(i,j,k) + el(3) * lapz(i,j,k)
     LH = el(1) * phh(1)      + el(2) * phh(2)      + el(3) * phh(3)
     PL =         pel(1,1)    +         pel(2,2)    +         pel(3,3)
     LX = el(1) * pel(1,1)    + el(2) * pel(2,1)    + el(3) * pel(3,1)
     LY = el(1) * pel(1,2)    + el(2) * pel(2,2)    + el(3) * pel(3,2)
     LZ = el(1) * pel(1,3)    + el(2) * pel(2,3)    + el(3) * pel(3,3)

 WX=gupxx(i,j,k)*phix(i,j,k) +gupxy(i,j,k)*phiy(i,j,k) +gupxz(i,j,k)*phiz(i,j,k)
 WY=gupxy(i,j,k)*phix(i,j,k) +gupyy(i,j,k)*phiy(i,j,k) +gupyz(i,j,k)*phiz(i,j,k)
 WZ=gupxz(i,j,k)*phix(i,j,k) +gupyz(i,j,k)*phiy(i,j,k) +gupzz(i,j,k)*phiz(i,j,k)

     gamx(i,j,k) = two*dtm1o3*a1*( two*hh*LA*el(1) + a1*LH*el(1) + &
                                    a1*hh*PL*el(1) + a1*hh*  LX )- fou*WX
     gamy(i,j,k) = two*dtm1o3*a1*( two*hh*LA*el(2) + a1*LH*el(2) + &
                                    a1*hh*PL*el(2) + a1*hh*  LY )- fou*WY
     gamz(i,j,k) = two*dtm1o3*a1*( two*hh*LA*el(3) + a1*LH*el(3) + &
                                    a1*hh*PL*el(3) + a1*hh*  LZ )- fou*WZ

!~~~~~~> determine g_{ij,k}
!        \tilde{g_{ij,k}} = 2 e^(-4\phi) [ \ell_i \ell_j H_{,k}
!                                        + H \ell_i @_k \ell_j
!                                        + H \ell_j @_k \ell_i ]
!                                        -4 \tilde{g_{ij}} \phi_{,k}

    gxxx(i,j,k)= two*idtm1o3*( el(1)*el(1)*phh(1) + hh*el(1)*pel(1,1)      + &
                               hh*el(1)*pel(1,1) )- fou*gxx(i,j,k)*phix(i,j,k)
    gxxy(i,j,k)= two*idtm1o3*( el(1)*el(1)*phh(2) + hh*el(1)*pel(2,1)      + &
                               hh*el(1)*pel(2,1) )- fou*gxx(i,j,k)*phiy(i,j,k)
    gxxz(i,j,k)= two*idtm1o3*( el(1)*el(1)*phh(3) + hh*el(1)*pel(3,1)      + &
                               hh*el(1)*pel(3,1) )- fou*gxx(i,j,k)*phiz(i,j,k)
  
    gxyx(i,j,k)= two*idtm1o3*( el(1)*el(2)*phh(1) + hh*el(1)*pel(1,2)      + &
                               hh*el(2)*pel(1,1) )- fou*gxy(i,j,k)*phix(i,j,k)
    gxyy(i,j,k)= two*idtm1o3*( el(1)*el(2)*phh(2) + hh*el(1)*pel(2,2)      + &
                               hh*el(2)*pel(2,1) )- fou*gxy(i,j,k)*phiy(i,j,k)
    gxyz(i,j,k)= two*idtm1o3*( el(1)*el(2)*phh(3) + hh*el(1)*pel(3,2)      + &
                               hh*el(2)*pel(3,1) )- fou*gxy(i,j,k)*phiz(i,j,k)
    
    gxzx(i,j,k)= two*idtm1o3*( el(1)*el(3)*phh(1) + hh*el(1)*pel(1,3)      + &
                               hh*el(3)*pel(1,1) )- fou*gxz(i,j,k)*phix(i,j,k)
    gxzy(i,j,k)= two*idtm1o3*( el(1)*el(3)*phh(2) + hh*el(1)*pel(2,3)      + &
                               hh*el(3)*pel(2,1) )- fou*gxz(i,j,k)*phiy(i,j,k)
    gxzz(i,j,k)= two*idtm1o3*( el(1)*el(3)*phh(3) + hh*el(1)*pel(3,3)      + &
                               hh*el(3)*pel(3,1) )- fou*gxz(i,j,k)*phiz(i,j,k)
  
    gyyx(i,j,k)= two*idtm1o3*( el(2)*el(2)*phh(1) + hh*el(2)*pel(1,2)      + &
                               hh*el(2)*pel(1,2) )- fou*gyy(i,j,k)*phix(i,j,k)
    gyyy(i,j,k)= two*idtm1o3*( el(2)*el(2)*phh(2) + hh*el(2)*pel(2,2)      + &
                               hh*el(2)*pel(2,2) )- fou*gyy(i,j,k)*phiy(i,j,k)
    gyyz(i,j,k)= two*idtm1o3*( el(2)*el(2)*phh(3) + hh*el(2)*pel(3,2)      + &
                               hh*el(2)*pel(3,2) )- fou*gyy(i,j,k)*phiz(i,j,k)

    gyzx(i,j,k)= two*idtm1o3*( el(2)*el(3)*phh(1) + hh*el(2)*pel(1,3)      + &
                               hh*el(3)*pel(1,2) )- fou*gyz(i,j,k)*phix(i,j,k)
    gyzy(i,j,k)= two*idtm1o3*( el(2)*el(3)*phh(2) + hh*el(2)*pel(2,3)      + &
                               hh*el(3)*pel(2,2) )- fou*gyz(i,j,k)*phiy(i,j,k)
    gyzz(i,j,k)= two*idtm1o3*( el(2)*el(3)*phh(3) + hh*el(2)*pel(3,3)      + &
                               hh*el(3)*pel(3,2) )- fou*gyz(i,j,k)*phiz(i,j,k)
  
    gzzx(i,j,k)= two*idtm1o3*( el(3)*el(3)*phh(1) + hh*el(3)*pel(1,3)      + &
                               hh*el(3)*pel(1,3) )- fou*gzz(i,j,k)*phix(i,j,k)
    gzzy(i,j,k)= two*idtm1o3*( el(3)*el(3)*phh(2) + hh*el(3)*pel(2,3)      + &
                               hh*el(3)*pel(2,3) )- fou*gzz(i,j,k)*phiy(i,j,k)
    gzzz(i,j,k)= two*idtm1o3*( el(3)*el(3)*phh(3) + hh*el(3)*pel(3,3)      + &
                               hh*el(3)*pel(3,3) )- fou*gzz(i,j,k)*phiz(i,j,k)
  
    end do
   end do
  end do
 
  return

  end subroutine convert_a

!-----------------------------
! analytical Gamma
!-----------------------------
  subroutine anagam(ex,                                  &
                       gxx,    gxy,    gxz,    gyy,    gyz,    gzz, &
                     gupxx,  gupxy,  gupxz,  gupyy,  gupyz,  gupzz, &
                    Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                    Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                    Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
                       gxxx,   gxyx,   gxzx,   gyyx,   gyzx,   gzzx,&
                       gxxy,   gxyy,   gxzy,   gyyy,   gyzy,   gzzy,&
                       gxxz,   gxyz,   gxzz,   gyyz,   gyzz,   gzzz )
  implicit none
!
! Input parameters:
!
  integer, dimension(3), intent(in) :: ex
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gxxx,gxyx,gxzx
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gxxy,gxyy,gxzy
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gxxz,gxyz,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gyyz,gyzz,gzzz
!
! Other variables
!
  real*8                             :: SYM, ANTI
  real*8                             :: ZERO, HALF, TWO
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT
  PARAMETER ( SYM = 1.0D0, ANTI = - 1.0D0, ZERO = 0.D0, HALF = 0.5D0 )
  PARAMETER ( TWO = 2.0D0 )
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2)

  Gamxxx =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz &
          - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / Gamxxx
  gupxy = - ( gxy * gzz - gyz * gxz ) / Gamxxx
  gupxz =   ( gxy * gyz - gyy * gxz ) / Gamxxx
  gupyy =   ( gxx * gzz - gxz * gxz ) / Gamxxx
  gupyz = - ( gxx * gyz - gxy * gxz ) / Gamxxx
  gupzz =   ( gxx * gyy - gxy * gxy ) / Gamxxx

!
! Get Connection coefficients
!
  Gamxxx =HALF*(gupxx *gxxx +gupxy *(TWO*gxyx -gxxy )+gupxz *(TWO*gxzx -gxxz ))
  Gamyxx =HALF*(gupxy *gxxx +gupyy *(TWO*gxyx -gxxy )+gupyz *(TWO*gxzx -gxxz ))
  Gamzxx =HALF*(gupxz *gxxx +gupyz *(TWO*gxyx -gxxy )+gupzz *(TWO*gxzx -gxxz ))

  Gamxyy =HALF*(gupxx *(TWO*gxyy -gyyx )+gupxy *gyyy +gupxz *(TWO*gyzy -gyyz ))
  Gamyyy =HALF*(gupxy *(TWO*gxyy -gyyx )+gupyy *gyyy +gupyz *(TWO*gyzy -gyyz ))
  Gamzyy =HALF*(gupxz *(TWO*gxyy -gyyx )+gupyz *gyyy +gupzz *(TWO*gyzy -gyyz ))

  Gamxzz =HALF*(gupxx *(TWO*gxzz -gzzx )+gupxy *(TWO*gyzz -gzzy )+gupxz *gzzz )
  Gamyzz =HALF*(gupxy *(TWO*gxzz -gzzx )+gupyy *(TWO*gyzz -gzzy )+gupyz *gzzz )
  Gamzzz =HALF*(gupxz *(TWO*gxzz -gzzx )+gupyz *(TWO*gyzz -gzzy )+gupzz *gzzz )

  Gamxxy =HALF*( gupxx * gxxy + gupxy * gyyx + gupxz *( gxzy + gyzx - gxyz )  )
  Gamyxy =HALF*( gupxy * gxxy + gupyy * gyyx + gupyz *( gxzy + gyzx - gxyz )  )
  Gamzxy =HALF*( gupxz * gxxy + gupyz * gyyx + gupzz *( gxzy + gyzx - gxyz )  )

  Gamxxz =HALF*( gupxx * gxxz + gupxy *( gxyz + gyzx - gxzy ) + gupxz * gzzx  )
  Gamyxz =HALF*( gupxy * gxxz + gupyy *( gxyz + gyzx - gxzy ) + gupyz * gzzx  )
  Gamzxz =HALF*( gupxz * gxxz + gupyz *( gxyz + gyzx - gxzy ) + gupzz * gzzx  )

  Gamxyz =HALF*( gupxx *( gxyz + gxzy - gyzx ) + gupxy * gyyz + gupxz * gzzy  )
  Gamyyz =HALF*( gupxy *( gxyz + gxzy - gyzx ) + gupyy * gyyz + gupyz * gzzy  )
  Gamzyz =HALF*( gupxz *( gxyz + gxzy - gyzx ) + gupyz * gyyz + gupzz * gzzy  )

  return

  end subroutine anagam

!-----------------------------------------------------------------------------
!
! Find non-boosted rotating Kerr initial data in quasi-isotropic coordinates
!
!-----------------------------------------------------------------------------
  subroutine kerr_initial_metric_quasiisotropic(ex, x, y, z, Rp,             &
       gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
       Axx, Axy, Axz, Ayy, Ayz, Azz, alpha, Bx,By,Bz, &
       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
       phi, sam, vx_boost)
  implicit none

!
! Input parameters:
!
  integer, dimension(3), intent(in)                  :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: x,y,z
  real*8,  intent(in)                                :: sam, vx_boost
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: trK,alpha,Bx,By,Bz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: phi
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: Rp
!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: k1,k2,k3,aa,hh,r2,w2,W,wx,wy,wz,tmp,tmq,tmr
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: f1o3=1.d0/3.d0,m0=1.d0,inbd=1.d0,hlf=0.5d0,fou=4.d0

  real*8 :: grr, grph, gthth, gphph, alp,Bph, dr_driso_alp
  real*8 :: guprr, guprph, gupthth, gupphph,psi4
  real*8 :: Arr, Arth, Arph, Athth, Athph, Aphph
  real*8 :: Krr, Krth,Krph,Kthth,Kthph,Kphph
  real*8 :: Sigma, r, sin2th, sinth, cos2th, costh, lorentz_gamma
  real*8 :: Lx_r, Ly_r, Lz_r, Lx_th, Ly_th, Lz_th, Lx_ph, Ly_ph, Lz_ph
  real*8 :: Lr_x, Lth_x, Lph_x, Lr_y, Lth_y, Lph_y, Lr_z, Lth_z, Lph_z
  real*8 :: r_l, x_l, y_l, z_l, Delta, dtdr, dpdr, rh_iso
  real*8 :: grrtemp, grphtemp, gththtemp, gphphtemp, gupxxtemp
  real*8 :: gupxxl, gupxyl, gupxzl, gupyyl, gupyzl, gupzzl, detg_l
  real*8 :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl
  real*8 :: Kxxl,Kxyl,Kxzl,Kyyl,Kyzl,Kzzl
  real*8  :: rh_out, rh_in, a, riso

!
! BL radius of the outer and inner horizons
!
  a = m0*sam
  rh_out = m0 + sqrt(m0*m0 - a*a)
  rh_in  = m0 - sqrt(m0*m0 - a*a)
  rh_iso = 0.5d0*sqrt(m0*m0 - a*a)

! Input translation
!
  imin = lbound(gxx,1)
  jmin = lbound(gxx,2)
  kmin = lbound(gxx,3)
  imax = ubound(gxx,1)
  jmax = ubound(gxx,2)
  kmax = ubound(gxx,3)

  lorentz_gamma = 1.d0/sqrt(1.d0-vx_boost*vx_boost)

!
! go to each gridpoint...
!
  do k = kmin, kmax
   do j = jmin, jmax
    do i = imin, imax

       riso = sqrt( (lorentz_gamma*X(i,j,k))**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
       z_l = Z(i,j,k)
       r = riso + m0 + 0.25d0*(m0*m0 - a*a)/riso
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
       alpha(i,j,k) = alp
        Bph = -2.d0*m0*a*r/aa

!-- the 3-metric (nonzero components)
!
       grr = Sigma/(r-rh_in)*( (riso+ rh_iso)**2 ) / riso**3
       grph = 0.d0
       gthth = Sigma
       gphph = aa*sin2th/Sigma
!
!--  K_{ij} in "isotropic" Boyer-Lindquist coordinates
!--         
!
       !!dr_driso_alp = sqrt(aa*riso/Sigma/(r-rh_in))*(riso+rh_iso)/(riso*riso)
       Krr = 0.d0
       Krth = 0.d0
       Krph = m0*a*sin2th/( Sigma*riso*riso*sqrt(aa*Sigma) ) * & 
              (3.d0*r**4 + 2.d0*a*a*r*r - a**4 - a*a*(r*r-a*a)*sin2th) * & 
              (riso + rh_iso)*sqrt( riso/(r-rh_in) )
       !!Krph = (m0*a*sin2th*(-a**4 + 2.d0*a*a*r2 + 3.d0*r2*r2 + &
       !!     a*a*(a - r)*(a + r)*sin2th))/Sigma/aa * dr_driso_alp
       Kthth = 0.d0
       Kthph = -2.d0*a*a*a*m0*r*costh*sinth*sin2th/(Sigma*sqrt(aa*Sigma)) * & 
               (riso - rh_iso)*sqrt( (r-rh_in)/riso )
       !!Kthph = -2.d0*a*a*a*r*costh*sinth*sin2th/Sigma * sqrt(Delta/Sigma/aa)
       !!if (riso .lt. rh_iso) Kthph = -Kthph
       Kphph = 0.d0
!
!-- the scalar extrinsic curvature
!
       trK(i,j,k) = 0.d0
!
!-- A_{ij} = K_{ij} - \gamma_{ij}K/3 
!--         
!
       Arr = Krr
       Arth = Krth
       Arph = Krph
       Athth = Kthth
       Athph = Kthph
       Aphph = Kphph

     !#############################################################################
     !Now convert everything to Cartesian
     !#############################################################################
       r_l = sqrt( (lorentz_gamma*x(i,j,k))**2 + y(i,j,k)**2 + z(i,j,k)**2)
       r2 = r_l * r_l

       x_l = x(i,j,k)*lorentz_gamma
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

       Bx(i,j,k) = Lx_ph*Bph
       By(i,j,k) = Ly_ph*Bph
       Bz(i,j,k) = 0.d0

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
       Axx(i,j,k) = Lr_x**2 * Arr + &
		    2.d0*Lr_x*Lth_x * Arth + &
                    2.d0*Lr_x*Lph_x * Arph + &
                    Lth_x**2 * Athth + &
		    2.d0*Lth_x*Lph_x * Athph + &
                    Lph_x**2 * Aphph

       Axy(i,j,k) = Lr_x*Lr_y * Arr + &
		    (Lr_x*Lth_y + Lr_y*Lth_x) * Arth + &
                    (Lr_x*Lph_y + Lr_y*Lph_x) * Arph + &
                    Lth_x*Lth_y * Athth + &
		    (Lth_x*Lph_y + Lth_y*Lph_x) * Athph + &
                    Lph_x*Lph_y * Aphph

       Axz(i,j,k) = Lr_x*Lr_z * Arr + &
		    (Lr_x*Lth_z + Lr_z*Lth_x) * Arth + &
                    (Lr_x*Lph_z + Lr_z*Lph_x) * Arph + &
                    Lth_x*Lth_z * Athth + &
		    (Lth_x*Lph_z + Lth_z*Lph_x) * Athph + &
                    Lph_x*Lph_z * Aphph

       Ayy(i,j,k) = Lr_y**2 * Arr + &
		    2.d0*Lr_y*Lth_y * Arth + &
                    2.d0*Lr_y*Lph_y * Arph + &
                    Lth_y**2 * Athth + &
		    2.d0*Lth_y*Lph_y * Athph + &
                    Lph_y**2 * Aphph

       Ayz(i,j,k) = Lr_y*Lr_z * Arr + &
		    (Lr_y*Lth_z + Lr_z*Lth_y) * Arth + &
                    (Lr_y*Lph_z + Lr_z*Lph_y) * Arph + &
                    Lth_y*Lth_z * Athth + &
		    (Lth_y*Lph_z + Lth_z*Lph_y) * Athph + &
                    Lph_y*Lph_z * Aphph

       Azz(i,j,k) = Lr_z**2 * Arr + &
		    2.d0*Lr_z*Lth_z * Arth + &
                    2.d0*Lr_z*Lph_z * Arph + &
                    Lth_z**2 * Athth + &
		    2.d0*Lth_z*Lph_z * Athph + &
                    Lph_z**2 * Aphph

       detg_l = gxxl * gyyl * gzzl + &
                gxyl * gyzl * gxzl + &
                gxzl * gxyl * gyzl &
                - gxzl * gyyl * gxzl &
                - gxyl * gxyl * gzzl &
                - gxxl * gyzl * gyzl

       phi(i,j,k) = log(detg_l)/12.d0
       psi4 = detg_l**(1.d0/3.d0)
       gxx(i,j,k) = gxxl/psi4
       gxy(i,j,k) = gxyl/psi4
       gxz(i,j,k) = gxzl/psi4
       gyy(i,j,k) = gyyl/psi4
       gyz(i,j,k) = gyzl/psi4
       gzz(i,j,k) = gzzl/psi4
       Axx(i,j,k) = Axx(i,j,k)/psi4
       Axy(i,j,k) = Axy(i,j,k)/psi4
       Axz(i,j,k) = Axz(i,j,k)/psi4
       Ayy(i,j,k) = Ayy(i,j,k)/psi4
       Ayz(i,j,k) = Ayz(i,j,k)/psi4
       Azz(i,j,k) = Azz(i,j,k)/psi4

    end do
   end do
  end do

  alpha = alpha - one

  gupxx =   ( gyy * gzz - gyz * gyz )
  gupxy = - ( gxy * gzz - gyz * gxz )
  gupxz =   ( gxy * gyz - gyy * gxz )
  gupyy =   ( gxx * gzz - gxz * gxz )
  gupyz = - ( gxx * gyz - gxy * gxz )
  gupzz =   ( gxx * gyy - gxy * gxy )

  return

end subroutine kerr_initial_metric_quasiisotropic
