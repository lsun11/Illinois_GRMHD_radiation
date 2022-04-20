!-----------------------------------------------------------------------------
!
! $Id: ks_initial_metric.f90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!
! Find non-boosted rotating Kerr-Schild metric initial data
!
!-----------------------------------------------------------------------------
  subroutine ks_initial_metric(ex, x, y, z,                        &
                    gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
                    Axx, Axy, Axz, Ayy, Ayz, Azz, alpha, Bx,By,Bz, &
                    sam, Symmetry)
  implicit none

!
! Input parameters:
!
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z
  real*8,                    intent(in) :: sam
  integer, intent(in)                   :: Symmetry
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

!~~~~~> Interface

  interface
    real*8 function radius(x,y,z,Symmetry)
      implicit none
      real*8,  intent(in) :: x,y,z
      integer, intent(in) :: Symmetry
    end function radius
  end interface

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
    r2 = radius(x(i,1,1),y(1,j,1),z(1,1,k),1)
    r2 = r2 * r2

!-- W^2 = 1/2 *[ r^2 - a^2 + \sqrt(( r^2 - a^2 )^2 + 4 a^2 z^2 )]

    w2 = r2 - sam * sam
    w2 = hlf *( w2 + sqrt( w2 * w2 + fou * sam * sam * z(1,1,k) * z(1,1,k) ))
     W = sqrt(w2)

    tmp = two * w2 + sam * sam - r2
    tmq =       w2 + sam * sam

     wx =   W * x(i,1,1) / tmp
     wy =   W * y(1,j,1) / tmp
     wz = tmq * z(1,1,k) /(tmp * W)

!-- H = M W^3 /( w^4 + a^2 z^2 )

    hh = m0 * w2 * W /( w2 * w2 + sam * sam * z(1,1,k) * z(1,1,k) )
 
!
!-- the null vector
!
     el(0) = one
     el(1) = ( W * x(i,1,1) + sam * y(1,j,1) )/ tmq
     el(2) = ( W * y(1,j,1) - sam * x(i,1,1) )/ tmq
     el(3) =       z(1,1,k)               / W

!
!-- @_\mu H
!

     tmp = w2 * w2 + sam * sam * z(1,1,k) * z(1,1,k)
     tmr = thr * sam * sam * z(1,1,k) * z(1,1,k) - w2 * w2

     phh(0) = zeo
     phh(1) = hh * tmr * wx /( W * tmp )
     phh(2) = hh * tmr * wy /( W * tmp )
     phh(3) = hh *(tmr * wz /  W - two * sam * sam * z(1,1,k) )/ tmp

!
!-- @_\mu L_\nu
!
       pel = zeo

       pel(1,1) = (     W + ( x(i,1,1) - two * W * el(1) )* wx )/ tmq
       pel(2,1) = (   sam + ( x(i,1,1) - two * W * el(1) )* wy )/ tmq
       pel(3,1) =           ( x(i,1,1) - two * W * el(1) )* wz  / tmq

       pel(1,2) = ( - sam + ( y(1,j,1) - two * W * el(2) )* wx )/ tmq
       pel(2,2) = (     W + ( y(1,j,1) - two * W * el(2) )* wy )/ tmq
       pel(3,2) =           ( y(1,j,1) - two * W * el(2) )* wz  / tmq

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

end subroutine ks_initial_metric
!-----------------------------------------------------------------------------
!
! Convert untilded metric to tilded metric with analytical value 
!
!-----------------------------------------------------------------------------
  subroutine convert_a(ex,    x,     y,     z,     phi,          &
                       gxx,   gxy,   gxz,   gyy,   gyz,   gzz,   & 
                       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                       gxxx,  gxxy,  gxxz,  gxyx,  gxyy,  gxyz,  &
                       gxzx,  gxzy,  gxzz,  gyyx,  gyyy,  gyyz,  &
                       gyzx,  gyzy,  gyzz,  gzzx,  gzzy,  gzzz,  &
                       gamx,  gamy,  gamz,  phix,  phiy,  phiz,  &
                       lapx,  lapy,  lapz, Symmetry)
  implicit none

!~~~~~> input variables

  integer, intent(in) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in)  :: x,y,z
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
  integer, intent(in)                                :: Symmetry

!~~~~~> local variable

  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,m,n
  real*8  :: a2, hh, rr, dtm, idtm1o3, dtm1o3, a1
  real*8  :: el(0:3),boo(0:3),phh(0:3),eta(0:3,0:3),pel(0:3,0:3)
  real*8  :: LA,LH,LX,LY,LZ,PL,WX,WY,WZ
  real*8, parameter :: f1o3=1.d0/3.d0, f1o12=1.d0/12.d0, m0=1.d0, thr=3.d0
  real*8, parameter :: zeo= 0.d0, one=1.d0, two=2.d0, fou=4.d0, f4o3=4.d0/3.d0
  real*8, parameter :: f1o6=1.d0/6.d0
 
!~~~~~> Interface

  interface
   real*8 function radius(x,y,z,Symmetry)
    implicit none
    real*8,  intent(in) :: x,y,z
    integer, intent(in) :: Symmetry
   end function 
  end interface

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

     rr = radius(x(i,1,1),y(1,j,1),z(1,1,k),1)
     hh = m0 / rr

!~~~~~~> the null vector
!   \ell_t = \gamma -\gamma^2 (v_x x + v_y y) / r
!   \ell_x = x / r - v_x \ell_t
!   \ell_y = y / r - v_y \ell_t
!   \ell_z = z / r

     el(0) = one
     el(1) = x(i,1,1) / rr
     el(2) = y(1,j,1) / rr
     el(3) = z(1,1,k) / rr

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
!---------------------------------------------------------------------
!
! radius
!
!---------------------------------------------------------------------

  real*8 function radius(x,y,z,Symmetry)
    implicit none
    real*8,  intent(in) :: x,y,z
    integer, intent(in) :: Symmetry
    integer, parameter  :: AXISYM = 4

    if (Symmetry == AXISYM) then
       radius = sqrt(x*x + z*z)
    else
       radius = sqrt(x*x + y*y + z*z)
    end if
  end function radius
