!-----------------------------------------------------------------------------
!
! Find non-boosted rotating Kerr-Schild metric initial data
! This is for the case when we're using the Boyer-Lindquist radius
!
!-----------------------------------------------------------------------------
  subroutine ks_initial_metric_bl(ex, x, y, z, Rp, dRdr,             &
       gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
       Axx, Axy, Axz, Ayy, Ayz, Azz, alpha, Bx,By,Bz, &
       sam, Symmetry)
  implicit none

!
! Input parameters:
!
  integer, dimension(3), intent(in)                  :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: x,y,z
  real*8,  intent(in)                                :: sam
  integer, intent(in)                                :: Symmetry
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: trK,alpha,Bx,By,Bz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: Rp,dRdr
!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: k1,k2,k3,aa,hh,r2,w2,W,wx,wy,wz,tmp,tmq,tmr
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: f1o3=1.d0/3.d0,m0=1.d0,inbd=1.d0,hlf=0.5d0,fou=4.d0

  real*8 :: grr, grph, gthth, gphph, Br
  real*8 :: guprr, guprph, gupthth, gupphph
  real*8 :: Arr, Arth, Arph, Athth, Athph, Aphph
  real*8 :: Sigma, r, sin2th, sinth, cos2th, costh
  real*8 :: Lx_r, Ly_r, Lz_r, Lx_th, Ly_th, Lz_th, Lx_ph, Ly_ph, Lz_ph
  real*8 :: Lr_x, Lth_x, Lph_x, Lr_y, Lth_y, Lph_y, Lr_z, Lth_z, Lph_z
  real*8 :: r_l, x_l, y_l, z_l, rgrid_l
  real*8 :: grrtemp, grphtemp, gththtemp, gphphtemp, gupxxtemp
  logical :: fish_to_phys
  real*8 :: gupxxl, gupxyl, gupxzl, gupyyl, gupyzl, gupzzl, detg_l
  integer :: itest,jtest,ktest


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
! go to each gridpoint...
!
  do k = kmin, kmax
   do j = jmin, jmax
    do i = imin, imax

    rgrid_l = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
    r = Rp(i,j,k)
    r2 = r * r
    z_l = z(i,j,k)*r/rgrid_l
    sin2th = one - z(i,j,k)**2/rgrid_l**2
    sinth = sqrt(sin2th)
    costh = z_l/r
    costh = z(i,j,k)/rgrid_l
    cos2th = costh**2
    Sigma = r2 + sam*sam*cos2th
!
!-- the lapse and shift (nonzero component)
!
     alpha(i,j,k) = 1.d0/sqrt(one+two*r/Sigma)
        Br = two*r/(Sigma + two*r)
!
!-- the 3-metric (nonzero components)
!
     grr = one+two*r/Sigma
     grph = -sam*sin2th*grr 
     gthth = Sigma
     gphph = sin2th*(Sigma + sam**2*sin2th*grr)

     guprr = (Sigma + sam**2*sin2th*grr)/(grr*Sigma)
     guprph = sam/Sigma
     gupthth = 1.d0/Sigma
     gupphph = 1.d0/(sin2th*Sigma)

!
!-- Store K_{ij} in A_{ij} 
!--         
!
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


!
!-- the scalar extrinsic curvature
!
     trK(i,j,k) = 4.d0/(Sigma*(Sigma+two*r)**3)*( &
          (Sigma+sam*sam*sin2th*(one+two*r/Sigma))* &
          (3.d0*r*Sigma + two*r2 + Sigma*Sigma)* &
          (one - two*r2/Sigma))/(two*alpha(i,j,k))
     trK(i,j,k) = trK(i,j,k) - 4.d0*sam*sam*sin2th/Sigma**2* &
          (one - two*r2/Sigma)/(two*alpha(i,j,k)) 
     trK(i,j,k) = trK(i,j,k) + 4.d0*r/(Sigma*(Sigma+two*r))* &
          (two*r + sam*sam*sin2th*((Sigma-two*r2)/Sigma**2))/(two*alpha(i,j,k))

!!$     trK(i,j,K) = guprr*Arr + two*guprph*Arph + &
!!$                  gupthth*Athth + gupphph*Aphph

!
!-- A_{ij} = K_{ij} - \gamma_{ij}K/3 
!--         
!
     Arr = Arr - f1o3*grr*trK(i,j,k)
     Arph = Arph - f1o3*grph*trK(i,j,k)
     Athth = Athth - f1o3*gthth*trK(i,j,k)
     Aphph = Aphph - f1o3*gphph*trK(i,j,k)

     !finally, transform to the fisheye system
     grr = grr*dRdr(i,j,k)**2
     grph = grph*dRdr(i,j,k)

     guprr = guprr/dRdr(i,j,k)**2
     guprph = guprph/dRdr(i,j,k)

     Arr = Arr*dRdr(i,j,k)**2
     Arth = Arth*dRdr(i,j,k)
     Arph = Arph*dRdr(i,j,k)

     Br = Br/dRdr(i,j,k)

     !#############################################################################
     !Now convert everything to Cartesian
     !#############################################################################

       r_l = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
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
       gxx(i,j,k) = Lr_x**2 * grr + &
                    2.d0*Lr_x*Lph_x * grph + &
                    Lth_x**2 * gthth + &
                    Lph_x**2 * gphph

       gxy(i,j,k) = Lr_x*Lr_y * grr + &
                    (Lr_x*Lph_y + Lr_y*Lph_x) * grph + &
                    Lth_x*Lth_y * gthth + &
                    Lph_x*Lph_y * gphph

       gxz(i,j,k) = Lr_x*Lr_z * grr + &
                    (Lr_x*Lph_z + Lr_z*Lph_x) * grph + &
                    Lth_x*Lth_z * gthth + &
                    Lph_x*Lph_z * gphph

       gyy(i,j,k) = Lr_y**2 * grr + &
                    2.d0*Lr_y*Lph_y * grph + &
                    Lth_y**2 * gthth + &
                    Lph_y**2 * gphph

       gyz(i,j,k) = Lr_y*Lr_z * grr + &
                    (Lr_y*Lph_z + Lr_z*Lph_y) * grph + &
                    Lth_y*Lth_z * gthth + &
                    Lph_y*Lph_z * gphph

       gzz(i,j,k) = Lr_z**2 * grr + &
                    2.d0*Lr_z*Lph_z * grph + &
                    Lth_z**2 * gthth + &
                    Lph_z**2 * gphph

       ! Transform the extrinsic curvature
       Axx(i,j,k) = Lr_x**2 * Arr + &
                    2.d0*Lr_x*Lph_x * Arph + &
                    Lth_x**2 * Athth + &
                    Lph_x**2 * Aphph

       Axy(i,j,k) = Lr_x*Lr_y * Arr + &
                    (Lr_x*Lph_y + Lr_y*Lph_x) * Arph + &
                    Lth_x*Lth_y * Athth + &
                    Lph_x*Lph_y * Aphph

       Axz(i,j,k) = Lr_x*Lr_z * Arr + &
                    (Lr_x*Lph_z + Lr_z*Lph_x) * Arph + &
                    Lth_x*Lth_z * Athth + &
                    Lph_x*Lph_z * Aphph

       Ayy(i,j,k) = Lr_y**2 * Arr + &
                    2.d0*Lr_y*Lph_y * Arph + &
                    Lth_y**2 * Athth + &
                    Lph_y**2 * Aphph

       Ayz(i,j,k) = Lr_y*Lr_z * Arr + &
                    (Lr_y*Lph_z + Lr_z*Lph_y) * Arph + &
                    Lth_y*Lth_z * Athth + &
                    Lph_y*Lph_z * Aphph

       Azz(i,j,k) = Lr_z**2 * Arr + &
                    2.d0*Lr_z*Lph_z * Arph + &
                    Lth_z**2 * Athth + &
                    Lph_z**2 * Aphph

    end do
   end do
  end do

  alpha = alpha - one

  return

end subroutine ks_initial_metric_bl


!-----------------------------------------------------------------------------
!
! Find non-boosted rotating Kerr-Schild metric initial data (fisheye version)
!
!-----------------------------------------------------------------------------
  subroutine ks_initial_metric_fisheye(ex, x, y, z, Rp, dRdr,              &
                    gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
                    Axx, Axy, Axz, Ayy, Ayz, Azz, alpha, Bx,By,Bz, &
                    sam, Symmetry)
  implicit none

!
! Input parameters:
!
  integer, dimension(3),     intent(in) :: ex
  real*8,                    intent(in) :: sam
  integer, intent(in)                   :: Symmetry
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: Rp, dRdr, x, y, z
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: trK,alpha,Bx,By,Bz
  real*8,  dimension(ex(1),ex(2),ex(3))              :: det
  real*8                                             :: detmin, detmax
!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: el(0:3),boo(0:3),phh(0:3),eta(0:3,0:3),pel(0:3,0:3)
  real*8  :: k1,k2,k3,aa,hh,r2,w2,W,wx,wy,wz,tmp,tmq,tmr
  real*8  :: rp2, xp, yp, zp
!  real*8,  dimension(ex(1),ex(2),ex(3)) :: qxx,qxy,qxz,qyy,qyz,qzz,q
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: f1o3=1.d0/3.d0,m0=1.d0,inbd=1.d0,hlf=0.5d0,fou=4.d0
  logical           :: fish_to_phys
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
!    Note, no longer using the radius function since this 
!    will give incorrect results off of the plane in Axisymmetry.
!    r2 = radius(x(i,j,k),y(i,j,k),z(i,j,k),1)
!    r2 = r2 * r2
    r2 = x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2
    rp2 = Rp(i,j,k)*Rp(i,j,k)
    xp = x(i,j,k)*sqrt(rp2/r2)
    yp = y(i,j,k)*sqrt(rp2/r2)
    zp = z(i,j,k)*sqrt(rp2/r2)

!-- W^2 = 1/2 *[ r^2 - a^2 + \sqrt(( r^2 - a^2 )^2 + 4 a^2 z^2 )]

    w2 = rp2 - sam * sam
    w2 = hlf *( w2 + sqrt( w2 * w2 + fou * sam*sam*zp*zp ))
     W = sqrt(w2)

    tmp = two * w2 + sam * sam - rp2
    tmq =       w2 + sam * sam

     wx =   W * xp / tmp
     wy =   W * yp / tmp
     wz = tmq * zp /(tmp * W)

!-- H = M W^3 /( w^4 + a^2 z^2 )

    hh = m0 * w2 * W /( w2 * w2 + sam*sam*zp*zp )

!
!-- the null vector
!
     el(0) = one
     el(1) = ( W * xp + sam * yp )/ tmq
     el(2) = ( W * yp - sam * xp )/ tmq
     el(3) =       zp               / W

!
!-- @_\mu H
!

     tmp = w2 * w2 + sam * sam * zp * zp
     tmr = thr * sam * sam * zp * zp - w2 * w2

     phh(0) = zeo
     phh(1) = hh * tmr * wx /( W * tmp )
     phh(2) = hh * tmr * wy /( W * tmp )
     phh(3) = hh *(tmr * wz /  W - two * sam * sam * zp )/ tmp

!
!-- @_\mu L_\nu
!
       pel = zeo

       pel(1,1) = (     W + ( xp - two * W * el(1) )* wx )/ tmq
       pel(2,1) = (   sam + ( xp - two * W * el(1) )* wy )/ tmq
       pel(3,1) =           ( xp - two * W * el(1) )* wz  / tmq

       pel(1,2) = ( - sam + ( yp - two * W * el(2) )* wx )/ tmq
       pel(2,2) = (     W + ( yp - two * W * el(2) )* wy )/ tmq
       pel(3,2) =           ( yp - two * W * el(2) )* wz  / tmq

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

     ! Note to fisheye users:  this is the physical A_ij, not the conformal A_ij

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

end subroutine ks_initial_metric_fisheye
