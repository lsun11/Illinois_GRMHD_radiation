!-----------------------------------------------------------------------------
!
! Find non-boosted rotating Kerr-Schild metric initial data
! This is for the case when were using the Boyer-Lindquist radius
!
!-----------------------------------------------------------------------------
  subroutine ks_initial_metric_shift_origin_bl(ex, x, y, z, &
       gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
       Axx, Axy, Axz, Ayy, Ayz, Azz, alpha, Bx,By,Bz, &
       sam, r0, Symmetry)
  implicit none
!
! Input parameters:
!
  integer, dimension(3), intent(in)                  :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in)  :: x,y,z
  real*8,  intent(in)                                :: sam,r0
  integer, intent(in)                                :: Symmetry
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: trK,alpha,Bx,By,Bz
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
  real*8 :: gupxxl, gupxyl, gupxzl, gupyyl, gupyzl, gupzzl, detg_l
  integer :: itest,jtest,ktest
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
    r = rgrid_l + r0
    r2 = r * r
    costh = z(i,j,k)/rgrid_l
    cos2th = costh*costh
    sin2th = 1.d0 - cos2th
    sinth = sqrt(sin2th)
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
     !#############################################################################
     !Now convert everything to Cartesian
     !#############################################################################
       r_l = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
       r2 = r_l * r_l
       x_l = x(i,j,k)
       y_l = y(i,j,k)
       z_l = z(i,j,k)
       Lx_r = x_l/r_l
       Ly_r = y_l/r_l
       Lz_r = z_l/r_l
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
    end do
   end do
  end do
  alpha = alpha - one
  return
end subroutine ks_initial_metric_shift_origin_bl
