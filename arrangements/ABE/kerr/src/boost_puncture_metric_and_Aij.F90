subroutine compute_boost_metric_Aij(m0,sam,x0_bh,vx_boost,ex,x,y,z, &
 			lapm1,betax,betay,betaz,phi,gxx,gxy,gxz, & 
                        gyy,gyz,gzz,gupxx,gupxy,gupxz,gupyy, & 
                        gupyz,gupzz,Axx,Axy,Axz,Ayy,Ayz,Azz,trK)
  implicit none
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z
  real*8,                    intent(in) :: m0,sam,x0_bh,vx_boost
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: trK,lapm1,phi,betax,betay,betaz
  integer :: i,j,k,imin,jmin,kmin,imax,jmax,kmax
  real*8 :: xl,yl,zl,f1o3,f1o2alp,riso2,lorentz_gamma,rh,rh2
  real*8 :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,detg,psi4,psim4
  real*8 :: gupxxl,gupxyl,gupxzl,gupyyl,gupyzl,gupzzl
  real*8 :: alpl,betaxl,betayl,betazl, a
  real*8 :: Kxxl,Kxyl,Kxzl,Kyyl,Kyzl,Kzzl,trKl
  real*8 :: Axxl,Axyl,Axzl,Ayyl,Ayzl,Azzl
  real*8 :: gxxt,gxyt,gxzt,gyyt,gyzt,gzzt
  real*8 :: beta_xx, beta_xy, beta_xz, beta_yy, beta_yz, beta_zz
!
  imin = lbound(gxx,1)
  jmin = lbound(gxx,2)
  kmin = lbound(gxx,3)
  imax = ubound(gxx,1)
  jmax = ubound(gxx,2)
  kmax = ubound(gxx,3)

  a = m0*sam 
  f1o3 = 1.d0/3.d0
  lorentz_gamma = 1.d0/sqrt(1.d0-vx_boost**2)
  ! Horizon radius in "isotropic" radial coordinate
  rh = 0.25d0*( m0 + sqrt(m0*m0 - a*a) )
  rh2 = rh*rh

  do k=kmin,kmax
     do j=jmin,jmax
        do i=imin,imax
	   xl = x(i,j,k) - x0_bh
	   yl = y(i,j,k)
	   zl = z(i,j,k)
           call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,xl,yl,zl, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, & 
                        gupxxl,gupxyl,gupxzl,gupyyl,gupyzl,gupzzl)
           lapm1(i,j,k) = alpl - 1.d0
	   betax(i,j,k) = betaxl
           betay(i,j,k) = betayl
           betaz(i,j,k) = betazl
    	   call dt_boost_gammaij(m0,a,vx_boost,xl,yl,zl,gxxt,gxyt,gxzt,gyyt,gyzt,gzzt)
           call cov_derv_beta(m0,a,vx_boost,xl,yl,zl,betaxl,betayl,betazl, &
                        beta_xx, beta_xy, beta_xz, beta_yy, beta_yz, beta_zz)

	   ! Compute Kij
           f1o2alp = 0.5d0/alpl
	   riso2 = (lorentz_gamma*xl)**2 + yl*yl + zl*zl
	   ! Note: we need to flip the sign of Kij inside the horizon
	   if (riso2 .lt. rh2) f1o2alp = -f1o2alp 
	   Kxxl = (beta_xx - gxxt)*f1o2alp
           Kxyl = (beta_xy - gxyt)*f1o2alp
           Kxzl = (beta_xz - gxzt)*f1o2alp
           Kyyl = (beta_yy - gyyt)*f1o2alp
           Kyzl = (beta_yz - gyzt)*f1o2alp
           Kzzl = (beta_zz - gzzt)*f1o2alp

	   trKl = gupxxl*Kxxl + 2.d0*gupxyl*Kxyl + 2.d0*gupxzl*Kxzl & 
		  + gupyyl*Kyyl + 2.d0*gupyzl*Kyzl + gupzzl*Kzzl
	   trK(i,j,k) = trKl
	   Axxl = Kxxl - f1o3*gxxl*trKl
           Axyl = Kxyl - f1o3*gxyl*trKl
           Axzl = Kxzl - f1o3*gxzl*trKl
           Ayyl = Kyyl - f1o3*gyyl*trKl
           Ayzl = Kyzl - f1o3*gyzl*trKl
           Azzl = Kzzl - f1o3*gzzl*trKl

	   ! Convert to "tilde"
	   detg = gxxl * gyyl * gzzl + &
            gxyl * gyzl * gxzl + &
            gxzl * gxyl * gyzl &
            - gxzl * gyyl * gxzl &
            - gxyl * gxyl * gzzl &
            - gxxl * gyzl * gyzl
	   phi(i,j,k) = log(detg)/12.d0
	   psi4 = detg**f1o3
	   psim4 = 1.d0/psi4
	   gxx(i,j,k) = gxxl*psim4
           gxy(i,j,k) = gxyl*psim4
           gxz(i,j,k) = gxzl*psim4
           gyy(i,j,k) = gyyl*psim4
           gyz(i,j,k) = gyzl*psim4
           gzz(i,j,k) = gzzl*psim4
           Axx(i,j,k) = Axxl*psim4
           Axy(i,j,k) = Axyl*psim4
           Axz(i,j,k) = Axzl*psim4
           Ayy(i,j,k) = Ayyl*psim4
           Ayz(i,j,k) = Ayzl*psim4
           Azz(i,j,k) = Azzl*psim4
           gupxx(i,j,k) = gupxxl*psi4
           gupxy(i,j,k) = gupxyl*psi4
           gupxz(i,j,k) = gupxzl*psi4
           gupyy(i,j,k) = gupyyl*psi4
           gupyz(i,j,k) = gupyzl*psi4
           gupzz(i,j,k) = gupzzl*psi4
	end do
     end do
  end do
end subroutine compute_boost_metric_Aij

subroutine dt_boost_gammaij(m0,a,vx_boost,x,y,z,gxxt,gxyt,gxzt,gyyt,gyzt,gzzt)
  implicit none
  real*8,  intent(in)                                :: m0,a,x,y,z,vx_boost
  real*8,  intent(out)                               :: gxxt,gxyt,gxzt,gyyt,gyzt,gzzt
  real*8 :: gamma,gammav,Ltt,Ltx,Lxx,Lxt, delta,r,f1o2d
  real*8 :: gttx0,gtxx0,gtyx0,gtzx0,gxxx0,gxyx0,gxzx0,gyyx0,gyzx0,gzzx0
  real*8 :: gtt0,gtx0,gty0,gtz0,gxx0,gxy0,gxz0,gyy0,gyz0,gzz0
  real*8 :: gtt1,gtx1,gty1,gtz1,gxx1,gxy1,gxz1,gyy1,gyz1,gzz1
  real*8 :: g4uptt0,g4uptx0,g4upty0,g4uptz0,gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0
!
  gamma = 1.d0/sqrt(1.d0-vx_boost*vx_boost)
  gammav = gamma*vx_boost
  ! \partial x^\mu / \partial x_boost^\nu
  Ltx = -gammav
  Lxx = gamma

  ! Compute x-derivative of the unboosted metric numerically 
  r = sqrt(x*x+y*y+z*z)
  delta = 1.d-4*r
  f1o2d = 0.5d0/delta
  call compute_kerr_puncture_metric(m0,a,gamma*x+delta,y,z,gtt1,gtx1,gty1,gtz1, &
                        gxx1,gxy1,gxz1,gyy1,gyz1,gzz1,g4uptt0,g4uptx0,g4upty0,g4uptz0, &
                        gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
  call compute_kerr_puncture_metric(m0,a,gamma*x-delta,y,z,gtt0,gtx0,gty0,gtz0, &
                        gxx0,gxy0,gxz0,gyy0,gyz0,gzz0,g4uptt0,g4uptx0,g4upty0,g4uptz0, &
                        gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
  gttx0 = (gtt1-gtt0)*f1o2d
  gtxx0 = (gtx1-gtx0)*f1o2d
  gtyx0 = (gty1-gty0)*f1o2d
  gtzx0 = (gtz1-gtz0)*f1o2d
  gxxx0 = (gxx1-gxx0)*f1o2d
  gxyx0 = (gxy1-gxy0)*f1o2d
  gxzx0 = (gxz1-gxz0)*f1o2d
  gyyx0 = (gyy1-gyy0)*f1o2d
  gyzx0 = (gyz1-gyz0)*f1o2d
  gzzx0 = (gzz1-gzz0)*f1o2d

  gxxt = -gammav*(Ltx*Ltx*gttx0 + 2.d0*Ltx*Lxx*gtxx0 + Lxx*Lxx*gxxx0)
  gxyt = -gammav*(Ltx*gtyx0 + Lxx*gxyx0)
  gxzt = -gammav*(Ltx*gtzx0 + Lxx*gxzx0)
  gyyt = -gammav*gyyx0
  gyzt = -gammav*gyzx0
  gzzt = -gammav*gzzx0

end subroutine dt_boost_gammaij

! Compute D_i beta_j + D_j beta_i for the boosted metric and shift 
!
subroutine cov_derv_beta(m0,a,vx_boost,x,y,z,betax,betay,betaz, &
                        beta_xx, beta_xy, beta_xz, beta_yy, beta_yz, beta_zz)
  implicit none
  real*8,  intent(in)                                :: m0,a,x,y,z,vx_boost
  real*8,  intent(in)				     :: betax,betay,betaz
  real*8,  intent(out)                               :: beta_xx, beta_xy, beta_xz, beta_yy, beta_yz, beta_zz
  real*8 :: beta_x0,beta_y0,beta_z0,beta_x1,beta_y1,beta_z1
  real*8 :: beta_xxt,beta_xyt,beta_xzt,beta_yxt,beta_yyt,beta_yzt,beta_zxt,beta_zyt,beta_zzt
  real*8 :: delta,r,f1o2d
  real*8 :: gxx0,gxy0,gxz0,gyy0,gyz0,gzz0,gxx1,gxy1,gxz1,gyy1,gyz1,gzz1
  real*8 :: gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0,alp0,betax0,betay0,betaz0
  real*8 :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8 :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8 :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8 :: gamxxx,gamxxy,gamxxz,gamxyy,gamxyz,gamxzz
  real*8 :: gamyxx,gamyxy,gamyxz,gamyyy,gamyyz,gamyzz
  real*8 :: gamzxx,gamzxy,gamzxz,gamzyy,gamzyz,gamzzz
!
  r = sqrt(x*x + y*y + z*z) 
  delta = 1.d-4*r
  f1o2d = 0.5d0/delta

  ! x-derivatives
  call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,x+delta,y,z,alp0,betax0,betay0,betaz0, &
                        gxx1,gxy1,gxz1,gyy1,gyz1,gzz1,gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
  beta_x1 = gxx1*betax0 + gxy1*betay0 + gxz1*betaz0
  beta_y1 = gxy1*betax0 + gyy1*betay0 + gyz1*betaz0
  beta_z1 = gxz1*betax0 + gyz1*betay0 + gzz1*betaz0

  call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,x-delta,y,z,alp0,betax0,betay0,betaz0, &
                        gxx0,gxy0,gxz0,gyy0,gyz0,gzz0,gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
  beta_x0 = gxx0*betax0 + gxy0*betay0 + gxz0*betaz0
  beta_y0 = gxy0*betax0 + gyy0*betay0 + gyz0*betaz0
  beta_z0 = gxz0*betax0 + gyz0*betay0 + gzz0*betaz0

  ! partial_x beta_i and partial_x g_ij
  beta_xxt = (beta_x1 - beta_x0)*f1o2d
  beta_yxt = (beta_y1 - beta_y0)*f1o2d
  beta_zxt = (beta_z1 - beta_z0)*f1o2d
  gxxx = (gxx1 - gxx0)*f1o2d
  gxyx = (gxy1 - gxy0)*f1o2d
  gxzx = (gxz1 - gxz0)*f1o2d
  gyyx = (gyy1 - gyy0)*f1o2d
  gyzx = (gyz1 - gyz0)*f1o2d
  gzzx = (gzz1 - gzz0)*f1o2d

  ! y-derivatives
  call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,x,y+delta,z,alp0,betax0,betay0,betaz0, &
                        gxx1,gxy1,gxz1,gyy1,gyz1,gzz1,gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
  beta_x1 = gxx1*betax0 + gxy1*betay0 + gxz1*betaz0
  beta_y1 = gxy1*betax0 + gyy1*betay0 + gyz1*betaz0
  beta_z1 = gxz1*betax0 + gyz1*betay0 + gzz1*betaz0

  call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,x,y-delta,z,alp0,betax0,betay0,betaz0, &
                        gxx0,gxy0,gxz0,gyy0,gyz0,gzz0,gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
  beta_x0 = gxx0*betax0 + gxy0*betay0 + gxz0*betaz0
  beta_y0 = gxy0*betax0 + gyy0*betay0 + gyz0*betaz0
  beta_z0 = gxz0*betax0 + gyz0*betay0 + gzz0*betaz0

  ! partial_y beta_i and partial_y g_ij
  beta_xyt = (beta_x1 - beta_x0)*f1o2d
  beta_yyt = (beta_y1 - beta_y0)*f1o2d
  beta_zyt = (beta_z1 - beta_z0)*f1o2d
  gxxy = (gxx1 - gxx0)*f1o2d
  gxyy = (gxy1 - gxy0)*f1o2d
  gxzy = (gxz1 - gxz0)*f1o2d
  gyyy = (gyy1 - gyy0)*f1o2d
  gyzy = (gyz1 - gyz0)*f1o2d
  gzzy = (gzz1 - gzz0)*f1o2d

  ! z-derivatives
  call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,x,y,z+delta,alp0,betax0,betay0,betaz0, &
                        gxx1,gxy1,gxz1,gyy1,gyz1,gzz1,gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
  beta_x1 = gxx1*betax0 + gxy1*betay0 + gxz1*betaz0
  beta_y1 = gxy1*betax0 + gyy1*betay0 + gyz1*betaz0
  beta_z1 = gxz1*betax0 + gyz1*betay0 + gzz1*betaz0

  call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,x,y,z-delta,alp0,betax0,betay0,betaz0, &
                        gxx0,gxy0,gxz0,gyy0,gyz0,gzz0,gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
  beta_x0 = gxx0*betax0 + gxy0*betay0 + gxz0*betaz0
  beta_y0 = gxy0*betax0 + gyy0*betay0 + gyz0*betaz0
  beta_z0 = gxz0*betax0 + gyz0*betay0 + gzz0*betaz0

  ! partial_z beta_i and partial_z g_ij
  beta_xzt = (beta_x1 - beta_x0)*f1o2d
  beta_yzt = (beta_y1 - beta_y0)*f1o2d
  beta_zzt = (beta_z1 - beta_z0)*f1o2d
  gxxz = (gxx1 - gxx0)*f1o2d
  gxyz = (gxy1 - gxy0)*f1o2d
  gxzz = (gxz1 - gxz0)*f1o2d
  gyyz = (gyy1 - gyy0)*f1o2d
  gyzz = (gyz1 - gyz0)*f1o2d
  gzzz = (gzz1 - gzz0)*f1o2d

  ! compute Gamma_{ijk}
  gamxxx = 0.5d0*gxxx
  gamxxy = 0.5d0*gxxy
  gamxxz = 0.5d0*gxxz
  gamxyy = gxyy - 0.5d0*gyyx
  gamxyz = 0.5d0*(gxyz + gxzy - gyzx)
  gamxzz = gxzz - 0.5d0*gzzx
  gamyxx = gxyx - 0.5d0*gxxy
  gamyxy = 0.5d0*gyyx
  gamyxz = 0.5d0*(gxyz + gyzx - gxzy)
  gamyyy = 0.5d0*gyyy
  gamyyz = 0.5d0*gyyz
  gamyzz = gyzz - 0.5d0*gzzy
  gamzxx = gxzx - 0.5d0*gxxz
  gamzxy = 0.5d0*(gxzy + gyzx - gxyz)
  gamzxz = 0.5d0*gzzx
  gamzyy = gyzy - 0.5d0*gyyz
  gamzyz = 0.5d0*gzzy
  gamzzz = 0.5d0*gzzz

  ! Compute D_i beta_j + D_j beta_i
  beta_xx = 2.d0*beta_xxt - 2.d0*(gamxxx*betax + gamyxx*betay + gamzxx*betaz)
  beta_xy = beta_xyt + beta_yxt - 2.d0*(gamxxy*betax + gamyxy*betay + gamzxy*betaz)
  beta_xz = beta_xzt + beta_zxt - 2.d0*(gamxxz*betax + gamyxz*betay + gamzxz*betaz)
  beta_yy = 2.d0*beta_yyt - 2.d0*(gamxyy*betax + gamyyy*betay + gamzyy*betaz)
  beta_yz = beta_yzt + beta_zyt - 2.d0*(gamxyz*betax + gamyyz*betay + gamzyz*betaz)
  beta_zz = 2.d0*beta_zzt - 2.d0*(gamxzz*betax + gamyzz*betay + gamzzz*betaz)

end subroutine cov_derv_beta

subroutine boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,x,y,z,alp,betax,betay,betaz, &
                        gxx,gxy,gxz,gyy,gyz,gzz,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
  implicit none
  real*8,  intent(in)                                :: m0,a,x,y,z,vx_boost
  real*8,  intent(out)                               :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  intent(out)                               :: gupxx,gupxy,gupxz,gupyy
  real*8,  intent(out)                               :: gupyz,gupzz,alp,betax,betay,betaz
  real*8 :: gtt0,gtx0,gty0,gtz0,gxx0,gxy0,gxz0,gyy0,gyz0,gzz0
  real*8 :: g4uptt0,g4uptx0,g4upty0,g4uptz0,gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0
  real*8 :: betax0,betay0,betaz0,g4upxx0,g4upxy0,g4upxz0,g4upyy0,g4upyz0,g4upzz0
  real*8 :: g4uptt,g4uptx,g4upty,g4uptz
  real*8 :: gamma,gammav,Ltt,Ltx,Lxx,Lxt
!
  gamma = 1.d0/sqrt(1.d0-vx_boost*vx_boost)
  gammav = gamma*vx_boost
  ! \partial x^\mu / \partial x_boost^\nu
  Ltx = -gammav
  Lxx = gamma

  call compute_kerr_puncture_metric(m0,a,gamma*x,y,z,gtt0,gtx0,gty0,gtz0, &
                        gxx0,gxy0,gxz0,gyy0,gyz0,gzz0,g4uptt0,g4uptx0,g4upty0,g4uptz0, &
                        gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)

  betax0 = -g4uptx0/g4uptt0
  betay0 = -g4upty0/g4uptt0
  betaz0 = -g4uptz0/g4uptt0
  g4upxx0 = gupxx0 + betax0*betax0*g4uptt0
  g4upxy0 = gupxy0 + betax0*betay0*g4uptt0
  g4upxz0 = gupxz0 + betax0*betaz0*g4uptt0
  g4upyy0 = gupyy0 + betay0*betay0*g4uptt0
  g4upyz0 = gupyz0 + betay0*betaz0*g4uptt0
  g4upzz0 = gupzz0 + betaz0*betaz0*g4uptt0

  gxx = Ltx*Ltx*gtt0 + 2.d0*Ltx*Lxx*gtx0 + Lxx*Lxx*gxx0
  gxy = Ltx*gty0 + Lxx*gxy0
  gxz = Ltx*gtz0 + Lxx*gxz0
  gyy = gyy0
  gyz = gyz0
  gzz = gzz0

  ! \partial x_boost^\mu / \partial x^nu
  Ltt = gamma
  Ltx = gammav
  Lxt = gammav
  Lxx = gamma

  g4uptt = Ltt*Ltt*g4uptt0 + 2.d0*Ltt*Ltx*g4uptx0 + Ltx*Ltx*g4upxx0
  g4uptx = Ltt*Lxt*g4uptt0 + (Ltt*Lxx+Ltx*Lxt)*g4uptx0 + Ltx*Lxx*g4upxx0
  g4upty = Ltt*g4upty0 + Ltx*g4upxy0
  g4uptz = Ltt*g4uptz0 + Ltx*g4upxz0

  alp = 1.d0/sqrt(-g4uptt)
  betax = -g4uptx/g4uptt
  betay = -g4upty/g4uptt
  betaz = -g4uptz/g4uptt

  gupxx = Lxt*Lxt*g4uptt0 + 2.d0*Lxt*Lxx*g4uptx0 + Lxx*Lxx*g4upxx0 - betax*betax*g4uptt
  gupxy = Lxt*g4upty0 + Lxx*g4upxy0 - betax*betay*g4uptt
  gupxz = Lxt*g4uptz0 + Lxx*g4upxz0 - betax*betaz*g4uptt
  gupyy = g4upyy0 - betay*betay*g4uptt
  gupyz = g4upyz0 - betay*betaz*g4uptt
  gupzz = g4upzz0 - betaz*betaz*g4uptt

end subroutine boost_kerr_puncture_alp_beta_gij_gupij

subroutine compute_kerr_puncture_metric(m0,a,x,y,z,gtt,gtx,gty,gtz, & 
                        gxx,gxy,gxz,gyy,gyz,gzz,g4uptt,g4uptx,g4upty,g4uptz, & 
                        gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
  implicit none
  real*8,  intent(in)                                :: m0,a,x,y,z
  real*8,  intent(out)                               :: gtt,gtx,gty,gtz,gxx,gyy,gzz
  real*8,  intent(out)                               :: gxy,gxz,gyz
  real*8,  intent(out)                               :: g4uptt,g4uptx,g4upty,g4uptz
  real*8,  intent(out)                               :: gupxx,gupxy,gupxz,gupyy
  real*8,  intent(out)                               :: gupyz,gupzz
  real*8 :: grr, gthth, gphph, alp,betaph,betax,betay,betaz
  real*8 :: beta_x, beta_y, beta_z, beta2
  real*8 :: guprr, gupthth, gupphph, r2
  real*8 :: Sigma, aa, r, sin2th, sinth, cos2th, costh, lorentz_gamma
  real*8 :: Lx_r, Ly_r, Lz_r, Lx_th, Ly_th, Lz_th, Lx_ph, Ly_ph, Lz_ph
  real*8 :: Lr_x, Lth_x, Lph_x, Lr_y, Lth_y, Lph_y, Lr_z, Lth_z, Lph_z
  real*8 :: r_l, x_l, y_l, z_l, Delta,pomega
  real*8  :: rh_out, rh_in, riso, dr_driso
!
  rh_out = m0 + sqrt(m0*m0 - a*a)
  rh_in  = m0 - sqrt(m0*m0 - a*a)
 
  riso = sqrt(x*x+y*y+z*z)
  r = riso*(1.d0 + 0.25d0*rh_out/riso)**2
  dr_driso = 1.d0 - (0.25d0*rh_out/riso)**2
  costh = z/riso
  cos2th = costh**2
  sin2th = 1.d0 - cos2th
  sinth = sqrt(sin2th)
  r2 = r * r
  Sigma = r2 + a*a*cos2th
  Delta = r2 - 2.d0*m0*r + a*a
  aa = (r2+a*a)**2 - a*a*Delta*sin2th
!
!-- the lapse and shift (nonzero component)
!
  alp = sqrt(Delta*Sigma/aa)
  betaph = -2.d0*m0*a*r/aa

!
!-- the 3-metric (nonzero components)
!
  if (riso .ge. 0.25d0*rh_out) then
     grr = Sigma/64.d0/(r-rh_in)*( (sqrt(r-rh_out)+sqrt(r))*(4.d0*riso+rh_out)/riso/riso )**2
  else
      grr = Sigma/64.d0/(r-rh_in)*( rh_out/(sqrt(r-rh_out)+sqrt(r))*(4.d0*riso+rh_out)/riso/riso )**2
  end if
  gthth = Sigma
  gphph = aa*sin2th/Sigma
  guprr = 1.d0/grr
  gupthth = 1.d0/gthth
  gupphph = 1.d0/gphph

  ! Transform to Cartesian components
  x_l = x
  y_l = y
  z_l = z
  r_l = riso
  r2 = riso*riso
  pomega = sqrt(x*x+y*y)
  Lx_r = x_l/r_l  ! = sin(th)cos(ph)
  Ly_r = y_l/r_l  ! = sin(th)sin(ph)
  Lz_r = z_l/r_l  ! = cos(ph)

  !!Lx_th = z_l*x_l/sqrt(r2 - z_l**2)
  !!Ly_th = z_l*y_l/sqrt(r2 - z_l**2)
  !!Lz_th = -sqrt(r2 - z_l**2)
  Lx_th = z_l*x_l/pomega
  Ly_th = z_l*y_l/pomega
  Lz_th = -pomega

  Lx_ph = -y_l
  Ly_ph = x_l
  Lz_ph = 0.d0

  Lr_x = Lx_r
  Lr_y = Ly_r
  Lr_z = Lz_r

  !!Lth_x = z_l*x_l/sqrt(r2-z_l**2)/(r_l*r_l)
  !!Lth_y = z_l*y_l/sqrt(r2-z_l**2)/(r_l*r_l)
  !!Lth_z = -sqrt(r2-z_l**2)/(r_l*r_l)
  Lth_x = z_l*x_l/(pomega*r2) 
  Lth_y = z_l*y_l/(pomega*r2)
  Lth_z = -pomega/r2

  Lph_x = -y_l/(pomega*pomega)
  Lph_y = x_l/(pomega*pomega)
  Lph_z = 0.d0

  betax = Lx_ph*betaph
  betay = Ly_ph*betaph
  betaz = 0.d0

  ! Transform the spatial metric
  gxx = Lr_x**2 * grr + &
               Lth_x**2 * gthth + &
               Lph_x**2 * gphph

  gxy = Lr_x*Lr_y * grr + &
               Lth_x*Lth_y * gthth + &
               Lph_x*Lph_y * gphph

  gxz = Lr_x*Lr_z * grr + &
               Lth_x*Lth_z * gthth + &
               Lph_x*Lph_z * gphph

  gyy = Lr_y**2 * grr + &
               Lth_y**2 * gthth + &
               Lph_y**2 * gphph

  gyz = Lr_y*Lr_z * grr + &
               Lth_y*Lth_z * gthth + &
               Lph_y*Lph_z * gphph

  gzz = Lr_z**2 * grr + &
               Lth_z**2 * gthth + &
               Lph_z**2 * gphph

  gupxx = Lx_r**2 * guprr + &
               Lx_th**2 * gupthth + &
               Lx_ph**2 * gupphph

  gupxy = Lx_r*Ly_r * guprr + &
               Lx_th*Ly_th * gupthth + &
               Lx_ph*Ly_ph * gupphph

  gupxz = Lx_r*Lz_r * guprr + &
               Lx_th*Lz_th * gupthth + &
               Lx_ph*Lz_ph * gupphph

  gupyy = Ly_r**2 * guprr + &
               Ly_th**2 * gupthth + &
               Ly_ph**2 * gupphph

  gupyz = Ly_r*Lz_r * guprr + &
               Ly_th*Lz_th * gupthth + &
               Ly_ph*Lz_ph * gupphph

  gupzz = Lz_r**2 * guprr + &
               Lz_th**2 * gupthth + &
               Lz_ph**2 * gupphph

  beta_x = gxx*betax + gxy*betay + gxz*betaz
  beta_y = gxy*betax + gyy*betay + gyz*betaz
  beta_z = gxz*betax + gyz*betay + gzz*betaz

  beta2 = betax*beta_x + betay*beta_y + betaz*beta_z

  gtt = beta2 - alp*alp
  gtx = beta_x
  gty = beta_y
  gtz = beta_z

  g4uptt = -1.d0/(alp*alp)
  g4uptx = -betax*g4uptt
  g4upty = -betay*g4uptt
  g4uptz = -betaz*g4uptt

end subroutine compute_kerr_puncture_metric
