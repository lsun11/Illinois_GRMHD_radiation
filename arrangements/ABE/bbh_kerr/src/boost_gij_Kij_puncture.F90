subroutine compute_boost_gij_Kij(m0,sam,x0,y0,z0,ux_boost,ex,x,y,z, &
 			lapm1,betax,betay,betaz,gxx,gxy,gxz, & 
                        gyy,gyz,gzz,Kxx,Kxy,Kxz,Kyy,Kyz,Kzz)
  implicit none
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z
  real*8,                    intent(in) :: m0,sam,x0,y0,z0,ux_boost
  real*8,  dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)) :: Kxx,Kxy,Kxz,Kyy,Kyz,Kzz
  real*8,  dimension(ex(1),ex(2),ex(3)) :: lapm1,betax,betay,betaz
  integer :: i,j,k,imin,jmin,kmin,imax,jmax,kmax
  real*8 :: xl,yl,zl,f1o3,f1o2alp,riso2,lorentz_gamma,rh,rh2
  real*8 :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, vx_boost,a
  real*8 :: gupxxl,gupxyl,gupxzl,gupyyl,gupyzl,gupzzl
  real*8 :: alpl,betaxl,betayl,betazl
  real*8 :: Kxxl,Kxyl,Kxzl,Kyyl,Kyzl,Kzzl
  real*8 :: gxxt,gxyt,gxzt,gyyt,gyzt,gzzt
  real*8 :: beta_xx, beta_xy, beta_xz, beta_yy, beta_yz, beta_zz
!
  imin = lbound(gxx,1)
  jmin = lbound(gxx,2)
  kmin = lbound(gxx,3)
  imax = ubound(gxx,1)
  jmax = ubound(gxx,2)
  kmax = ubound(gxx,3)
 
  f1o3 = 1.d0/3.d0
  lorentz_gamma = sqrt(1.d0+ux_boost**2)
  vx_boost = ux_boost / lorentz_gamma
  a = m0*sam
  ! Horizon radius in "isotropic" radial coordinate
  rh = 0.25d0*( m0 + sqrt(m0*m0 - a*a) )
  rh2 = rh*rh

  do k=kmin,kmax
     do j=jmin,jmax
        do i=imin,imax
	   xl = x(i,j,k) - x0
	   yl = y(i,j,k) - y0
	   zl = z(i,j,k) - z0
           call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,xl,yl,zl, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, & 
                        gupxxl,gupxyl,gupxzl,gupyyl,gupyzl,gupzzl)
           lapm1(i,j,k) = lapm1(i,j,k) + (alpl - 1.d0)
	   betax(i,j,k) = betax(i,j,k) + betaxl
           betay(i,j,k) = betay(i,j,k) + betayl
           betaz(i,j,k) = betaz(i,j,k) + betazl
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

	   gxx(i,j,k) = gxx(i,j,k) + (gxxl - 1.d0)
           gxy(i,j,k) = gxy(i,j,k) + gxyl
           gxz(i,j,k) = gxz(i,j,k) + gxzl
           gyy(i,j,k) = gyy(i,j,k) + (gyyl - 1.d0)
           gyz(i,j,k) = gyz(i,j,k) + gyzl
           gzz(i,j,k) = gzz(i,j,k) + (gzzl - 1.d0)
           Kxx(i,j,k) = Kxx(i,j,k) + Kxxl
           Kxy(i,j,k) = Kxy(i,j,k) + Kxyl
           Kxz(i,j,k) = Kxz(i,j,k) + Kxzl
           Kyy(i,j,k) = Kyy(i,j,k) + Kyyl
           Kyz(i,j,k) = Kyz(i,j,k) + Kyzl
           Kzz(i,j,k) = Kzz(i,j,k) + Kzzl
	end do
     end do
  end do
end subroutine compute_boost_gij_Kij
