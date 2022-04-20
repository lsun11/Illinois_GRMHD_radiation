! Compute \tilde{Gamma}^i by numerically differentiating the analytic values of 
! \tilde{gamma}^{ij}
!
subroutine setgamma_kerr_puncture(m0,a,x0_bh,vx_boost,ex,x,y,z,Gammax,Gammay,Gammaz)
  implicit none
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z
  real*8,                    intent(in) :: m0,a,vx_boost,x0_bh
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Gammax,Gammay,Gammaz
  integer :: i,j,k,imin,jmin,kmin,imax,jmax,kmax
  real*8 :: xl,yl,zl,rl,delta,psi4,f1o2d
  real*8 :: alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl
  real*8 :: gupxx1,gupxy1,gupxz1,gupyy1,gupyz1,gupzz1
  real*8 :: gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0
  real*8 :: gamx,gamy,gamz
!
  imin = lbound(Gammax,1)
  jmin = lbound(Gammax,2)
  kmin = lbound(Gammax,3)
  imax = ubound(Gammax,1)
  jmax = ubound(Gammax,2)
  kmax = ubound(Gammax,3)

  do k=kmin,kmax
     do j=jmin,jmax
        do i=imin,imax
           xl = x(i,j,k) - x0_bh
           yl = y(i,j,k)
           zl = z(i,j,k)
	   rl = sqrt(xl*xl + yl*yl + zl*zl)
	   delta = 1.d-4*(rl+1.d0)
	   f1o2d = 0.5d0/delta

           ! x-derivatives
           call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,xl+delta,yl,zl, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx1,gupxy1,gupxz1,gupyy1,gupyz1,gupzz1)
	   psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
	   gupxx1 = psi4*gupxx1
           gupxy1 = psi4*gupxy1
           gupxz1 = psi4*gupxz1
           !!gupyy1 = psi4*gupyy1
           !!gupyz1 = psi4*gupyz1
           !!gupzz1 = psi4*gupzz1

           call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,xl-delta,yl,zl, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
           psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
           gupxx0 = psi4*gupxx0
           gupxy0 = psi4*gupxy0
           gupxz0 = psi4*gupxz0
           !!gupyy0 = psi4*gupyy0
           !!gupyz0 = psi4*gupyz0
           !!gupzz0 = psi4*gupzz0

           gamx = (gupxx0-gupxx1)*f1o2d
           gamy = (gupxy0-gupxy1)*f1o2d
           gamz = (gupxz0-gupxz1)*f1o2d

           ! y-derivatives
           call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,xl,yl+delta,zl, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx1,gupxy1,gupxz1,gupyy1,gupyz1,gupzz1)
           psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
           !!gupxx1 = psi4*gupxx1
           gupxy1 = psi4*gupxy1
           !!gupxz1 = psi4*gupxz1
           gupyy1 = psi4*gupyy1
           gupyz1 = psi4*gupyz1
           !!gupzz1 = psi4*gupzz1

           call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,xl,yl-delta,zl, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
           psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
           !!gupxx0 = psi4*gupxx0
           gupxy0 = psi4*gupxy0
           !!gupxz0 = psi4*gupxz0
           gupyy0 = psi4*gupyy0
           gupyz0 = psi4*gupyz0
           !!gupzz0 = psi4*gupzz0

           gamx = gamx + (gupxy0-gupxy1)*f1o2d
           gamy = gamy + (gupyy0-gupyy1)*f1o2d
           gamz = gamz + (gupyz0-gupyz1)*f1o2d

           ! z-derivatives
           call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,xl,yl,zl+delta, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx1,gupxy1,gupxz1,gupyy1,gupyz1,gupzz1)
           psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
           !!gupxx1 = psi4*gupxx1
           !!gupxy1 = psi4*gupxy1
           gupxz1 = psi4*gupxz1
           !!gupyy1 = psi4*gupyy1
           gupyz1 = psi4*gupyz1
           gupzz1 = psi4*gupzz1

           call boost_kerr_puncture_alp_beta_gij_gupij(m0,a,vx_boost,xl,yl,zl-delta, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
           psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
           !!gupxx0 = psi4*gupxx0
           !!gupxy0 = psi4*gupxy0
           gupxz0 = psi4*gupxz0
           !!gupyy0 = psi4*gupyy0
           gupyz0 = psi4*gupyz0
           gupzz0 = psi4*gupzz0

           gamx = gamx + (gupxz0-gupxz1)*f1o2d
           gamy = gamy + (gupyz0-gupyz1)*f1o2d
           gamz = gamz + (gupzz0-gupzz1)*f1o2d

           Gammax(i,j,k) = gamx
           Gammay(i,j,k) = gamy
           Gammaz(i,j,k) = gamz

        end do
     end do
  end do
end subroutine setgamma_kerr_puncture

subroutine setgamma_kerr_quasiisotropic(m0,a,x0_bh,vx_boost,ex,x,y,z,Gammax,Gammay,Gammaz)
  implicit none
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z
  real*8,                    intent(in) :: m0,a,vx_boost,x0_bh
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Gammax,Gammay,Gammaz
  integer :: i,j,k,imin,jmin,kmin,imax,jmax,kmax
  real*8 :: xl,yl,zl,rl,delta,psi4,f1o2d
  real*8 :: alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl
  real*8 :: gupxx1,gupxy1,gupxz1,gupyy1,gupyz1,gupzz1
  real*8 :: gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0
  real*8 :: gamx,gamy,gamz
!
  imin = lbound(Gammax,1)
  jmin = lbound(Gammax,2)
  kmin = lbound(Gammax,3)
  imax = ubound(Gammax,1)
  jmax = ubound(Gammax,2)
  kmax = ubound(Gammax,3)

  do k=kmin,kmax
     do j=jmin,jmax
        do i=imin,imax
           xl = x(i,j,k) - x0_bh
           yl = y(i,j,k)
           zl = z(i,j,k)
	   rl = sqrt(xl*xl + yl*yl + zl*zl)
	   delta = 1.d-4*(rl+1.d0)
	   f1o2d = 0.5d0/delta

           ! x-derivatives
           call boost_kerr_quasiisotropic_alp_beta_gij_gupij(m0,a,vx_boost,xl+delta,yl,zl, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx1,gupxy1,gupxz1,gupyy1,gupyz1,gupzz1)
	   psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
	   gupxx1 = psi4*gupxx1
           gupxy1 = psi4*gupxy1
           gupxz1 = psi4*gupxz1
           !!gupyy1 = psi4*gupyy1
           !!gupyz1 = psi4*gupyz1
           !!gupzz1 = psi4*gupzz1

           call boost_kerr_quasiisotropic_alp_beta_gij_gupij(m0,a,vx_boost,xl-delta,yl,zl, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
           psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
           gupxx0 = psi4*gupxx0
           gupxy0 = psi4*gupxy0
           gupxz0 = psi4*gupxz0
           !!gupyy0 = psi4*gupyy0
           !!gupyz0 = psi4*gupyz0
           !!gupzz0 = psi4*gupzz0

           gamx = (gupxx0-gupxx1)*f1o2d
           gamy = (gupxy0-gupxy1)*f1o2d
           gamz = (gupxz0-gupxz1)*f1o2d

           ! y-derivatives
           call boost_kerr_quasiisotropic_alp_beta_gij_gupij(m0,a,vx_boost,xl,yl+delta,zl, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx1,gupxy1,gupxz1,gupyy1,gupyz1,gupzz1)
           psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
           !!gupxx1 = psi4*gupxx1
           gupxy1 = psi4*gupxy1
           !!gupxz1 = psi4*gupxz1
           gupyy1 = psi4*gupyy1
           gupyz1 = psi4*gupyz1
           !!gupzz1 = psi4*gupzz1

           call boost_kerr_quasiisotropic_alp_beta_gij_gupij(m0,a,vx_boost,xl,yl-delta,zl, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
           psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
           !!gupxx0 = psi4*gupxx0
           gupxy0 = psi4*gupxy0
           !!gupxz0 = psi4*gupxz0
           gupyy0 = psi4*gupyy0
           gupyz0 = psi4*gupyz0
           !!gupzz0 = psi4*gupzz0

           gamx = gamx + (gupxy0-gupxy1)*f1o2d
           gamy = gamy + (gupyy0-gupyy1)*f1o2d
           gamz = gamz + (gupyz0-gupyz1)*f1o2d

           ! z-derivatives
           call boost_kerr_quasiisotropic_alp_beta_gij_gupij(m0,a,vx_boost,xl,yl,zl+delta, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx1,gupxy1,gupxz1,gupyy1,gupyz1,gupzz1)
           psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
           !!gupxx1 = psi4*gupxx1
           !!gupxy1 = psi4*gupxy1
           gupxz1 = psi4*gupxz1
           !!gupyy1 = psi4*gupyy1
           gupyz1 = psi4*gupyz1
           gupzz1 = psi4*gupzz1

           call boost_kerr_quasiisotropic_alp_beta_gij_gupij(m0,a,vx_boost,xl,yl,zl-delta, &
                        alpl,betaxl,betayl,betazl,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                        gupxx0,gupxy0,gupxz0,gupyy0,gupyz0,gupzz0)
           psi4 = (gxxl * gyyl * gzzl + gxyl * gyzl * gxzl + &
                        gxzl * gxyl * gyzl - gxzl * gyyl * gxzl &
                      - gxyl * gxyl * gzzl - gxxl * gyzl * gyzl)**(1.d0/3.d0)
           !!gupxx0 = psi4*gupxx0
           !!gupxy0 = psi4*gupxy0
           gupxz0 = psi4*gupxz0
           !!gupyy0 = psi4*gupyy0
           gupyz0 = psi4*gupyz0
           gupzz0 = psi4*gupzz0

           gamx = gamx + (gupxz0-gupxz1)*f1o2d
           gamy = gamy + (gupyz0-gupyz1)*f1o2d
           gamz = gamz + (gupzz0-gupzz1)*f1o2d

           Gammax(i,j,k) = gamx
           Gammay(i,j,k) = gamy
           Gammaz(i,j,k) = gamz

        end do
     end do
  end do
end subroutine setgamma_kerr_quasiisotropic
