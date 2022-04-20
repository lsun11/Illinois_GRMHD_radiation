!-----------------------------------------------------------------------------
!
! $Id: compute_bwh.f90
!
!-----------------------------------------------------------------------------
!
! Compute matter source terms
!
!-----------------------------------------------------------------------------
! Corotational case
!-----------------------------------------------------------------------------
subroutine compute_cbwh(ex, Omega, Omega_Frame, q_atmos, n, &
     X, Y, Z, &
     qi, lapse, phi, shiftx, shifty, shiftz, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     S, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     rho, rho_star, tau, st_x, st_y, st_z, &
     P, w, vx, vy, vz, &
     rho_b, u0, h);
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8                                   :: Omega, Omega_Frame, q_atmos, n
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: qi,lapse,phi
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx,shifty,shiftz
  real*8, dimension(ex(1),ex(2),ex(3))     :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3))     :: S,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3))     :: Sxx,Sxy,Sxz,Syy,Syz,Szz
  real*8, dimension(ex(1),ex(2),ex(3))     :: rho, rho_star,tau,st_x,st_y,st_z
  real*8, dimension(ex(1),ex(2),ex(3))     :: P,w,vx,vy,vz,rho_b,u0,h
!
! Other variables:
!
  real*8                     :: rho0,rhoi,Press,ut,alpha,Gam,U
  real*8                     :: fac2,psi,psi4,psi6,psi8
  real*8                     :: dOmega,r
  real*8                     :: pgxx,pgxy,pgxz,pgyy,pgyz,pgzz
  real*8                     :: q,v2,gamma2,enth
  real*8                     :: ux,uy,uz,u_x,u_y,u_z,betax,betay
  integer                    :: i, j, k, index
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8                     :: ZERO, ONE, THREE, Tiny
  parameter(ZERO = 0.0D0, ONE = 1.D0, THREE = 3.D0, Tiny = 1.D-10)
!
! compute
!
  imin = lbound(phi,1)
  jmin = lbound(phi,2)
  kmin = lbound(phi,3)
  imax = ubound(phi,1)
  jmax = ubound(phi,2)
  kmax = ubound(phi,3)
  Gam   = ONE + ONE/n
  dOmega = Omega - Omega_Frame
  do i = imin, imax
     do j = jmin, jmax
        do k = kmin, kmax
!  NOTE: There is no atmosphere (yet)
           if (qi(i,j,k).lt.Tiny) then
                   shiftx(i,j,k) = shiftx(i,j,k) - Omega_Frame*Y(1,j,1)
                   shifty(i,j,k) = shifty(i,j,k) + Omega_Frame*X(i,1,1)
                   h(i,j,k)         = ONE
                   rho(i,j,k)       = ZERO
                   S(i,j,k)         = ZERO
                   Sx(i,j,k)        = ZERO
                   Sy(i,j,k)        = ZERO
                   Sz(i,j,k)        = ZERO
                   Sxx(i,j,k)       = ZERO
                   Sxy(i,j,k)       = ZERO
                   Sxz(i,j,k)       = ZERO
                   Syy(i,j,k)       = ZERO
                   Syz(i,j,k)       = ZERO
                   Szz(i,j,k)       = ZERO
                   vx(i,j,k)        = -shiftx(i,j,k)
                   vy(i,j,k)        = -shifty(i,j,k)
                   vz(i,j,k)        = ZERO
                   rho_star(i,j,k)  = ZERO
                   st_x(i,j,k)      = ZERO
                   st_y(i,j,k)      = ZERO
                   st_z(i,j,k)      = ZERO
                   tau(i,j,k)       = ZERO
                   P(i,j,k)         = ZERO
                   w(i,j,k)         = ZERO
                   u0(i,j,k)        = ONE/(lapse(i,j,k) + ONE)
                   rho_b(i,j,k)     = ZERO
           else
                q     =  qi(i,j,k)
                   alpha = lapse(i,j,k) + ONE
                   psi  = exp(phi(i,j,k))
                   psi4 = psi**4
                   psi6 = psi4*psi*psi
                   psi8 = psi4*psi4
!                  NOTE: Since u^0 (ut) is the same in both frames 
!                (inertial and rotating) we calculate it
!                in the inertial frame
                   v2 = (psi4/(alpha**2))* &
                        ( (-Omega*Y(1,j,1)+shiftx(i,j,k) )**2 &
                        + ( Omega*X(i,1,1)+shifty(i,j,k) )**2 &
                        +              shiftz(i,j,k)  **2 )
                   gamma2 = ONE/(ONE-v2)
                   ut = sqrt(gamma2)/alpha
!                  Now we add the (Omega x r) part to the Shift vector
!                  and calculate u^i
                   shiftx(i,j,k) = shiftx(i,j,k) - Omega_Frame*Y(1,j,1)
                   shifty(i,j,k) = shifty(i,j,k) + Omega_Frame*X(i,1,1)
                   ux    = -ut * dOmega * Y(1,j,1)
                   uy    =  ut * dOmega * X(i,1,1)
                   U     = alpha*ut
                   pgxx  = psi4*gxx(i,j,k)
                   pgxy  = psi4*gxy(i,j,k)
                   pgxz  = psi4*gxz(i,j,k)
                   pgyy  = psi4*gyy(i,j,k)
                   pgyz  = psi4*gyz(i,j,k)
                   pgzz  = psi4*gzz(i,j,k)
                   betax = shiftx(i,j,k)
                   betay = shifty(i,j,k)
                   u_x   = pgxx*(ux + ut*betax) + pgxy*(uy + ut*betay)
                   u_y   = pgxy*(ux + ut*betax) + pgyy*(uy + ut*betay)
                   u_z   = pgxz*(ux + ut*betax) + pgyz*(uy + ut*betay)
                   rho0 = q**n
                   rhoi = n * q**(n+ONE)
                   Press = q**(n+ONE)
                   enth       = (rho0 + rhoi + Press)/rho0
                   fac2       = rho0*enth*U*U
                   h(i,j,k)         = enth
                   rho(i,j,k)       = fac2 - Press
                   S(i,j,k)         = THREE*Press + enth*rho0*(U*U - ONE)
                   Sx(i,j,k)        = rho0*enth*U*u_x
                   Sy(i,j,k)        = rho0*enth*U*u_y
                   Sz(i,j,k)        = rho0*enth*U*u_z
                   Sxx(i,j,k)       = Press*pgxx + Sx(i,j,k)*Sx(i,j,k)/fac2
                   Sxy(i,j,k)       = Press*pgxy + Sx(i,j,k)*Sy(i,j,k)/fac2
                   Sxz(i,j,k)       = Press*pgxz + Sx(i,j,k)*Sz(i,j,k)/fac2
                   Syy(i,j,k)       = Press*pgyy + Sy(i,j,k)*Sy(i,j,k)/fac2
                   Syz(i,j,k)       = Press*pgyz + Sy(i,j,k)*Sz(i,j,k)/fac2
                   Szz(i,j,k)       = Press*pgzz + Sz(i,j,k)*Sz(i,j,k)/fac2
                   vx(i,j,k)        = ux/ut
                   vy(i,j,k)        = uy/ut
                   vz(i,j,k)        = ZERO
                   rho_star(i,j,k)  = ut*alpha*psi6*rho0
                   st_x(i,j,k)      = rho_star(i,j,k)*enth*u_x
                   st_y(i,j,k)      = rho_star(i,j,k)*enth*u_y
                   st_z(i,j,k)      = rho_star(i,j,k)*enth*u_z
                   tau(i,j,k)       = (U*enth - ONE)*rho_star(i,j,k) - psi6*Press
                   P(i,j,k)         = Press
                   w(i,j,k)         = rho_star(i,j,k)*U
                   u0(i,j,k)        = ut
                   rho_b(i,j,k)     = rho0
           endif
         end do
     end do
  end do
end subroutine compute_cbwh
!-----------------------------------------------------------------------------
! Irrotational case
!-----------------------------------------------------------------------------
subroutine compute_ibwh(ex, Omega_Frame, rho_b_atm, n, &
     X, Y, Z, &
     qi, lapse, phi, shiftx, shifty, shiftz, &
     ux, uy, uz, ut, &
     rho, S, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     rho_star, tau, st_x, st_y, st_z, &
     P, w, vx, vy, vz, rho_b, u0, h);

  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8                                   :: Omega_Frame, rho_b_atm, n
  real*8, dimension(ex(1),ex(2),ex(3))     :: qi,lapse,phi,X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx,shifty,shiftz
  real*8, dimension(ex(1),ex(2),ex(3))     :: ux,uy,uz,ut
  real*8, dimension(ex(1),ex(2),ex(3))     :: rho,S,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3))     :: Sxx,Sxy,Sxz,Syy,Syz,Szz
  real*8, dimension(ex(1),ex(2),ex(3))     :: rho_star,tau,st_x,st_y,st_z
  real*8, dimension(ex(1),ex(2),ex(3))     :: P,w,vx,vy,vz,rho_b,u0,h
!
! Other variables:
!
  real*8                     :: rho0,taui,Press,enth,enth1,Gam,fac,alpha
  real*8                     :: psi,psi4,psi6,psi8,rhoi,fac2
  real*8                     :: r,Tiny,q,u_u,u_w,U,uux,uuy,uuz
  real*8                     :: gijuiuj,au0m1
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8, dimension(3)       :: ui
  real*8                     :: ZERO, HALF, ONE, TWO, THREE
  parameter(ZERO=0.D0,HALF=0.5D0,ONE=1.D0,TWO=2.D0,THREE=3.D0)

!
! compute
!
  imin = lbound(rho,1)
  jmin = lbound(rho,2)
  kmin = lbound(rho,3)
  imax = ubound(rho,1)
  jmax = ubound(rho,2)
  kmax = ubound(rho,3)

  Gam = ONE + ONE/n
  TINY = rho_b_atm**(1.d0/n)

  do i = imin, imax 
     do j = jmin, jmax
        do k = kmin, kmax

!  NOTE: There is no atmosphere (yet)
	
	   if (qi(i,j,k).lt.Tiny) then		! Outside the stars

           	shiftx(i,j,k) = shiftx(i,j,k) - Omega_Frame*Y(1,j,1)
           	shifty(i,j,k) = shifty(i,j,k) + Omega_Frame*X(i,1,1)

                ut(i,j,k) 	 = ONE/(lapse(i,j,k) + ONE)
           	h(i,j,k)         = ONE
           	rho(i,j,k)       = ZERO
           	S(i,j,k)         = ZERO
           	Sx(i,j,k)        = ZERO
           	Sy(i,j,k)        = ZERO
           	Sz(i,j,k)        = ZERO
           	Sxx(i,j,k)       = ZERO
           	Sxy(i,j,k)       = ZERO
           	Sxz(i,j,k)       = ZERO
           	Syy(i,j,k)       = ZERO
           	Syz(i,j,k)       = ZERO
           	Szz(i,j,k)       = ZERO
           	vx(i,j,k)       = -shiftx(i,j,k)
           	vy(i,j,k)       = -shifty(i,j,k)
           	vz(i,j,k)       = ZERO
           	rho_star(i,j,k)  = ZERO
           	st_x(i,j,k)      = ZERO
           	st_y(i,j,k)      = ZERO
           	st_z(i,j,k)      = ZERO
           	tau(i,j,k)       = ZERO
           	P(i,j,k)         = ZERO
           	w(i,j,k)         = ZERO
           	u0(i,j,k)        = ut(i,j,k)
           	rho_b(i,j,k)     = ZERO
	   else					! Inside the stars 

		q     =  qi(i,j,k)

           	alpha = lapse(i,j,k) + ONE

           	psi   = exp(phi(i,j,k))
           	psi4  = psi**4
           	psi6  = psi4*psi*psi
           	psi8  = psi4*psi4

                shiftx(i,j,k)  = shiftx(i,j,k) - Omega_Frame*Y(1,j,1)
                shifty(i,j,k)  = shifty(i,j,k) + Omega_Frame*X(i,1,1)

!	IMPORTANT: Note that 'ui' is u_i (lower index) and thus, it is
!		   the same in both the Inertial and Rotating Frame

                ui(1) = ux(i,j,k)
                ui(2) = uy(i,j,k)
                ui(3) = uz(i,j,k)

                u_u   = ui(1)*ui(1) + ui(2)*ui(2) + ui(3)*ui(3)
!                u_w   = ui(1)*shiftx(i,j,k) + ui(2)*shifty(i,j,k) + ui(3)*shiftz(i,j,k)
!		fac   = alpha*sqrt(ONE + u_u/psi4/alpha**2)
!
!                u0(i,j,k)	= (ONE/(alpha**2)) * (TWO*u_w + fac)
!                ut(i,j,k) 	= u0(i,j,k)
!
!                U     = alpha*u0(i,j,k)

		gijuiuj = u_u/psi4
	        au0m1 = gijuiuj/(ONE + sqrt(ONE + gijuiuj))
		U = au0m1 + ONE
	        u0(i,j,k) = U/alpha
	        ut(i,j,k) = u0(i,j,k)

!           	uux    = (ONE/(alpha**2)) * (shiftx(i,j,k)*fac+ui(1)/psi4)
!           	uuy    = (ONE/(alpha**2)) * (shifty(i,j,k)*fac+ui(2)/psi4)
	       
		uux = ui(1)/psi4 - shiftx(i,j,k)*u0(i,j,k)
		uuy = ui(2)/psi4 - shifty(i,j,k)*u0(i,j,k)
  		uuz = ui(3)/psi4 - shiftz(i,j,k)*u0(i,j,k)

           	rho0  = q**n
           	rhoi  = n * q**(n+ONE)
           	Press = q**(n+ONE)

  		enth1 = (rhoi + Press)/rho0
           	enth  = ONE + enth1
           	fac2  = rho0*enth*U*U

           	h(i,j,k)	= enth
           	rho(i,j,k)	= fac2 - Press

           	S(i,j,k)	= THREE*Press + enth*rho0*(U*U - ONE)
           	Sx(i,j,k)	= rho0*enth*U*ui(1)
           	Sy(i,j,k)	= rho0*enth*U*ui(2)
           	Sz(i,j,k)	= rho0*enth*U*ui(3)

           	Sxx(i,j,k)	= Press*psi4 + Sx(i,j,k)*Sx(i,j,k)/fac2
           	Sxy(i,j,k)	=              Sx(i,j,k)*Sy(i,j,k)/fac2
           	Sxz(i,j,k)	=              Sx(i,j,k)*Sz(i,j,k)/fac2
           	Syy(i,j,k)	= Press*psi4 + Sy(i,j,k)*Sy(i,j,k)/fac2
           	Syz(i,j,k)	=              Sy(i,j,k)*Sz(i,j,k)/fac2
           	Szz(i,j,k)	= Press*psi4 + Sz(i,j,k)*Sz(i,j,k)/fac2

           	vx(i,j,k)	= uux/u0(i,j,k)
           	vy(i,j,k)	= uuy/u0(i,j,k)
           	vz(i,j,k)	= uuz/u0(i,j,k)

           	rho_star(i,j,k)  = u0(i,j,k)*alpha*psi6*rho0
           	st_x(i,j,k)      = rho_star(i,j,k)*enth*ui(1)
           	st_y(i,j,k)      = rho_star(i,j,k)*enth*ui(2)
           	st_z(i,j,k)      = rho_star(i,j,k)*enth*ui(3)
!           	eps_star         = (rhoi/rho0)*(alpha*u0(i,j,k)*psi6)**(Gam-ONE)
!           	e_star(i,j,k)    = (rho_star(i,j,k)*eps_star)**(ONE/Gam)
                tau(i,j,k)       = (au0m1*enth + enth1)*rho_star(i,j,k) - psi6*Press
           	P(i,j,k)         = Press
           	w(i,j,k)         = rho_star(i,j,k)*U
           	rho_b(i,j,k)     = rho0
           end if
!           write(*,*) "hello.",i,j,k,vx(i,j,k),uux,u0(i,j,k),ui(1),psi4,shiftx(i,j,k)
!	if ((i.gt.21).and.(i.lt.23).and.		&
!	    (j.gt.29).and.(j.lt.31).and.		&
!	    (k.eq.1)) then
!		write(*,*) 'compute_cbwh i,j,k: ',i,j,k,X(i,1,1),Y(1,j,1)
!		write(*,*) 'compute_cbwh rho_b(i,j,k): ',rho_b(i,j,k),rho0
!		write(*,*) 'compute_cbwh rho_star(i,j,k): ',rho_star(i,j,k)
!		write(*,*) ' '
!	endif

         end do
     end do
  end do

end subroutine compute_ibwh
