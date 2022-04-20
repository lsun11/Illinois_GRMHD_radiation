!-----------------------------------------------------------------------------
!
! $Id: compute_bns.f90,v 1.15 2000/11/22 17:43:17 mduez Exp $
!
!-----------------------------------------------------------------------------
!
! Compute matter source terms
!
!-----------------------------------------------------------------------------
! Corotational case
!-----------------------------------------------------------------------------
subroutine compute_bns(ex, re, C, Omega, ntot, Nx, Ny, X, Y, Z, qi, &
     lapse, phi, shiftx, shifty, shiftz, &
     rho, S, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz);
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8                                   :: re,C,Omega
  integer                                  :: ntot,Nx,Ny
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: qi,lapse,phi
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx,shifty,shiftz
  real*8, dimension(ex(1),ex(2),ex(3))     :: rho,S,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3))     :: Sxx,Sxy,Sxz,Syy,Syz,Szz
!
! Other variables:
!
  real*8                     :: rho0,rhoi,P,n,ut2,fac,alpha,stuff
  real*8                     :: psi,psi4,psi8
  real*8                     :: q,gamma,v2,gamma2,fac2
  integer                    :: i, j, k, index
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8, dimension(3)       :: w,xi
  real*8                     :: HALF, ONE, TWO, THREE, FOUR, EIGHT, PI
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0, FOUR = 4.D0, EIGHT = 8.D0)
  PI = acos(-ONE)
!
! compute
!
  imin = lbound(rho,1)
  jmin = lbound(rho,2)
  kmin = lbound(rho,3)
  imax = ubound(rho,1)
  jmax = ubound(rho,2)
  kmax = ubound(rho,3)
  n = ONE
  do i = imin, imax
     do j = jmin, jmax
        do k = kmin, kmax
           alpha = lapse(i,j,k) + ONE
           q = qi(i,j,k)
           psi = exp(phi(i,j,k))
           psi4 = psi*psi*psi*psi
           psi8 = psi4*psi4
           w(1) = -shiftx(i,j,k)
           w(2) = -shifty(i,j,k)
           w(3) = -shiftz(i,j,k)
           xi(1) = Y(1,j,1)
           xi(2) = -X(i,1,1)
           xi(3) = 0.D0
           v2 = (psi4/(alpha**2))*((Omega*Y(1,j,1)-w(1))**2 &
                + (w(3))**TWO + (Omega*X(i,1,1)+w(2))**2)
           gamma2 = ONE/(ONE-v2)
! Note:  gamma = 1/(1-v**2)**1/2
           rho0 = q**n
           rhoi = n * q**(n+ONE)
           P = q**(n+ONE)
           ut2 = gamma2/alpha**2
           rho(i,j,k) = (rho0 + rhoi + P)*(gamma2) - P
           S(i,j,k) = (rho0 + rhoi + P)*(gamma2)*v2 + THREE*P
           fac = q**n * ((ONE+(ONE+n)*q)/(alpha)) * gamma2 * psi4
           Sx(i,j,k) = fac * (Omega * xi(1) - w(1))
           Sy(i,j,k) = fac * (Omega * xi(2) - w(2))
           Sz(i,j,k) = fac * (Omega * xi(3) - w(3))
           fac = ut2 * psi8 * (rho0 + rhoi + P)
           fac2 = psi4 * P
           Sxx(i,j,k) = fac * ((Omega*Omega)*xi(1)*xi(1) - Omega*xi(1)*w(1) &
                - Omega*w(1)*xi(1) + w(1)*w(1)) + fac2
           Sxy(i,j,k) = fac * ((Omega*Omega)*xi(1)*xi(2) - Omega*xi(1)*w(2) &
                - Omega*w(1)*xi(2) + w(1)*w(2))
           Sxz(i,j,k) = fac * ((Omega*Omega)*xi(1)*xi(3) - Omega*xi(1)*w(3) &
                - Omega*w(1)*xi(3) + w(1)*w(3))
           Syy(i,j,k) = fac * ((Omega*Omega)*xi(2)*xi(2) - Omega*xi(2)*w(2) &
                - Omega*w(2)*xi(2) + w(2)*w(2)) + fac2
           Syz(i,j,k) = fac * ((Omega*Omega)*xi(2)*xi(3) - Omega*xi(2)*w(3) &
                - Omega*w(2)*xi(3) + w(2)*w(3))
           Szz(i,j,k) = fac * ((Omega*Omega)*xi(3)*xi(3) - Omega*xi(3)*w(3) &
                - Omega*w(3)*xi(3) + w(3)*w(3)) + fac2
         end do
     end do
  end do
end subroutine compute_bns
!-----------------------------------------------------------------------------
! Irrotational case
!-----------------------------------------------------------------------------
subroutine compute_irr(ex, X, Y, Z, qi, &
     lapse, phi, shiftx, shifty, shiftz, &
     ux, uy, uz, ut, &
     rho, S, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz);
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  integer                                  :: ntot,Nx,Ny
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: qi,lapse,phi
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx,shifty,shiftz
  real*8, dimension(ex(1),ex(2),ex(3))     :: ux,uy,uz,ut
  real*8, dimension(ex(1),ex(2),ex(3))     :: rho,S,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3))     :: Sxx,Sxy,Sxz,Syy,Syz,Szz
!
! Other variables:
!
  real*8                     :: rho0,eps,Press,enth,n,Gam,fac,alpha
  real*8                     :: psi,psi4,psi8
  real*8                     :: u0,w_w,u_u,u_w,U_s,u_x,u_y,u_z
  integer                    :: i, j, k, index
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8, dimension(3)       :: w,ui
  real*8                     :: HALF, ONE, TWO, ZERO, FOUR, EIGHT, PI
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0)
  parameter(FOUR = 4.D0, EIGHT = 8.D0)
  PI = acos(-ONE)
!
! compute
!
  imin = lbound(rho,1)
  jmin = lbound(rho,2)
  kmin = lbound(rho,3)
  imax = ubound(rho,1)
  jmax = ubound(rho,2)
  kmax = ubound(rho,3)
  n = ONE
  Gam = ONE + ONE/n
  do i = imin, imax
     do j = jmin, jmax
        do k = kmin, kmax
           alpha = lapse(i,j,k) + ONE
           rho0 = qi(i,j,k)
           if (rho0>ZERO) then
              psi = exp(phi(i,j,k))
              psi4 = psi*psi*psi*psi
              psi8 = psi4*psi4
              Press = rho0**Gam
              eps = (rho0**(ONE/n))/(Gam - ONE)
              enth = ONE + eps + Press/rho0
              w(1) = -shiftx(i,j,k)
              w(2) = -shifty(i,j,k)
              w(3) = -shiftz(i,j,k)
              ui(1) = ux(i,j,k)
              ui(2) = uy(i,j,k)
              ui(3) = uz(i,j,k)
              u_u = ui(1)*ui(1) + ui(2)*ui(2) + ui(3)*ui(3)
              u_w = ui(1)*w(1) + ui(2)*w(2) + ui(3)*w(3)
              u0 = (ONE/(alpha**2))*(-TWO*u_w + alpha*sqrt(ONE + u_u/psi4/alpha**2))
              U_s = alpha*u0
              fac = rho0*enth*U_s**2
              ut(i,j,k) = u0
              rho(i,j,k) = fac - Press
              S(i,j,k) = 3.D0*Press + rho0*enth*(U_s**2 - ONE)
              Sx(i,j,k) = rho0*enth*U_s*ui(1)
              Sy(i,j,k) = rho0*enth*U_s*ui(2)
              Sz(i,j,k) = rho0*enth*U_s*ui(3)
              Sxx(i,j,k) = Press*psi4 + Sx(i,j,k)*Sx(i,j,k)/(fac)
              Sxy(i,j,k) = Sx(i,j,k)*Sy(i,j,k)/(fac)
              Sxz(i,j,k) = Sx(i,j,k)*Sz(i,j,k)/(fac)
              Syy(i,j,k) = Press*psi4 + Sy(i,j,k)*Sy(i,j,k)/(fac)
              Syz(i,j,k) = Sy(i,j,k)*Sz(i,j,k)/(fac)
              Szz(i,j,k) = Press*psi4 + Sz(i,j,k)*Sz(i,j,k)/(fac)
           else
              ut(i,j,k) = ONE/alpha
              rho(i,j,k) = ZERO
              S(i,j,k) = ZERO
              Sx(i,j,k) = ZERO
              Sy(i,j,k) = ZERO
              Sz(i,j,k) = ZERO
              Sxx(i,j,k) = ZERO
              Sxy(i,j,k) = ZERO
              Sxz(i,j,k) = ZERO
              Syy(i,j,k) = ZERO
              Syz(i,j,k) = ZERO
              Szz(i,j,k) = ZERO
           end if
         end do
     end do
  end do
end subroutine compute_irr
!
!-----------------------------------------------------------------------------
!
! Integrate the rest mass integrand over the grid.  Since the functions are
! already centered, the integral is just a sum times dX dY dZ.
!
!-----------------------------------------------------------------------------
subroutine bns_rest_mass(ex, dmass, &
     X, Y, Z, &
     q, &
     phi, lapse, &
     adjimin, adjimax, &
     adjjmin, adjjmax, &
     adjkmin, adjkmax,C)
!  use tensor
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: q
  real*8, dimension(ex(1),ex(2),ex(3))        :: phi, lapse
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  real*8                                      :: C
!
! output:
! 
  real*8                                      :: dmass
!
! Other variables:
!
  real*8                             :: dX, dY, dZ
  real*8                             :: dV
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  real*8                             :: ZERO, ONE, SIX
  real*8                             :: u0, rho0,gamma11,n,alpha
  parameter( ZERO = 0.D0 )
  parameter( ONE = 1.D0 )
  parameter( SIX = 6.D0 )
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(phi,1) - adjimin
  jmin = lbound(phi,2) - adjjmin
  kmin = lbound(phi,3) - adjkmin
  imax = ubound(phi,1) - adjimax
  jmax = ubound(phi,2) - adjjmax
  kmax = ubound(phi,3) - adjkmax
!
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
  dV = dX * dY * dZ
!
! Set up integration
!
        write(*,*) ' We are in bns_rest_mass and C is:',C
  n = ONE
  dmass = ZERO
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           alpha = lapse(i,j,k) + ONE
           gamma11 = alpha * (q(i,j,k)*(n+ONE) + ONE)/(ONE + C)
           rho0 = q(i,j,k)**n
           u0 = gamma11/alpha
!
! Integration sum
!
           dmass = dmass + dV * u0 * &
                alpha * exp(SIX * phi(i,j,k)) * rho0
!!$           dmass = dmass + dV * u0 * &
!!$                alpha * rho0
        end do
     end do
  end do
  return
end subroutine bns_rest_mass
!--------------------------------------------------------
! Same as above, but used for irrotational stars.
!--------------------------------------------------------
subroutine irr_rest_mass(ex, dmass, &
     X, Y, Z, &
     rho_b, &
     phi, alpha, u0, &
     adjimin, adjimax, &
     adjjmin, adjjmax, &
     adjkmin, adjkmax)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_b
  real*8, dimension(ex(1),ex(2),ex(3))        :: phi, alpha, u0
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
!
! output:
! 
  real*8                                      :: dmass
!
! Other variables:
!
  real*8                             :: dX, dY, dZ
  real*8                             :: dV
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  real*8                             :: ZERO, ONE, SIX
  parameter( ZERO = 0.D0 )
  parameter( ONE = 1.D0 )
  parameter( SIX = 6.D0 )
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(phi,1) - adjimin
  jmin = lbound(phi,2) - adjjmin
  kmin = lbound(phi,3) - adjkmin
  imax = ubound(phi,1) - adjimax
  jmax = ubound(phi,2) - adjjmax
  kmax = ubound(phi,3) - adjkmax
!
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
  dV = dX * dY * dZ
!
! Set up integration
!
  dmass = ZERO
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
!
! Integration sum
!
           dmass = dmass + dV * u0(i,j,k) * &
                (alpha(i,j,k) + ONE) * exp(SIX * phi(i,j,k)) * rho_b(i,j,k)
        end do
     end do
  end do
  return
end subroutine irr_rest_mass
!
!-----------------------------------------------------------------------------
!
! Compute Aij using formula (13) of gr-qc/9709026
!
!-----------------------------------------------------------------------------
subroutine kset_bns(ex, X, Y, Z,  &
     Axx, Axy, Axz, Ayy, Ayz, Azz, betax, betay, betaz,  &
     phi, lapse, Symmetry)
  implicit none
  interface
     subroutine gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_oct
     subroutine gderivs_eq(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_eq
  end interface
!
! Input parameters:
!
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: Axx, Axy, Axz, Ayy, Ayz, Azz
  real*8, dimension(ex(1),ex(2),ex(3))        :: betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3))        :: phi, lapse
  integer                                            :: Symmetry
!
! Other variables:
!
  integer                                     :: imin, jmin, kmin
  real*8, dimension(ex(1),ex(2),ex(3))        :: betaxx, betaxy, betaxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: betayx, betayy, betayz
  real*8, dimension(ex(1),ex(2),ex(3))        :: betazx, betazy, betazz
  real*8, dimension(ex(1),ex(2),ex(3))        :: div_beta
  real*8                                      :: dX, dY, dZ
  real*8                                      :: psi, kxxout
  real*8                                      :: F1o3, F1o6, ONE, TWO, FOUR, ZERO
  real*8                                      :: F2o3, SIX, EIGHT, HALF, PI
  real*8                                      :: SYM, ANTI
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  parameter ( ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0, F1o6 = 1.D0/6.D0 )
  parameter ( ZERO = 0.D0, F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0, SIX = 6.D0 )
  parameter ( SYM = 1.D0, ANTI = - 1.D0, EIGHT = 8.D0, HALF = 0.5D0 )
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3)
  PI = acos(-ONE)
  imin = lbound(phi,1)
  jmin = lbound(phi,2)
  kmin = lbound(phi,3)
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
!-----------------------------------------------------------------------------
! Compute first derivatives of the shift
!-----------------------------------------------------------------------------
  if(Symmetry==OCTANT) then
     call gderivs_oct(ex,betax,betaxx,betaxy,betaxz,dX,dY,dZ,ANTI,SYM, SYM)
     call gderivs_oct(ex,betay,betayx,betayy,betayz,dX,dY,dZ,SYM, ANTI,SYM)
     call gderivs_oct(ex,betaz,betazx,betazy,betazz,dX,dY,dZ,SYM, SYM, ANTI)
  else if (Symmetry == EQUATORIAL) then
     call gderivs_eq(ex,betax,betaxx,betaxy,betaxz,dX,dY,dZ,ANTI,SYM, SYM)
     call gderivs_eq(ex,betay,betayx,betayy,betayz,dX,dY,dZ,SYM, ANTI,SYM)
     call gderivs_eq(ex,betaz,betazx,betazy,betazz,dX,dY,dZ,SYM, SYM, ANTI)
  else if (Symmetry == PI_SYMM) then
     write(*,*) "PI_SYMMETRY not supported!"
     stop
  end if
  div_beta = betaxx + betayy + betazz
!----------------------------------------------------------
! Output
!----------------------------------------------------------
  Axx = (TWO * betaxx - F2o3 * div_beta) / (TWO * (lapse + ONE) )
  Axy = (betaxy + betayx) / (TWO * (lapse + ONE) )
  Axz = (betaxz + betazx) / (TWO * (lapse + ONE) )
  Ayy = (TWO * betayy - F2o3 * div_beta) / (TWO * (lapse + ONE) )
  Ayz = (betayz + betazy) / (TWO * (lapse + ONE) )
  Azz = (TWO * betazz - F2o3 * div_beta) / (TWO * (lapse + ONE) )
end subroutine kset_bns
