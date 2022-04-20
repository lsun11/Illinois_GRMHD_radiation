!-----------------------------------------------------------------------------
!
! Compute b^2/rho0 and b^2/P
!
!-----------------------------------------------------------------------------
!
subroutine compute_b2p_b2rho(ext,sbt,sbx,sby,sbz,rho_b,P,alpham1,betax,betay, &
                betaz, phi, gxx,gxy,gxz,gyy,gyz,gzz, b2rho,b2p,b2,Symmetry)
 implicit none
 integer, dimension(3)                                  :: ext
 real*8, dimension(ext(1),ext(2),ext(3))                :: sbt,sbx,sby,sbz
 real*8, dimension(ext(1),ext(2),ext(3))                :: rho_b,P,alpham1,phi
 real*8, dimension(ext(1),ext(2),ext(3))                :: betax,betay,betaz
 real*8, dimension(ext(1),ext(2),ext(3))                :: gxx,gxy,gxz,gyy,gyz,gzz
 real*8, dimension(ext(1),ext(2),ext(3))                :: b2rho,b2p,b2
 integer                                                :: Symmetry
 integer, parameter                                        :: AXISYM = 4
!
 b2 = ( -((alpham1+1.d0)*sbt)**2 + exp(4.d0*phi)*( &
        gxx*(sbx + betax*sbt)**2 + &
        2.d0*gxy*(sbx + betax*sbt)*(sby + betay*sbt) + &
        2.d0*gxz*(sbx + betax*sbt)*(sbz + betaz*sbt) + &
        gyy*(sby + betay*sbt)**2 + &
        2.d0*gyz*(sby + betay*sbt)*(sbz + betaz*sbt) + &
        gzz*(sbz + betaz*sbt)**2 )) / (4.d0*acos(-1.d0)) / (alpham1+1.d0)**2
 b2rho = b2/rho_b
 b2p   = b2/P
 if (Symmetry==AXISYM) then
    b2rho(:,1,:) = 0.d0
    b2rho(:,3,:) = 0.d0
    b2p(:,1,:) = 0.d0
    b2p(:,3,:) = 0.d0
 end if
end subroutine compute_b2p_b2rho
subroutine compute_b2p_b2rho_old(ext,sbt,sbx,sby,sbz,rho_b,P,alpham1,betax,betay, &
                betaz, phi, gxx,gxy,gxz,gyy,gyz,gzz, b2rho,b2p,b2,Epsilon,n)
 implicit none
 integer, dimension(3)                                  :: ext
 real*8, dimension(ext(1),ext(2),ext(3))                :: sbt,sbx,sby,sbz
 real*8, dimension(ext(1),ext(2),ext(3))                :: rho_b,P,alpham1,phi
 real*8, dimension(ext(1),ext(2),ext(3))                :: betax,betay,betaz
 real*8, dimension(ext(1),ext(2),ext(3))                :: gxx,gxy,gxz,gyy,gyz,gzz
 real*8, dimension(ext(1),ext(2),ext(3))                :: b2rho,b2p,b2,Epsilon
 real*8                                                 :: n
 real*8, parameter :: TINY = 1.d-10
!
 b2 = ( -((alpham1+1.d0)*sbt)**2 + exp(4.d0*phi)*( &
        gxx*(sbx + betax*sbt)**2 + &
        2.d0*gxy*(sbx + betax*sbt)*(sby + betay*sbt) + &
        2.d0*gxz*(sbx + betax*sbt)*(sbz + betaz*sbt) + &
        gyy*(sby + betay*sbt)**2 + &
        2.d0*gyz*(sby + betay*sbt)*(sbz + betaz*sbt) + &
        gzz*(sbz + betaz*sbt)**2 )) / (4.d0*acos(-1.d0)) / (alpham1+1.d0)**2
 Epsilon = b2 + rho_b + (n+1)*P
 b2rho = b2/(rho_b+TINY)
 b2p   = b2/(P+TINY)
end subroutine compute_b2p_b2rho_old
!
!-----------------------------------------------------------------------------
!
! Integrand for the EM field energy: (E^2 + B^2)/8pi
!
!-----------------------------------------------------------------------------
!
subroutine em_energy(ext, demenergy, &
     X, Y, Z, phi, gxx,gxy,gxz,gyy,gyz,gzz, Ex,Ey,Ez, &
     Bx,By,Bz, Symmetry, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                        :: ext
  real*8, dimension(ext(1),ext(2),ext(3))      :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3))      :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))      :: Ex,Ey,Ez, Bx,By,Bz
  integer                                      :: adjimin, adjjmin, adjkmin
  integer                                      :: adjimax, adjjmax, adjkmax
  integer                                      :: Symmetry
! output:
  real*8                                       :: demenergy
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - ONE )
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
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  demenergy = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           if (Symmetry == AXISYM) then
              demenergy = demenergy + exp(10.d0*phi(i,j,k)) * &
                ( gxx(i,j,k)*(Ex(i,j,k)**2 + Bx(i,j,k)**2) + &
                2.d0*gxy(i,j,k)*(Ex(i,j,k)*Ey(i,j,k) + Bx(i,j,k)*By(i,j,k)) + &
                2.d0*gxz(i,j,k)*(Ex(i,j,k)*Ez(i,j,k) + Bx(i,j,k)*Bz(i,j,k)) + &
                gyy(i,j,k)*(Ey(i,j,k)**2 + By(i,j,k)**2) + &
                2.d0*gyz(i,j,k)*(Ey(i,j,k)*Ez(i,j,k) + By(i,j,k)*Bz(i,j,k)) + &
                gzz(i,j,k)*(Ez(i,j,k)**2 + Bz(i,j,k)**2) ) * abs(X(i,1,1))
           else
              demenergy = demenergy + exp(10.d0*phi(i,j,k)) * &
                ( gxx(i,j,k)*(Ex(i,j,k)**2 + Bx(i,j,k)**2) + &
                2.d0*gxy(i,j,k)*(Ex(i,j,k)*Ey(i,j,k) + Bx(i,j,k)*By(i,j,k)) + &
                2.d0*gxz(i,j,k)*(Ex(i,j,k)*Ez(i,j,k) + Bx(i,j,k)*Bz(i,j,k)) + &
                gyy(i,j,k)*(Ey(i,j,k)**2 + By(i,j,k)**2) + &
                2.d0*gyz(i,j,k)*(Ey(i,j,k)*Ez(i,j,k) + By(i,j,k)*Bz(i,j,k)) + &
                gzz(i,j,k)*(Ez(i,j,k)**2 + Bz(i,j,k)**2) )
           end if
        end do
     end do
  end do
  demenergy = demenergy * dV/(8.d0*PI)
  return
end subroutine em_energy
!
! Integrand for the EM field energy: (E^2 + B^2)/8pi
!
!-----------------------------------------------------------------------------
!
subroutine em_energy2(ext, demenergy, &
     X, Y, Z, alpham1, vx,vy,vz, shiftx, shifty, shiftz, &
     phi, gxx,gxy,gxz,gyy,gyz,gzz, &
     sbt, sbx, sby, sbz, Symmetry, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                        :: ext
  real*8, dimension(ext(1),ext(2),ext(3))      :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3))      :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))      :: u0, alpham1, shiftx, shifty, shiftz
  real*8, dimension(ext(1),ext(2),ext(3))      :: sbt, sbx, sby, sbz
  real*8, dimension(ext(1),ext(2),ext(3))      :: er, el, au0m1, vx, vy, vz
  integer                                      :: adjimin, adjjmin, adjkmin
  integer                                      :: adjimax, adjjmax, adjkmax
  integer                                      :: Symmetry
! output:
  real*8                                       :: demenergy
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  real*8                             :: b2
  real*8, parameter                     :: fac = 0.99d0
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - ONE )
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
 dX = X(imin+1,1,1)-X(imin,1,1)
 dY = Y(1,jmin+1,1)-Y(1,jmin,1)
 dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  end if
! Compute al*u0-1
  er = exp(4.d0*phi)*(gxx*(vx + shiftx)**2 + &
       2.d0*gxy*(vx + shiftx)*(vy + shifty) +         &
       2.d0*gxz*(vx + shiftx)*(vz + shiftz) +         &
       gyy*(vy + shifty)**2 +                           &
       2.d0*gyz*(vy + shifty)*(vz + shiftz) +         &
       gzz*(vz + shiftz)**2 )/(alpham1+1.d0)**2
! *** Check for superluminal velocity ***
  where (er .gt. 1.d0)
    vx = (vx + shiftx)*sqrt(fac/er)-shiftx
    vy = (vy + shifty)*sqrt(fac/er)-shifty
    vz = (vz + shiftz)*sqrt(fac/er)-shiftz
    er = fac
  end where
! ***************************************
  el = sqrt(1.d0-er)
  au0m1 = er/el/(1.d0+el)
  u0 = (au0m1+1.d0)/(alpham1+1.d0)
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  demenergy = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           b2 = -(alpham1(i,j,k)+1.d0)**2*(sbt(i,j,k))**2 + &
                exp(4.d0*phi(i,j,k))*( &
                gxx(i,j,k)*(shiftx(i,j,k)*sbt(i,j,k)+sbx(i,j,k))**2 + &
                2.d0*gxy(i,j,k)*(shiftx(i,j,k)*sbt(i,j,k)+sbx(i,j,k))*(shifty(i,j,k)*sbt(i,j,k)+sby(i,j,k)) + &
                2.d0*gxz(i,j,k)*(shiftx(i,j,k)*sbt(i,j,k)+sbx(i,j,k))*(shiftz(i,j,k)*sbt(i,j,k)+sbz(i,j,k)) + &
                gyy(i,j,k)*(shifty(i,j,k)*sbt(i,j,k)+sby(i,j,k))**2 + &
                2.d0*gyz(i,j,k)*(shifty(i,j,k)*sbt(i,j,k)+sby(i,j,k))*(shiftz(i,j,k)*sbt(i,j,k)+sbz(i,j,k)) + &
                gzz(i,j,k)*(shiftz(i,j,k)*sbt(i,j,k)+sbz(i,j,k))**2)
           b2 = b2/(alpham1(i,j,k)+1.d0)
           if (Symmetry == AXISYM) then
              demenergy = demenergy + u0(i,j,k)*exp(6.d0*phi(i,j,k))*b2*abs(X(i,1,1))
           else
              demenergy = demenergy + u0(i,j,k)*exp(6.d0*phi(i,j,k))*b2
           end if
        end do
     end do
  end do
  demenergy = demenergy * dV/(8.d0*pi)
  return
end subroutine em_energy2
!
! Integrate |b_{\phi}| over proper volume
!
subroutine sb_phi(ex, ib_phi, &
     X, Y, Z, alpham1, phi, sbt, sbx, sby, sbz, &
     betax, betay, betaz, &
     gxy, gyy, gyz, Symmetry, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: phi,alpham1
  real*8, dimension(ex(1),ex(2),ex(3))        :: sbx, sby, sbz, sbt
  real*8, dimension(ex(1),ex(2),ex(3))        :: gxy, gyy, gyz
  real*8, dimension(ex(1),ex(2),ex(3))        :: betax, betay, betaz
  real*8                                      :: b_phi
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
! output:
  real*8                                      :: ib_phi
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - one )
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
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  ib_phi = 0.d0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           b_phi = X(i,1,1)*exp(4.0*phi(i,j,k))*( &
                   gxy(i,j,k)*(betax(i,j,k)*sbt(i,j,k) + sbx(i,j,k)) + &
                   gyy(i,j,k)*(betay(i,j,k)*sbt(i,j,k) + sby(i,j,k)) + &
                   gyz(i,j,k)*(betaz(i,j,k)*sbt(i,j,k) + sbz(i,j,k)))
           if (Symmetry == AXISYM) then
              ib_phi = ib_phi + X(i,1,1)*exp(6.d0*phi(i,j,k))*abs(b_phi)/(1.d0+alpham1(i,j,k))
           else
              ib_phi = ib_phi + exp(6.d0*phi(i,j,k))*abs(b_phi)/(1.d0+alpham1(i,j,k))
           end if
        end do
     end do
  end do
  ib_phi = ib_phi * dV/sqrt(4.d0*pi)
  return
end subroutine sb_phi
!
! Integrate |B_{\phi}| over proper volume
!
subroutine b_phi(ex, iB_phi, rho_star, &
     X, Y, Z, phi, Bx, By, Bz, gxy, gyy, gyz, Symmetry, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: Bx,By,Bz,phi, rho_star
  real*8, dimension(ex(1),ex(2),ex(3))        :: gxy, gyy, gyz
  real*8                                      :: Bsubphi
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
! output:
  real*8                                      :: iB_phi
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - one )
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
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  iB_phi = 0.d0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           Bsubphi = X(i,1,1)*exp(4.d0*phi(i,j,k))*(gxy(i,j,k)*Bx(i,j,k) + &
                            gyy(i,j,k)*By(i,j,k) + gyz(i,j,k)*Bz(i,j,k))
           if (Symmetry == AXISYM) then
              iB_phi = iB_phi + X(i,1,1)*rho_star(i,j,k)*abs(Bsubphi)
           else
              iB_phi = iB_phi + rho_star(i,j,k)*abs(Bsubphi)
           end if
        end do
     end do
  end do
  iB_phi = iB_phi * dV
  return
end subroutine b_phi
!
! Compute the Alfven speed
!
subroutine compute_va(ext,sbt,sbx,sby,sbz,rho_b,P,alpham1,phi, &
        betax,betay,betaz,gxx,gxy,gxz,gyy,gyz,gzz,va,n)
 implicit none
 integer, dimension(3)                                  :: ext
 real*8, dimension(ext(1),ext(2),ext(3))                :: sbt,sbx,sby,sbz
 real*8, dimension(ext(1),ext(2),ext(3))                :: rho_b,P,alpham1,phi
 real*8, dimension(ext(1),ext(2),ext(3))                :: betax,betay,betaz
 real*8, dimension(ext(1),ext(2),ext(3))                :: gxx,gxy,gxz,gyy,gyz,gzz
 real*8, dimension(ext(1),ext(2),ext(3))                :: va,b2
 real*8                                                        :: n
!
 b2 = ( -((alpham1+1.d0)*sbt)**2 + exp(4.d0*phi)*( &
        gxx*(sbx + betax*sbt)**2 + &
        2.d0*gxy*(sbx + betax*sbt)*(sby + betay*sbt) + &
        2.d0*gxz*(sbx + betax*sbt)*(sbz + betaz*sbt) + &
        gyy*(sby + betay*sbt)**2 + &
        2.d0*gyz*(sby + betay*sbt)*(sbz + betaz*sbt) + &
        gzz*(sbz + betaz*sbt)**2 )) / (4.d0*acos(-1.d0)) / (alpham1+1.d0)**2
 va = sqrt( b2/(rho_b + (n+1.d0)*P + b2) )
end subroutine compute_va
!-------------------------------------------------------------------
!  Compute F(x) = Integrate[Bz(u,0),{u,0,x}] (Assuming axisymmetry)
!-------------------------------------------------------------------
!
subroutine int_bz(ext,X,Z,Bz,F)
 implicit none
 integer, dimension(3)                                  :: ext
 real*8, dimension(ext(1),ext(2),ext(3))                :: X,Z
 real*8, dimension(ext(1),ext(2),ext(3))                :: Bz
 real*8, dimension(ext(1))                              :: F
 integer                                                :: i,j,k,imin,imax
 integer                                                :: kmin
 real*8                                                 :: dX, dZ
!
 imin = lbound(Bz,1)
 imax = ubound(Bz,1)
 kmin = lbound(Bz,3)
 dX = X(imin+1,1,1)-X(imin,1,1)
 dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
 if (Z(1,1,kmin) .gt. dZ) then
    F = 0.d0
 else
    j = 2
    k = 1
    F(imin) = Bz(imin,j,k)
    do i=imin+1,imax
       F(i) = F(i-1) + Bz(i,j,k)
    end do
    F = F * dX
 end if
end subroutine int_bz
!
!-----------------------------------------------------------------------------
!
! Compute sqrt(B \sqrt{gamma}) = (\gamma_{ij} B^i B^j * \gamma)^{1/4}
!
!-----------------------------------------------------------------------------
!
subroutine b2_1o4(ext, B2p1o4, phi, gxx,gxy,gxz,gyy,gyz,gzz,Bx,By,Bz)
  implicit none
! Input parameters:
  integer, dimension(3)                        :: ext
  real*8, dimension(ext(1),ext(2),ext(3))      :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))      :: Bx,By,Bz
! output:
  real*8, dimension(ext(1),ext(2),ext(3))      :: B2p1o4
!
  B2p1o4 = exp(4.d0*phi) * ( gxx*Bx**2 + 2.d0*gxy*Bx*By + 2.d0*gxz*Bx*Bz + &
                gyy*By**2 + 2.d0*gyz*By*Bz + gzz*Bz**2)**(0.25d0)
end subroutine b2_1o4
subroutine compute_wavelength(ext,sbt,sbx,sby,sbz,rho_b,P,alpham1,phi, &
        betax,betay,betaz,gxx,gxy,gxz,gyy,gyz,gzz,n,X,Z,vy)
 implicit none
 integer, dimension(3)                                  :: ext
 real*8, dimension(ext(1),ext(2),ext(3))                :: X,Z
 real*8, dimension(ext(1),ext(2),ext(3))                :: sbt,sbx,sby,sbz
 real*8, dimension(ext(1),ext(2),ext(3))                :: rho_b,P,alpham1,phi
 real*8, dimension(ext(1),ext(2),ext(3))                :: betax,betay,betaz
 real*8, dimension(ext(1),ext(2),ext(3))                :: gxx,gxy,gxz,gyy,gyz,gzz
 real*8, dimension(ext(1),ext(2),ext(3))                :: b2,vy, va
 real*8                                                 :: n, dZ, omg, tpi,lam,b2op
 integer                                                :: i,imin,imax,kmin
!
 imin = lbound(P,1)
 imax = ubound(P,1)
 kmin = lbound(P,1)
 dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
 tpi = 2.d0*acos(-1.d0)
 b2 = ( -((alpham1+1.d0)*sbt)**2 + exp(4.d0*phi)*( &
        gxx*(sbx + betax*sbt)**2 + &
        2.d0*gxy*(sbx + betax*sbt)*(sby + betay*sbt) + &
        2.d0*gxz*(sbx + betax*sbt)*(sbz + betaz*sbt) + &
        gyy*(sby + betay*sbt)**2 + &
        2.d0*gyz*(sby + betay*sbt)*(sbz + betaz*sbt) + &
        gzz*(sbz + betaz*sbt)**2 )) / (4.d0*acos(-1.d0)) / (alpham1+1.d0)**2
 va = sqrt( b2/(rho_b + (n+1.d0)*P + b2) )
 if (Z(1,1,kmin) .lt. dZ) then
    do i=imin,imax
       omg = vy(i,2,kmin)/X(i,1,1)
       lam = tpi/omg*va(i,2,kmin)
       b2op = b2(i,2,kmin)/2.d0/P(i,2,kmin)
       write(*,*) X(i,1,1),lam,b2op
    end do
 end if
end subroutine compute_wavelength
!
! Find the wavelength of the fastest-growing mode
!
subroutine diagnose_mri(ext,X,rho_b,P,Bz,vy,n)
 implicit none
 integer, dimension(3)                                  :: ext
 real*8, dimension(ext(1),ext(2),ext(3))                :: X
 real*8, dimension(ext(1),ext(2),ext(3))                :: rho_b,P,Bz,vy
 real*8                                                 :: omega,kappa2,s2
 real*8                                                 :: lambda,va, kvamax
 real*8                                                 :: omega_max,tau
 real*8                                                 :: n, PI, dX2, cs
 real*8                                                 :: drho_dx, rat
 integer                                                :: i,j,k,imin,jmin,kmin
 integer                                                :: imax,jmax,kmax
 PI = acos(-1.0)
!
 imin = lbound(P,1)
 imax = ubound(P,1)
 jmin = lbound(P,2)
 jmax = ubound(P,2)
 kmin = lbound(P,3)
 kmax = ubound(P,3)
 dX2 = 2.d0*(X(imin+1,1,1)-X(imin,1,1))
 j = jmin+1
 k = kmin
 do i=imin+1,imax-1
    kappa2 = (vy(i+1,j,k)**2*X(i+1,1,1)**2 - vy(i-1,j,k)**2*X(i-1,1,1)**2)/dX2
    kappa2 = kappa2/(X(i,1,1)**3)
    omega = vy(i,j,k)/X(i,1,1)
    s2 = 4.d0*omega**2 - kappa2
    kvamax = s2/4.d0*(s2+2.d0*kappa2)/(s2+kappa2)
    if (kvamax .lt. 0.d0) kvamax = 0.d0
    kvamax = sqrt(max(kvamax,1.d-100))
    va = Bz(i,j,k)/sqrt(4.d0*PI*(rho_b(i,j,k)+(n+1.0)*P(i,j,k)))
    lambda = 2.d0*PI*abs(va)/kvamax
    omega_max = sqrt( max( s2*s2/4.d0/(s2+kappa2), 1.d-100) )
    tau = 2.d0*PI/omega_max
    if (rho_b(i,j,k) .ne. 0.d0) then
       drho_dx = (P(i+1,j,k) - P(i-1,j,k))/dX2
       rat = drho_dx*lambda/P(i,j,k)
    else
       rat = 0.d0
    end if
    cs = sqrt((1.d0+1.d0/n)*P(i,j,k)/(rho_b(i,j,k)+(n+1.0)*P(i,j,k)))
    write(*,100) X(i,1,1),rho_b(i,j,k), kappa2, lambda, rat, tau
 end do
100 format(6E13.5)
end subroutine diagnose_mri
subroutine int_xb(ext,X,Y,phi,gxx,gxy,gxz,gyy,gyz,gzz,Bx,By,Bz,xB)
 implicit none
 integer, dimension(3)                                :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y
 real*8, dimension(ext(1),ext(2),ext(3))        :: phi,gxx,gxy,gxz,gyy,gyz,gzz
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bx,By,Bz, xB
 integer                                        :: i,imin,imax,j,jmin,jmax
 integer                                        :: k,kmin,kmax
 real*8                                                :: sq_gpp
!
 imin = lbound(Bx,1)
 imax = ubound(Bx,1)
 jmin = lbound(Bx,2)
 jmax = ubound(Bx,2)
 kmin = lbound(Bx,3)
 kmax = ubound(Bx,3)
 do k = kmin,kmax
    do j = jmin,jmax
       do i = imin,imax
          sq_gpp = exp(2.d0*phi(i,j,k)) * sqrt(gxx(i,j,k)*Y(1,j,1)**2 &
                 - 2.d0*X(i,1,1)*Y(1,j,1) * gxy(i,j,k) + gyy(i,j,k)*X(i,1,1)**2)
          xB(i,j,k) = sq_gpp * exp(2.d0*phi(i,j,k)) * sqrt(  &
             gxx(i,j,k)*Bx(i,j,k)**2 + 2.d0*gxy(i,j,k)*Bx(i,j,k)*By(i,j,k) + &
             2.d0*gxz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) + gyy(i,j,k)*By(i,j,k)**2 + &
             2.d0*gyz(i,j,k)*By(i,j,k)*Bz(i,j,k) + gzz(i,j,k)*Bz(i,j,k)**2 )
       end do
    end do
 end do
end subroutine int_xb
!
!-----------------------------------------------------------------------------
!
! Compute the loss of angular momentum (assuming axisymmetry!) along the 
! cartesian x or z directions due to outflow of matter.
!
!-----------------------------------------------------------------------------
!
subroutine outf_integral(ext, dout_angm, m, &
     X, Y, Z, Symmetry, rho_b, h, vx,vy,vz, &
     alpham1, shiftx, shifty, shiftz, &
     phi, gxx,gxy,gxz,gyy,gyz,gzz, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax, &
     glob_xmax, glob_zmax)
  implicit none
! Input parameters:
  integer, dimension(3)                        :: ext
  real*8, dimension(ext(1),ext(2),ext(3))      :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3))      :: rho_b, h
  real*8, dimension(ext(1),ext(2),ext(3))      :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))      :: u0, alpham1, shiftx, shifty, shiftz
  real*8, dimension(ext(1),ext(2),ext(3))      :: er, el, au0m1, vx, vy, vz
  integer                                      :: adjimin, adjjmin, adjkmin
  integer                                      :: adjimax, adjjmax, adjkmax
  integer                                      :: Symmetry, m
  real*8                                       :: glob_xmax, glob_zmax
! output:
  real*8                                       :: dout_angm
! Other variables:
  real*8                                       :: dX, dY, dZ, dV
  real*8                                       :: pi,FOUR,ONE
  integer                                      :: i,j,k
  integer                                      :: imin, jmin, kmin
  integer                                      :: imax, jmax, kmax
  integer                                      :: AXISYM
  real*8, parameter                               :: fac = 0.99d0
  real*8                                       :: u_y
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - ONE )
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
 dX = X(imin+1,1,1)-X(imin,1,1)
 dY = Y(1,jmin+1,1)-Y(1,jmin,1)
 dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
  end if
  if (m==1) then
     dV = 4.d0*pi*dX
  else if (m==3) then
     dV = 4.d0*pi*dZ
  end if
! Compute al*u0-1
  er = exp(4.d0*phi)*(gxx*(vx + shiftx)**2 + &
       2.d0*gxy*(vx + shiftx)*(vy + shifty) +         &
       2.d0*gxz*(vx + shiftx)*(vz + shiftz) +         &
       gyy*(vy + shifty)**2 +                           &
       2.d0*gyz*(vy + shifty)*(vz + shiftz) +         &
       gzz*(vz + shiftz)**2 )/(alpham1+1.d0)**2
! *** Check for superluminal velocity ***
  where (er .gt. 1.d0)
    vx = (vx + shiftx)*sqrt(fac/er)-shiftx
    vy = (vy + shifty)*sqrt(fac/er)-shifty
    vz = (vz + shiftz)*sqrt(fac/er)-shiftz
    er = fac
  end where
! ***************************************
  el = sqrt(1.d0-er)
  au0m1 = er/el/(1.d0+el)
  u0 = (au0m1+1.d0)/(alpham1+1.d0)
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  dout_angm = 0.D0
  if (m==1) then
     k = kmax
     j = jmin
     if (abs(Z(1,1,kmax) - glob_zmax) .lt. 1.5d0*dZ) then
        do i = imin, imax
           u_y = exp(4.d0*phi(i,j,k))*u0(i,j,k)*(gxy(i,j,k)*(vx(i,j,k)+shiftx(i,j,k)) + &
                                                 gyy(i,j,k)*(vy(i,j,k)+shifty(i,j,k)) + &
                                                 gyz(i,j,k)*(vz(i,j,k)+shiftz(i,j,k)))
           dout_angm = dout_angm + X(i,1,1)**2*(alpham1(i,j,k)+1.d0)* &
                       exp(6.d0*phi(i,j,k))*rho_b(i,j,k)*h(i,j,k)* &
                       u0(i,j,k)*vz(i,j,k)*u_y
        end do
     end if
  else if (m==3) then
     i = imax
     j = jmin
     if (abs(X(imax,1,1) - glob_xmax) .lt. 1.5d0*dX) then
        do k = kmin, kmax
           u_y = exp(4.d0*phi(i,j,k))*u0(i,j,k)*(gxy(i,j,k)*(vx(i,j,k)+shiftx(i,j,k)) + &
                                              gyy(i,j,k)*(vy(i,j,k)+shifty(i,j,k)) + &
                                              gyz(i,j,k)*(vz(i,j,k)+shiftz(i,j,k)))
           dout_angm = dout_angm + X(i,1,1)**2*(alpham1(i,j,k)+1.d0)* &
                       exp(6.d0*phi(i,j,k))*rho_b(i,j,k)*h(i,j,k)* &
                       u0(i,j,k)*vx(i,j,k)*u_y
        end do
     end if
  else
     write(*,*) 'm out of range!'
  end if
  dout_angm = -dout_angm * dV
  return
end subroutine outf_integral
!
!-----------------------------------------------------------------------------
!
! Compute the loss of angular momentum (assuming axisymmetry!) along the 
! cartesian x or z directions due to outflow of matter.
!
!-----------------------------------------------------------------------------
!
subroutine emout_integral(ext, dout_angm, m, &
     X, Y, Z, Symmetry, vx,vy,vz, sbt,sbx,sby,sbz, &
     alpham1, shiftx, shifty, shiftz, &
     phi, gxx,gxy,gxz,gyy,gyz,gzz, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax, &
     glob_xmax, glob_zmax)
  implicit none
! Input parameters:
  integer, dimension(3)                        :: ext
  real*8, dimension(ext(1),ext(2),ext(3))      :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3))      :: sbt, sbx, sby, sbz
  real*8, dimension(ext(1),ext(2),ext(3))      :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))      :: u0, alpham1, shiftx, shifty, shiftz
  real*8, dimension(ext(1),ext(2),ext(3))      :: er, el, au0m1, vx, vy, vz
  integer                                      :: adjimin, adjjmin, adjkmin
  integer                                      :: adjimax, adjjmax, adjkmax
  integer                                      :: Symmetry, m
  real*8                                       :: glob_xmax, glob_zmax
! output:
  real*8                                       :: dout_angm
! Other variables:
  real*8                                       :: dX, dY, dZ, dV
  real*8                                       :: pi,FOUR,ONE
  integer                                      :: i,j,k
  integer                                      :: imin, jmin, kmin
  integer                                      :: imax, jmax, kmax
  integer                                      :: AXISYM
  real*8, parameter                               :: fac = 0.99d0
  real*8                                       :: u_y, b_y, b2
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - ONE )
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
 dX = X(imin+1,1,1)-X(imin,1,1)
 dY = Y(1,jmin+1,1)-Y(1,jmin,1)
 dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
  end if
  if (m==1) then
     dV = 4.d0*pi*dX
  else if (m==3) then
     dV = 4.d0*pi*dZ
  end if
! Compute al*u0-1
  er = exp(4.d0*phi)*(gxx*(vx + shiftx)**2 +          &
       2.d0*gxy*(vx + shiftx)*(vy + shifty) +         &
       2.d0*gxz*(vx + shiftx)*(vz + shiftz) +         &
       gyy*(vy + shifty)**2 +                         &
       2.d0*gyz*(vy + shifty)*(vz + shiftz) +         &
       gzz*(vz + shiftz)**2 )/(alpham1+1.d0)**2
! *** Check for superluminal velocity ***
  where (er .gt. 1.d0)
    vx = (vx + shiftx)*sqrt(fac/er)-shiftx
    vy = (vy + shifty)*sqrt(fac/er)-shifty
    vz = (vz + shiftz)*sqrt(fac/er)-shiftz
    er = fac
  end where
! ***************************************
  el = sqrt(1.d0-er)
  au0m1 = er/el/(1.d0+el)
  u0 = (au0m1+1.d0)/(alpham1+1.d0)
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  dout_angm = 0.D0
  if (m==1) then
     k = kmax
     j = jmin
     if (abs(Z(1,1,kmax)-glob_zmax) .lt. 1.5d0*dZ) then
        do i = imin, imax
           u_y = exp(4.d0*phi(i,j,k))*u0(i,j,k)*(gxy(i,j,k)*(vx(i,j,k)+shiftx(i,j,k)) + &
                                                 gyy(i,j,k)*(vy(i,j,k)+shifty(i,j,k)) + &
                                                 gyz(i,j,k)*(vz(i,j,k)+shiftz(i,j,k)))
           b_y = exp(4.d0*phi(i,j,k))*(gxy(i,j,k)*(sbt(i,j,k)*shiftx(i,j,k)+sbx(i,j,k)) + &
                                    gyy(i,j,k)*(sbt(i,j,k)*shifty(i,j,k)+sby(i,j,k)) + &
                                    gyz(i,j,k)*(sbt(i,j,k)*shiftz(i,j,k)+sbz(i,j,k)))
           b2 = ( -((alpham1(i,j,k)+1.d0)*sbt(i,j,k))**2 + exp(4.d0*phi(i,j,k))*( &
                gxx(i,j,k)*(sbx(i,j,k) + shiftx(i,j,k)*sbt(i,j,k))**2 + &
                2.d0*gxy(i,j,k)*(sbx(i,j,k) + shiftx(i,j,k)*sbt(i,j,k))*(sby(i,j,k) + shifty(i,j,k)*sbt(i,j,k)) + &
                2.d0*gxz(i,j,k)*(sbx(i,j,k) + shiftx(i,j,k)*sbt(i,j,k))*(sbz(i,j,k) + shiftz(i,j,k)*sbt(i,j,k)) + &
                gyy(i,j,k)*(sby(i,j,k) + shifty(i,j,k)*sbt(i,j,k))**2 + &
                2.d0*gyz(i,j,k)*(sby(i,j,k) + shifty(i,j,k)*sbt(i,j,k))*(sbz(i,j,k) + shiftz(i,j,k)*sbt(i,j,k)) + &
                gzz(i,j,k)*(sbz(i,j,k) + shiftz(i,j,k)*sbt(i,j,k))**2 )) / (4.d0*acos(-1.d0)) / (alpham1(i,j,k)+1.d0)**2
!!$           if (b2 .gt. 0.d0) then 
!!$              write(*,*) b2 = ,  b2,  at i = , i
!!$           end if 
!!$            
!!$           if (sbz(i,j,k) .gt. 0.d0) then
!!$              write(*,*) sbz = , sbz(i,j,k),  at i = , i
!!$           end if
           dout_angm = dout_angm + X(i,1,1)**2*(alpham1(i,j,k)+1.d0)* &
                exp(6.d0*phi(i,j,k))*(b2*u0(i,j,k)*vz(i,j,k)*u_y - sbz(i,j,k)*b_y / (4.d0*acos(-1.d0)) )
        end do
     end if
  else if (m==3) then
     i = imax
     j = jmin
     if (abs(X(imax,1,1) - glob_xmax) .lt. 1.5d0*dX) then
        do k = kmin, kmax
           u_y = exp(4.d0*phi(i,j,k))*u0(i,j,k)*(gxy(i,j,k)*(vx(i,j,k)+shiftx(i,j,k)) + &
                                              gyy(i,j,k)*(vy(i,j,k)+shifty(i,j,k)) + &
                                              gyz(i,j,k)*(vz(i,j,k)+shiftz(i,j,k)))
           b_y = exp(4.d0*phi(i,j,k))*(gxy(i,j,k)*(sbt(i,j,k)*shiftx(i,j,k)+sbx(i,j,k)) + &
                                    gyy(i,j,k)*(sbt(i,j,k)*shifty(i,j,k)+sby(i,j,k)) + &
                                    gyz(i,j,k)*(sbt(i,j,k)*shiftz(i,j,k)+sbz(i,j,k)))
           b2 = ( -((alpham1(i,j,k)+1.d0)*sbt(i,j,k))**2 + exp(4.d0*phi(i,j,k))*( &
                gxx(i,j,k)*(sbx(i,j,k) + shiftx(i,j,k)*sbt(i,j,k))**2 + &
                2.d0*gxy(i,j,k)*(sbx(i,j,k) + shiftx(i,j,k)*sbt(i,j,k))*(sby(i,j,k) + shifty(i,j,k)*sbt(i,j,k)) + &
                2.d0*gxz(i,j,k)*(sbx(i,j,k) + shiftx(i,j,k)*sbt(i,j,k))*(sbz(i,j,k) + shiftz(i,j,k)*sbt(i,j,k)) + &
                gyy(i,j,k)*(sby(i,j,k) + shifty(i,j,k)*sbt(i,j,k))**2 + &
                2.d0*gyz(i,j,k)*(sby(i,j,k) + shifty(i,j,k)*sbt(i,j,k))*(sbz(i,j,k) + shiftz(i,j,k)*sbt(i,j,k)) + &
                gzz(i,j,k)*(sbz(i,j,k) + shiftz(i,j,k)*sbt(i,j,k))**2 )) / (4.d0*acos(-1.d0)) / (alpham1(i,j,k)+1.d0)**2
           dout_angm = dout_angm + X(i,1,1)**2*(alpham1(i,j,k)+1.d0)* &
                exp(6.d0*phi(i,j,k))*(b2*u0(i,j,k)*vx(i,j,k)*u_y - sbx(i,j,k)*b_y / (4.d0*acos(-1.d0)) )
        end do
     end if
  else
     write(*,*) 'm out of range!'
  end if
  dout_angm = -dout_angm * dV
  return
end subroutine emout_integral
!
! Assuming axisymmetry, find the angular momentum contributed by the matter
!
subroutine jtot_axi_hydro(ex, dQ, &
     X, Y, Z, Symmetry, st_y, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: st_y
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
! output:
  real*8                                      :: dQ
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - one )
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(st_y,1) - adjimin
  jmin = lbound(st_y,2) - adjjmin
  kmin = lbound(st_y,3) - adjkmin
  imax = ubound(st_y,1) - adjimax
  jmax = ubound(st_y,2) - adjjmax
  kmax = ubound(st_y,3) - adjkmax
!
 dX = X(imin+1,1,1)-X(imin,1,1)
 dY = Y(1,jmin+1,1)-Y(1,jmin,1)
 dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  else
     dQ = 0.d0
     return
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
 dQ = 0.d0
 j = jmin
  do k = kmin, kmax
     do i = imin, imax
        dQ = dQ + X(i,1,1)**2 * st_y(i,j,k)
     end do
  end do
  dQ = dQ * dV
  return
end subroutine jtot_axi_hydro
!
! Assuming axisymmetry, find the angular momentum contributed by the EM fields
!
subroutine jtot_axi_em(ext, dQ, &
     X, Y, Z, Symmetry, phi, Bx, Bz, Ex, Ez, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ext
  real*8, dimension(ext(1),ext(2),ext(3))     :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3))     :: phi, Bx, Bz, Ex, Ez
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
! output:
  real*8                                      :: dQ
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - one )
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
 dX = X(imin+1,1,1)-X(imin,1,1)
 dY = Y(1,jmin+1,1)-Y(1,jmin,1)
 dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  else
     dQ = 0.d0
     return
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
 dQ = 0.d0
 j = jmin
  do k = kmin, kmax
     do i = imin, imax
        dQ = dQ + X(i,1,1)**2 * exp(12.d0*phi(i,j,k))*(Ez(i,j,k)*Bx(i,j,k)- &
                                                   Ex(i,j,k)*Bz(i,j,k))/(8.d0*PI)
     end do
  end do
  dQ = dQ * dV
  return
end subroutine jtot_axi_em
!
! Assuming axisymmetry, calculate the time derivative of the Ixx
! for quadrupole approx wave extraction in 
! Shibata & Sekiguchi 68, 104020 (2003)
!
subroutine dtixx(ex, dQ, X, Y, Z, rho_star, vx, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax, Symmetry)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_star, vx
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
! output:
  real*8                                      :: dQ
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - one )
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(rho_star,1) - adjimin
  jmin = lbound(rho_star,2) - adjjmin
  kmin = lbound(rho_star,3) - adjkmin
  imax = ubound(rho_star,1) - adjimax
  jmax = ubound(rho_star,2) - adjjmax
  kmax = ubound(rho_star,3) - adjkmax
!
 dX = X(imin+1,1,1)-X(imin,1,1)
 dY = Y(1,jmin+1,1)-Y(1,jmin,1)
 dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  else
     dQ = 0.d0
     return
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
 dQ = 0.d0
 j = jmin
  do k = kmin, kmax
     do i = imin, imax
        dQ = dQ + 2.d0*X(i,1,1)*rho_star(i,j,k)*X(i,1,1)*vx(i,j,k)
     end do
  end do
  dQ = dQ * dV
  return
end subroutine dtixx
!
! Assuming axisymmetry, calculate the time derivative of the Ixx
! for quadrupole approx wave extraction in 
! Shibata & Sekiguchi 68, 104020 (2003)
!
subroutine dtizz(ex, dQ, X, Y, Z, rho_star, vz, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax, Symmetry)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_star, vz
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
! output:
  real*8                                      :: dQ
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - one )
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(rho_star,1) - adjimin
  jmin = lbound(rho_star,2) - adjjmin
  kmin = lbound(rho_star,3) - adjkmin
  imax = ubound(rho_star,1) - adjimax
  jmax = ubound(rho_star,2) - adjjmax
  kmax = ubound(rho_star,3) - adjkmax
!
 dX = X(imin+1,1,1)-X(imin,1,1)
 dY = Y(1,jmin+1,1)-Y(1,jmin,1)
 dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  else
     dQ = 0.d0
     return
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
 dQ = 0.d0
 j = jmin
  do k = kmin, kmax
     do i = imin, imax
        dQ = dQ + 2.d0*X(i,1,1)*rho_star(i,j,k)*Z(1,1,k)*vz(i,j,k)
     end do
  end do
  dQ = dQ * dV
  return
end subroutine dtizz
