!
!-----------------------------------------------------------------------------
!
! Integrand for the EM field energy: (E^2 + B^2)/8pi
!
!-----------------------------------------------------------------------------
!
subroutine em_energy_excise_disk(ext, demenergy, &
     X, Y, Z, phi, gxx,gxy,gxz,gyy,gyz,gzz, Ex,Ey,Ez, &
     Bx,By,Bz, Symmetry, mskf, x_plus, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                        :: ext
  real*8, dimension(ext(1),ext(2),ext(3))      :: X,Y,Z, mskf
  real*8, dimension(ext(1),ext(2),ext(3))      :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))      :: Ex,Ey,Ez, Bx,By,Bz
  integer                                      :: adjimin, adjjmin, adjkmin
  integer                                      :: adjimax, adjjmax, adjkmax
  integer                                      :: Symmetry
  real*8                                       :: x_plus
! output:
  real*8                                       :: demenergy
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,pi,FOUR,ONE, temp
  integer                            :: i,j,k, ilow, klow
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

  ilow = imin
  klow = kmin

!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  demenergy = 0.D0
  do k = klow, kmax
     do j = jmin, jmax
        do i = ilow, imax
           if (Symmetry == AXISYM) then
              temp = exp(10.d0*phi(i,j,k)) * &
                   ( gxx(i,j,k)*(Ex(i,j,k)**2 + Bx(i,j,k)**2) + &
                   2.d0*gxy(i,j,k)*(Ex(i,j,k)*Ey(i,j,k) + Bx(i,j,k)*By(i,j,k)) + &
                   2.d0*gxz(i,j,k)*(Ex(i,j,k)*Ez(i,j,k) + Bx(i,j,k)*Bz(i,j,k)) + &
                   gyy(i,j,k)*(Ey(i,j,k)**2 + By(i,j,k)**2) + &
                   2.d0*gyz(i,j,k)*(Ey(i,j,k)*Ez(i,j,k) + By(i,j,k)*Bz(i,j,k)) + &
                   gzz(i,j,k)*(Ez(i,j,k)**2 + Bz(i,j,k)**2) ) * abs(X(i,1,1))
              temp = temp*mskf(i,j,k)

              demenergy = demenergy + temp

           else
              temp = exp(10.d0*phi(i,j,k)) * &
                   ( gxx(i,j,k)*(Ex(i,j,k)**2 + Bx(i,j,k)**2) + &
                   2.d0*gxy(i,j,k)*(Ex(i,j,k)*Ey(i,j,k) + Bx(i,j,k)*By(i,j,k)) + &
                   2.d0*gxz(i,j,k)*(Ex(i,j,k)*Ez(i,j,k) + Bx(i,j,k)*Bz(i,j,k)) + &
                   gyy(i,j,k)*(Ey(i,j,k)**2 + By(i,j,k)**2) + &
                   2.d0*gyz(i,j,k)*(Ey(i,j,k)*Ez(i,j,k) + By(i,j,k)*Bz(i,j,k)) + &
                   gzz(i,j,k)*(Ez(i,j,k)**2 + Bz(i,j,k)**2) )

                 demenergy = demenergy + temp
           end if
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  demenergy = demenergy * dV/(8.d0*PI)
  return
end subroutine em_energy_excise_disk
