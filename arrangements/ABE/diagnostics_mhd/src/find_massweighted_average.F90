!-----------------------------------------------------------------------------
!
! Integrate R*rho_star to get <R>
!
!-----------------------------------------------------------------------------
subroutine wavg(ex, dR, &
     X, Y, Z, w, f, Symmetry, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: w, f
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
! output:
  real*8                                      :: dR
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,R
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer, parameter                 :: AXISYM=4
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(f,1) - adjimin
  jmin = lbound(f,2) - adjjmin
  kmin = lbound(f,3) - adjkmin
  imax = ubound(f,1) - adjimax
  jmax = ubound(f,2) - adjjmax
  kmax = ubound(f,3) - adjkmax
!
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
  dV = dX * dY * dZ
   if (Symmetry==AXISYM) then
     dV=dX*dZ*4.d0*acos(-1.d0)
     jmin=2
     jmax=2
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  dR = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           if (Symmetry==AXISYM) then
              dR = dR + dV  * f(i,j,k) * w(i,j,k) * X(i,1,1)
           else
              dR = dR + dV  * f(i,j,k) * w(i,j,k)
!              write(*,*) "i,j,k = ", i,j,k,f(i,j,k),w(i,j,k),imax,jmax,kmax,adjkmax
           end if
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  return
end subroutine wavg
