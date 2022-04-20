!----------------------------------------------------------------------------
!
! $Id: mass_integral.F90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!

!-----------------------------------------------------------------------------
!
! Integrate the rest mass integrand over the grid.  Since the functions are
! already centered, the integral is just a sum times dX dY dZ.
!
!-----------------------------------------------------------------------------
subroutine rest_mass_integral(ex, dmass, &
     X, Y, Z, &
     rho_b, &
     phi, alpha, u0, Symmetry, &
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
  integer                                     :: Symmetry
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
           dmass = dmass + u0(i,j,k) * &
                (alpha(i,j,k) + ONE) * exp(SIX * phi(i,j,k)) * rho_b(i,j,k)


        end do ! i-loop
     end do ! j-loop
  end do ! j-loop

  dmass = dV * dmass

  return

end subroutine rest_mass_integral

!-----------------------------------------------------------------------------
!
! Integrate the rest mass integrand over the grid.  Since the functions are
! already centered, the integral is just a sum times dX dY dZ.  This 
! function integrates over rho_star, unlike the previous function...
!
!-----------------------------------------------------------------------------
subroutine mass_0_integral(ex, dmass, &
     X, Y, Z, &
     rho_star, Symmetry, &
     adjimin, adjimax, &
     adjjmin, adjjmax, &
     adjkmin, adjkmax,proc)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_star
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax,proc
  integer                                     :: Symmetry
  integer                                     :: AXISYM
  parameter(AXISYM = 4)

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
  real*8                             :: ZERO, FOUR, ONE, SIX, PI
  parameter( ZERO = 0.D0, FOUR = 4.D0, ONE = 1.D0 )
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

!  imin = lbound(phi,1) - adjimin
!  jmin = lbound(phi,2) - adjjmin
!  kmin = lbound(phi,3) - adjkmin
!  imax = ubound(phi,1) - adjimax
!  jmax = ubound(phi,2) - adjjmax
!  kmax = ubound(phi,3) - adjkmax

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
! Set up integration
!

  dmass = ZERO
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax

           if (Symmetry == AXISYM) then
             dmass = dmass + X(i,1,1)*rho_star(i,j,k)
           else 
             dmass = dmass + rho_star(i,j,k)
           endif

        end do ! i-loop
     end do ! j-loop
  end do ! k-loop

  dmass = dV * dmass 

  return

end subroutine mass_0_integral


!-----------------------------------------------------------------------------
!
! Integrate the total proper volume of the grid.  
!
!-----------------------------------------------------------------------------
subroutine proper_volume(ex, prop_volume, &
     X, Y, Z, &
     rho_star, rho_b, &
     adjimin, adjimax, &
     adjjmin, adjjmax, &
     adjkmin, adjkmax)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_star, rho_b
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
!
! output:
! 
  real*8                                      :: prop_volume
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
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
  dV = dX * dY * dZ
!
! Set up integration
!

  prop_volume = ZERO
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax

           prop_volume = prop_volume + rho_star(i,j,k)/rho_b(i,j,k)

        end do ! i-loop
     end do ! j-loop
  end do ! j-loop

  prop_volume= dV * prop_volume

  return

end subroutine proper_volume

