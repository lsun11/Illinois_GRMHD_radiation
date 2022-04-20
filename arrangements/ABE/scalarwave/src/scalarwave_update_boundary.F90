!-----------------------------------------------------------------------------
!
! $Id: update_boundary.F90,v 1.2 2006/03/10 22:32:24 zetienne Exp $
!
!-----------------------------------------------------------------------------
!
! Update outer boundaries
!
!-----------------------------------------------------------------------------
subroutine scalarwave_update_boundary(ex, X, Y, Z, dT, &
     phi_old,phi_new,phi_dot_old,phi_dot_new,Symmetry)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8                                   :: dT
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi_old,phi_new,X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi_dot_old,phi_dot_new
  integer                                  :: Symmetry
!
! Other variables:
!
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku,interp_order,ill,iuu
  integer                    :: OCTANT, AXISYM, LINEAR, CUBIC
  real*8                     :: Courant,dX,dY,dZ,Delx,Dely,Delz,f_int
  real*8                     :: HALF, ONE, TWO
  real*8                     :: xp,phip,D,xll,xl,xu,xuu,pll,pl,pu,puu,dx3
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
  parameter(OCTANT = 0, AXISYM = 1)
  parameter(LINEAR = 0, CUBIC = 1)
!
! interpolation order
  interp_order = CUBIC

! Next, properly set the min. array values, depending on whether we are near an 
! axis.  This properly accounts for the ghost zones.
!  if(x(lbound(phi_old,1),1,1) < 0) then
!     imin = lbound(phi_old,1)+1
!  else 
!     imin = lbound(phi_old,1)
!  endif
!  if(y(1,lbound(phi_old,2),1) < 0) then
!     jmin = lbound(phi_old,2)+1
!  else 
!     jmin = lbound(phi_old,2)
!  endif
!  if(z(1,1,lbound(phi_old,3)) < 0) then
!     kmin = lbound(phi_old,3)+1
!  else 
!     kmin = lbound(phi_old,3)
!  endif

  imin = lbound(phi_old,1)
  jmin = lbound(phi_old,2)
  kmin = lbound(phi_old,3)

  imax = ubound(phi_old,1)
  jmax = ubound(phi_old,2)
  kmax = ubound(phi_old,3)
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
  !write (*,*) X(imin), Y(jmin), Z(kmin)
  !  write (*,*) imin,jmin,kmin,imax,jmax,kmax
  !write (*,*) Symmetry

  if(Symmetry==OCTANT) then
  !write (*,*) "octant"
  !
  ! x = max
  !
     i= imax
     do k = kmin, kmax - 1
        do j = jmin, jmax - 1
           call scalarwave_weights(X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
           call scalarwave_interp(ex,phi_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
           phi_new(i,j,k) = f_int

!           write (*,*) X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz
        enddo
     enddo
  !
  ! z = zmax
  !
     k = kmax
     do j = jmin, jmax
        do i = imin, imax
           call scalarwave_weights(X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
           call scalarwave_interp(ex,phi_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
           phi_new(i,j,k) = f_int
        enddo
     enddo
  !
  ! y = ymax
  !
     j = jmax
     do k = kmin, kmax - 1
        do i = imin, imax
           call scalarwave_weights(X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
           call scalarwave_interp(ex,phi_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
           phi_new(i,j,k) = f_int
        enddo
     enddo

  else if(Symmetry==AXISYM) then
  !
  ! y = ymin, y = ymax
  !
  !   write (*,*) "axisymmetry"
   imin = lbound(phi_old,1)
   jmin = lbound(phi_old,2)
   kmin = lbound(phi_old,3)
     if(interp_order==LINEAR) then
        do k = kmin, kmax
           do i = imin, imax-1
              xp = sqrt(X(i,1,1)*X(i,1,1) + dY*dY)
              D  = (xp - X(i,1,1))/dX
              phip = (1.d0-D)*phi_new(i,2,k) + D*phi_new(i+1,2,k)
              phi_new(i,1,k) = phip
              phi_new(i,3,k) = phip
           enddo
        enddo
     else if(interp_order==CUBIC) then
        dx3 = dX * dX * dX
        do k = kmin, kmax
           do i = imin, imax-1
              xp = sqrt(X(i,1,1)*X(i,1,1) + dY*dY)
              il = i
!              if(il.lt.imin+1) il = imin+1
!              if(il.gt.imax-2) il = imax-2
              if(il.lt.imin+2) il = imin+2
              if(il.gt.imax-2) il = imax-2
              ill = il-1
              iu  = il+1
              iuu = il+2
              xll = X(ill,1,1)
              xl  = xll + dX
              xu  = xl  + dX
              xuu = xu  + dX
              pll = phi_new(ill,2,k)
              pl  = phi_new(il,2,k)
              pu  = phi_new(iu,2,k)
              puu = phi_new(iuu,2,k)
              phip = (xp-xl)*(xp-xu)*(xp-xuu) *pll / (-6.d0 * dx3) &
                   + (xp-xll)*(xp-xu)*(xp-xuu)*pl  / ( 2.d0 * dx3) &
                   + (xp-xll)*(xp-xl)*(xp-xuu)*pu  / (-2.d0 * dx3) &
                   + (xp-xll)*(xp-xl)*(xp-xu) *puu / ( 6.d0 * dx3)
              phi_new(i,1,k) = phip
              phi_new(i,3,k) = phip
           enddo
        enddo
     end if
  !
  ! x = max
  !
     i = imax
     do j = jmin, jmax - 1
        do k = kmin, kmax - 1
           call scalarwave_weights(X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
           call scalarwave_interp2d(ex,phi_old,f_int,il,iu,j,kl,ku,Delx,Delz)
           phi_new(i,j,k) = f_int
        enddo
     enddo
  !
  ! z = zmax
  !
     k = kmax
     do i = imin, imax
        do j = jmin, jmax
           call scalarwave_weights(X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
           call scalarwave_interp2d(ex,phi_old,f_int,il,iu,j,kl,ku,Delx,Delz)
           phi_new(i,j,k) = f_int
        enddo
     enddo
  end if

end subroutine scalarwave_update_boundary
!
!-----------------------------------------------------------------------------
!
! Interpolation for point i,j,k
!
!-----------------------------------------------------------------------------
!
subroutine scalarwave_weights(X,Y,Z,dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax,kmin,kmax, &
     il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
  implicit none
  real*8, intent(in)   :: x,y,z,dX,dY,dZ,dT
  integer, intent(in)  :: i,j,k
  integer, intent(in)  :: imin,imax,jmin,jmax,kmin,kmax
  integer, intent(out) :: il,iu,jl,ju,kl,ku
  real*8, intent(out)  :: Delx,Dely,Delz
!
  real*8               :: xo,yo,zo,xl,xu,yl,yu,zl,zu,r,fac,ONE
  parameter ( ONE = 1.D0 )
!
! find intersection of characteristic with previous timeslice
! 
  r = sqrt(X**2 + Y**2 + Z**2)
  fac = ONE - dT/r
  xo = X*fac
  yo = Y*fac
  zo = Z*fac
!
! find indices for interpolation
!        
  if (xo < X) then
     il = i - 1
     xl = X - dX
  else
     il = i
     xl = X
  end if
  if (yo < Y) then
     jl = j - 1
     yl = Y - dY
  else
     jl = j
     yl = Y
  end if
  if (zo < Z) then
     kl = k - 1
     zl = Z - dZ
  else
     kl = k
     zl = Z
  end if
!
! make sure that interpolation molecule does not extent beyond grid:
!
  if (il < imin) then 
     il = imin
     xl = X
  end if
  if (jl < jmin) then 
     jl = jmin
     yl = Y
  end if
  if (kl < kmin) then
     kl = kmin
     zl = Z
  end if
  if (il > imax - 1) then 
     il = imax - 1
     xl = X - dX
  end if
  if (jl > jmax - 1) then 
     jl = jmax - 1
     yl = Y - dY
  end if
  if (kl > kmax - 1) then
     kl = kmax - 1
     zl = Z - dZ
  end if
!
! get 'upper' values
!
  iu = il + 1
  ju = jl + 1
  ku = kl + 1
  xu = xl + dX
  yu = yl + dY
  zu = zl + dZ
!
! find coefficients for interpolation
!
  Delx = (xo - xl)/(xu - xl)
  Dely = (yo - yl)/(yu - yl)
  Delz = (zo - zl)/(zu - zl)

  if ((xu - xl).eq.0 .or. (yu - yl).eq.0 .or. (zu - zl).eq.0) then
     write (*,*) "Uh oh"
  endif
  return
end subroutine scalarwave_weights
!
!-----------------------------------------------------------------------------
!
! Interpolate function f using weights Delx, Dely and Delz
!
!-----------------------------------------------------------------------------
!
subroutine scalarwave_interp(ex,f,f_int,il,iu,jl,ju,kl,ku,Dx,Dy,Dz)
  implicit none
!
! Input:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f
  real*8                                   :: f_int
  integer                                  :: il,iu,jl,ju,kl,ku
  real*8                                   :: Dx,Dy,Dz
!
! Other variables:
!
  real*8      :: l
  parameter ( l = 1.D0 )
!
!
  f_int = (l - Dx)*( (l - Dy)*( (l - Dz)*f(il,jl,kl) + Dz*f(il,jl,ku) )   &
                        + Dy *( (l - Dz)*f(il,ju,kl) + Dz*f(il,ju,ku) ) ) &
             + Dx *( (l - Dy)*( (l - Dz)*f(iu,jl,kl) + Dz*f(iu,jl,ku) )   &
                        + Dy *( (l - Dz)*f(iu,ju,kl) + Dz*f(iu,ju,ku) ) )
  return
end subroutine scalarwave_interp
!
!-----------------------------------------------------------------------------
!
! 2D Interpolate function f using weights Delx and Delz
!
!-----------------------------------------------------------------------------
!
subroutine scalarwave_interp2d(ex,f,f_int,il,iu,j,kl,ku,Dx,Dz)
  implicit none
!
! Input:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f
  real*8                                   :: f_int
  integer                                  :: il,iu,j,kl,ku
  real*8                                   :: Dx,Dz
!
! Other variables:
!
  real*8      :: l
  parameter ( l = 1.D0 )
!
!
  f_int = (l - Dx)*( (l - Dz)*f(il,j,kl) + Dz*f(il,j,ku) ) &
             + Dx *( (l - Dz)*f(iu,j,kl) + Dz*f(iu,j,ku) )
  return
end subroutine scalarwave_interp2d


