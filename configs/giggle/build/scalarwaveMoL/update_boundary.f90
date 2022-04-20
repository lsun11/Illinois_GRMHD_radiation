!-----------------------------------------------------------------------------
! Update outer boundaries
!-----------------------------------------------------------------------------
subroutine scalarwaveMoL_update_boundary_lowlevel(ex, fake_ex,X, Y, Z, &
     dX,dY,dZ,dT, &
     phi_old,phi_new,phi_dot_old,phi_dot_new, phi_stagger_old,phi_stagger_new, &
     phi_dot_stagger_old,phi_dot_stagger_new, &
     Symmetry,scawave_have_bdry_min,scawave_have_bdry_max)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ex,fake_ex,scawave_have_bdry_min,scawave_have_bdry_max
  real*8                                   :: dX,dY,dZ,dT
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi_old,phi_new,X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi_dot_old,phi_dot_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi_stagger_old,phi_stagger_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi_dot_stagger_old,phi_dot_stagger_new
  integer                                  :: Symmetry
  !
  ! Other variables:
  !
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku,interp_order,ill,iuu
  integer                    :: OCTANT,EQUATORIAL, AXISYM, LINEAR, CUBIC
  real*8                     :: Courant,Delx,Dely,Delz,f_int,r
  real*8                     :: HALF, ONE, TWO, yl,zl
  real*8                     :: xp,phip,D,xll,xl,xu,xuu,pll,pl,pu,puu,dx3
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
  parameter(EQUATORIAL = 1,OCTANT = 2, AXISYM = 4)
  parameter(LINEAR = 0, CUBIC = 1)
  !
  ! interpolation order
  interp_order = CUBIC
  imin = ex(1)-fake_ex(1)+1
  jmin = ex(2)-fake_ex(2)+1
  kmin = ex(3)-fake_ex(3)+1
  imax = fake_ex(1)
  jmax = fake_ex(2)
  kmax = fake_ex(3)
  if(Symmetry.ne.AXISYM) then
     ! x = xmin
     !
     if(scawave_have_bdry_min(1)==1) then
        i= imin
        do k = kmin, kmax - 1
           do j = jmin+1, jmax - 1
              call scalarwave_weights(X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              call scalarwave_interp(ex,phi_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              r = sqrt(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k) + Z(i,j,k)*Z(i,j,k))
              phi_new(i,j,k) = f_int * (1.D0-dT/r)
              call scalarwave_interp(ex,phi_dot_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              phi_dot_new(i,j,k) = f_int * (1.D0-dT/r)
              !!xl = X(i,j,k)
              !!yl = Y(i,j,k) + 0.5d0*dY
              !!zl = Z(i,j,k) + 0.5d0*dZ
              !!call scalarwave_weights(xl,yl,zl,dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
              !!     kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!call scalarwave_interp(ex,phi_stagger_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!r = sqrt(xl*xl + yl*yl + zl*zl)
              !!phi_stagger_new(i,j,k) = f_int * (1.D0-dT/r)
              !!call scalarwave_interp(ex,phi_dot_stagger_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!phi_dot_stagger_new(i,j,k) = f_int * (1.D0-dT/r)
              ! *** TEST ***
              ! doing quadratic extrapolation:
              phi_stagger_new(i,j,k) = 3.d0*phi_stagger_new(i+1,j,k) - 3.d0*phi_stagger_new(i+2,j,k) + phi_stagger_new(i+3,j,k)
              phi_dot_stagger_new(i,j,k) = 3.d0*phi_dot_stagger_new(i+1,j,k) - 3.d0*phi_dot_stagger_new(i+2,j,k) + phi_dot_stagger_&
  &new(i+3,j,k)
              !!phi_stagger_new(i,j,k) = phi_stagger_old(i,j,k)
              !!phi_dot_stagger_new(i,j,k) = phi_dot_stagger_old(i,j,k)
              ! ************
           enddo
        enddo
     end if
     ! y = ymin
     !
     if(scawave_have_bdry_min(2)==1) then
        j = jmin
        do k = kmin, kmax - 1
           do i = imin, imax
              call scalarwave_weights(X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              call scalarwave_interp(ex,phi_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              r = sqrt(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k) + Z(i,j,k)*Z(i,j,k))
              phi_new(i,j,k) = f_int * (1.D0-dT/r)
              call scalarwave_interp(ex,phi_dot_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              phi_dot_new(i,j,k) = f_int * (1.D0-dT/r)
              !!xl = X(i,j,k)
              !!yl = Y(i,j,k) + 0.5d0*dY
              !!zl = Z(i,j,k) + 0.5d0*dZ
              !!call scalarwave_weights(xl,yl,zl,dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
              !!     kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!call scalarwave_interp(ex,phi_stagger_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!r = sqrt(xl*xl + yl*yl + zl*zl)
              !!phi_stagger_new(i,j,k) = f_int * (1.D0-dT/r)
              !!call scalarwave_interp(ex,phi_dot_stagger_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!phi_dot_stagger_new(i,j,k) = f_int * (1.D0-dT/r)
              ! *** TEST ***
              ! doing quadratic extrapolation:
              phi_stagger_new(i,j,k) = 3.d0*phi_stagger_new(i,j+1,k) - 3.d0*phi_stagger_new(i,j+2,k) + phi_stagger_new(i,j+3,k)
              phi_dot_stagger_new(i,j,k) = 3.d0*phi_dot_stagger_new(i,j+1,k) - 3.d0*phi_dot_stagger_new(i,j+2,k) + phi_dot_stagger_&
  &new(i,j+3,k)
              !!phi_stagger_new(i,j,k) = phi_stagger_old(i,j,k)
              !!phi_dot_stagger_new(i,j,k) = phi_dot_stagger_old(i,j,k)
              ! ************
           enddo
        enddo
     end if
     ! x = xmax
     !
     if(scawave_have_bdry_max(1)==1) then
        i= imax
        do k = kmin, kmax - 1
           do j = jmin, jmax - 1
              call scalarwave_weights(X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              call scalarwave_interp(ex,phi_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              r = sqrt(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k) + Z(i,j,k)*Z(i,j,k))
              phi_new(i,j,k) = f_int * (1.D0-dT/r)
              call scalarwave_interp(ex,phi_dot_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              phi_dot_new(i,j,k) = f_int * (1.D0-dT/r)
              !!xl = X(i,j,k)
              !!yl = Y(i,j,k) + 0.5d0*dY
              !!zl = Z(i,j,k) + 0.5d0*dZ
              !!call scalarwave_weights(xl,yl,zl,dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
              !!     kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!call scalarwave_interp(ex,phi_stagger_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!r = sqrt(xl*xl + yl*yl + zl*zl)
              !!phi_stagger_new(i,j,k) = f_int * (1.D0-dT/r)
              !!call scalarwave_interp(ex,phi_dot_stagger_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!phi_dot_stagger_new(i,j,k) = f_int * (1.D0-dT/r)
              ! *** TEST ***
              ! doing quadratic extrapolation:
              phi_stagger_new(i,j,k) = 3.d0*phi_stagger_new(i-1,j,k) - 3.d0*phi_stagger_new(i-2,j,k) + phi_stagger_new(i-3,j,k)
              phi_dot_stagger_new(i,j,k) = 3.d0*phi_dot_stagger_new(i-1,j,k) - 3.d0*phi_dot_stagger_new(i-2,j,k) + phi_dot_stagger_&
  &new(i-3,j,k)
              !!phi_stagger_new(i,j,k) = phi_stagger_old(i,j,k)
              !!phi_dot_stagger_new(i,j,k) = phi_dot_stagger_old(i,j,k)
              ! ************
          enddo
        enddo
     end if
     ! z = zmax
     !
     if(scawave_have_bdry_max(3)==1) then
        k = kmax
        do j = jmin, jmax
           do i = imin, imax
              call scalarwave_weights(X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              call scalarwave_interp(ex,phi_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              r = sqrt(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k) + Z(i,j,k)*Z(i,j,k))
              phi_new(i,j,k) = f_int * (1.D0-dT/r)
              call scalarwave_interp(ex,phi_dot_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              phi_dot_new(i,j,k) = f_int * (1.D0-dT/r)
              !!xl = X(i,j,k)
              !!yl = Y(i,j,k) + 0.5d0*dY
              !!zl = Z(i,j,k) + 0.5d0*dZ
              !!call scalarwave_weights(xl,yl,zl,dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
              !!     kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!call scalarwave_interp(ex,phi_stagger_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!r = sqrt(xl*xl + yl*yl + zl*zl)
              !!phi_stagger_new(i,j,k) = f_int * (1.D0-dT/r)
              !!call scalarwave_interp(ex,phi_dot_stagger_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!phi_dot_stagger_new(i,j,k) = f_int * (1.D0-dT/r)
              ! *** TEST ***
              ! doing quadratic extrapolation:
              phi_stagger_new(i,j,k) = 3.d0*phi_stagger_new(i,j,k-1) - 3.d0*phi_stagger_new(i,j,k-2) + phi_stagger_new(i,j,k-3)
              phi_dot_stagger_new(i,j,k) = 3.d0*phi_dot_stagger_new(i,j,k-1) - 3.d0*phi_dot_stagger_new(i,j,k-2) + phi_dot_stagger_&
  &new(i,j,k-3)
              !!phi_stagger_new(i,j,k) = phi_stagger_old(i,j,k)
              !!phi_dot_stagger_new(i,j,k) = phi_dot_stagger_old(i,j,k)
              ! ************
           enddo
        enddo
     end if
     ! y = ymax
     !
     if(scawave_have_bdry_max(2)==1) then
        j = jmax
        do k = kmin, kmax - 1
           do i = imin, imax
              call scalarwave_weights(X(i,1,1),Y(1,j,1),Z(1,1,k),dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              call scalarwave_interp(ex,phi_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              r = sqrt(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k) + Z(i,j,k)*Z(i,j,k))
              phi_new(i,j,k) = f_int * (1.D0-dT/r)
              call scalarwave_interp(ex,phi_dot_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              phi_dot_new(i,j,k) = f_int * (1.D0-dT/r)
              !!xl = X(i,j,k)
              !!yl = Y(i,j,k) + 0.5d0*dY
              !!zl = Z(i,j,k) + 0.5d0*dZ
              !!call scalarwave_weights(xl,yl,zl,dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
              !!     kmin,kmax,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!call scalarwave_interp(ex,phi_stagger_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!r = sqrt(xl*xl + yl*yl + zl*zl)
              !!phi_stagger_new(i,j,k) = f_int * (1.D0-dT/r)
              !!call scalarwave_interp(ex,phi_dot_stagger_old,f_int,il,iu,jl,ju,kl,ku,Delx,Dely,Delz)
              !!phi_dot_stagger_new(i,j,k) = f_int * (1.D0-dT/r)
              ! *** TEST ***
              ! doing quadratic extrapolation:
              phi_stagger_new(i,j,k) = 3.d0*phi_stagger_new(i,j-1,k) - 3.d0*phi_stagger_new(i,j-2,k) + phi_stagger_new(i,j-3,k)
              phi_dot_stagger_new(i,j,k) = 3.d0*phi_dot_stagger_new(i,j-1,k) - 3.d0*phi_dot_stagger_new(i,j-2,k) + phi_dot_stagger_&
  &new(i,j-3,k)
              !!phi_stagger_new(i,j,k) = phi_stagger_old(i,j,k)
              !!phi_dot_stagger_new(i,j,k) = phi_dot_stagger_old(i,j,k)
              ! ************
           enddo
        enddo
     end if
  end if
  if(Symmetry==AXISYM) then
     write(*,*) "Symmetry==4 (axisymmetry) is not yet supported with Carpet!"
     ! y = ymin, y = ymax
     !
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
           r = sqrt(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k) + Z(i,j,k)*Z(i,j,k))
           phi_new(i,j,k) = f_int * (1.D0-dT/r)
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
           r = sqrt(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k) + Z(i,j,k)*Z(i,j,k))
           phi_new(i,j,k) = f_int * (1.D0-dT/r)
        enddo
     enddo
  end if
end subroutine scalarwaveMoL_update_boundary_lowlevel
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
  ! get upper values
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
