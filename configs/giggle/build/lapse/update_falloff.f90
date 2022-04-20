!-----------------------------------------------------------------------------
!
! $Id: update_falloff.f90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!
!-----------------------------------------------------------------------------
!
! Update outer boundaries with x^n1 y^n2 z^n3 / r^m falloff condition
!
!-----------------------------------------------------------------------------
subroutine update_falloff(ex, X, Y, Z, f, &
     n1, n2, n3, m , Symmetry)
implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: f
  integer                                  :: n1,n2,n3,m
  integer                                  :: Symmetry
!
! Other variables:
!
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku
  integer                    :: NO_SYMM, EQUATORIAL, OCTANT
  integer                    :: PI_SYMM, AXISYM
  real*8                     :: Courant,dX,dY,dZ
  real*8                     :: xu,yu,zu,ru,rl,xl,yl,zl
  real*8                     :: HALF, ONE, TWO
  real*8, dimension(ex(1),ex(2),ex(3))     :: fold
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2)
  parameter(PI_SYMM = 3, AXISYM = 4)
!
! where the action is...
!
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  jmax = ubound(f,2)
  kmax = ubound(f,3)
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
!
  fold = f
!
! x = max
!
  i = imax
  il = i - 1
  xu = X(imax,1,1)
  xl = X(imax - 1,1,1)
  do j = jmin, jmax
     do k = kmin,  kmax
        yu = Y(1,j,1)
        zu = Z(1,1,k)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xl**2 + yu**2 + zu**2)
        f(i,j,k) = fold(il,j,k)*(rl**m / ru**m) * (xu**n1 / xl**n1)
     enddo
  enddo
!
! y = max
!
  if(Symmetry.ne.AXISYM) then
     j = jmax
     jl = j - 1
     yu = Y(1,j,1)
     yl = Y(1,jl,1)
     do i = imin, imax
        do k = kmin,  kmax
           xu = X(i,1,1)
           zu = Z(1,1,k)
           ru = sqrt(xu**2 + yu**2 + zu**2)
           rl = sqrt(xu**2 + yl**2 + zu**2)
           f(i,j,k) = fold(i,jl,k)*(rl**m / ru**m) * (yu**n2 / yl**n2)
        enddo
     enddo
  end if
!
! z = max
!
  k = kmax
  kl = k - 1
  zu = Z(1,1,k)
  zl = Z(1,1,kl)
  do i = imin, imax - 1
     do j = jmin, jmax - 1
        xu = X(i,1,1)
        yu = Y(1,j,1)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xu**2 + yu**2 + zl**2)
        f(i,j,k) = fold(i,j,kl)*(rl**m / ru**m) * (zu**n3 / zl**n3)
     enddo
  enddo
!
! x = max, y = ymax
!
  if(Symmetry.ne.AXISYM) then
     i = imax
     il = i - 1
     xu = X(imax,1,1)
     xl = X(imax - 1,1,1)
     j = jmax
     jl = j - 1
     yu = Y(1,j,1)
     yl = Y(1,jl,1)
     do k = kmin,  kmax
        zu = Z(1,1,k)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xl**2 + yl**2 + zu**2)
        f(i,j,k) = fold(il,jl,k)*(rl**m / ru**m) * (xu**n1 / xl**n1) * (yu**n2 / yl**n2)
     enddo
  end if
  if (Symmetry == EQUATORIAL.or.Symmetry == NO_SYMM.or.Symmetry == PI_SYMM) then
!
! x = min
!
     i = imin
     il = imin + 1
     xu = X(imin,1,1)
     xl = X(il,1,1)
     do j = jmin, jmax
        do k = kmin, kmax
           yu = Y(1,j,1)
           zu = Z(1,1,k)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xl**2 + yu**2 + zu**2)
        f(i,j,k) = fold(il,j,k)*(rl**m / ru**m) * (xu**n1 / xl**n1)
     enddo
  enddo
  endif
  if (Symmetry == EQUATORIAL.or.Symmetry == NO_SYMM) then
!
! y = min
!
  j = jmin
  jl = jmin + 1
  yu = Y(1,jmin,1)
  yl = Y(1,jl,1)
  do i = imin, imax
     do k = kmin, kmax
        xu = X(i,1,1)
        zu = Z(1,1,k)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xu**2 + yl**2 + zu**2)
        f(i,j,k) = fold(i,jl,k)*(rl**m / ru**m) * (yu**n2 / yl**n2)
     enddo
  enddo
  endif
  if (Symmetry == EQUATORIAL.or.Symmetry == NO_SYMM.or.Symmetry == PI_SYMM) then
!
! x = min, y = ymax
!
  i = imin
  il = i + 1
  xu = X(i,1,1)
  xl = X(il,1,1)
  j = jmax
  jl = j - 1
  yu = Y(1,j,1)
  yl = Y(1,jl,1)
  do k = kmin,  kmax
     zu = Z(1,1,k)
     ru = sqrt(xu**2 + yu**2 + zu**2)
     rl = sqrt(xl**2 + yl**2 + zu**2)
     f(i,j,k) = fold(il,jl,k)*(rl**m / ru**m) * (xu**n1 / xl**n1) * (yu**n2 / yl**n2)
  enddo
  endif
  if (Symmetry == EQUATORIAL.or.Symmetry == NO_SYMM) then
!
! x = max, y = ymin
!
  i = imax
  il = i - 1
  xu = X(i,1,1)
  xl = X(il,1,1)
  j = jmin
  jl = j + 1
  yu = Y(1,j,1)
  yl = Y(1,jl,1)
  do k = kmin,  kmax
     zu = Z(1,1,k)
     ru = sqrt(xu**2 + yu**2 + zu**2)
     rl = sqrt(xl**2 + yl**2 + zu**2)
     f(i,j,k) = fold(il,jl,k)*(rl**m / ru**m) * (xu**n1 / xl**n1) * (yu**n2 / yl**n2)
  enddo
!
! x = min, y = ymin
!
  i = imin
  il = i + 1
  xu = X(i,1,1)
  xl = X(il,1,1)
  j = jmin
  jl = j + 1
  yu = Y(1,j,1)
  yl = Y(1,jl,1)
  do k = kmin,  kmax
     zu = Z(1,1,k)
     ru = sqrt(xu**2 + yu**2 + zu**2)
     rl = sqrt(xl**2 + yl**2 + zu**2)
     f(i,j,k) = fold(il,jl,k)*(rl**m / ru**m) * (xu**n1 / xl**n1) * (yu**n2 / yl**n2)
  enddo
!
end if
  if (Symmetry == NO_SYMM) then
!
! z = min
!
  k = kmin
  kl = k + 1
  zu = Z(1,1,k)
  zl = Z(1,1,kl)
  do i = imin + 1, imax
     do j = jmin + 1, jmax
        xu = X(i,1,1)
        yu = Y(1,j,1)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xu**2 + yu**2 + zl**2)
        f(i,j,k) = fold(i,j,kl)*(rl**m / ru**m) * (zu**n3 / zl**n3)
     enddo
  enddo
end if
!
end subroutine update_falloff
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Fix first derivatives using assumed falloff behavior
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine dfix(ex,f,fx,fy,fz,X,Y,Z,dX,dY,dZ,mx,my,mz,mr,Symmetry)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz,X,Y,Z
  real*8                                   :: dX,dY,dZ
  integer                                  :: mx,my,mz,mr,Symmetry
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm1,dym1,dzm1
  real*8               :: xl,xu,yl,yu,zl,zu,rl,ru
  real*8, parameter    :: HALF = 0.5D0, TWO = 2.D0
  integer              :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3)
!
! Input translation
!
  imin = lbound(f,1)
  imax = ubound(f,1)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
  dxm1 = HALF/dX
  dym1 = HALF/dY
  dzm1 = HALF/dZ
!
  i = imax
  xu = X(imax,1,1)+dX
  xl = X(imax,1,1)
  do j = jmin, jmax
     do k = kmin,  kmax
        yu = Y(1,j,1)
        zu = Z(1,1,k)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xl**2 + yu**2 + zu**2)
        fx(i,j,k) = dxm1 * (f(i,j,k)*((rl/ru)**mr)*((xu/xl)**mx) - f(i-1,j,k))
     enddo
  enddo
!
  j = jmax
  yu = Y(1,jmax,1)+dY
  yl = Y(1,jmax,1)
  do i = imin, imax
     do k = kmin,  kmax
        xu = X(i,1,1)
        zu = Z(1,1,k)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xu**2 + yl**2 + zu**2)
        fy(i,j,k) = dym1 * (f(i,j,k)*((rl/ru)**mr)*((yu/yl)**my) - f(i,j-1,k))
     enddo
  enddo
!
  k = kmax
  zu = Z(1,1,kmax)+dZ
  zl = Z(1,1,kmax)
  do i = imin, imax
     do j = jmin,  jmax
        xu = X(i,1,1)
        yu = Y(1,j,1)
!        yu = Y(k)
        ru = sqrt(xu**2 + yu**2 + zl**2)
        rl = sqrt(xu**2 + yu**2 + zu**2)
        fz(i,j,k) = dzm1 * (f(i,j,k)*((rl/ru)**mr)*((zu/zl)**mz) - f(i,j,k-1))
     enddo
  enddo
!
  if (Symmetry == EQUATORIAL.or.Symmetry == NO_SYMM.or.Symmetry == PI_SYMM) then
!
     i = imin
     xu = X(imin,1,1)-dX
     xl = X(imin,1,1)
     do j = jmin, jmax
        do k = kmin,  kmax
           yu = Y(1,j,1)
           zu = Z(1,1,k)
           ru = sqrt(xu**2 + yu**2 + zu**2)
           rl = sqrt(xl**2 + yu**2 + zu**2)
           fx(i,j,k) = -dxm1 * (f(i,j,k)*((rl/ru)**mr)*((xu/xl)**mx) - f(i+1,j,k))
        enddo
     enddo
  end if
!
  if (Symmetry == EQUATORIAL.or.Symmetry == NO_SYMM) then
!
     j = jmin
     yu = Y(1,jmin,1)-dY
     yl = Y(1,jmin,1)
     do i = imin, imax
        do k = kmin,  kmax
           xu = X(i,1,1)
           zu = Z(1,1,k)
           ru = sqrt(xu**2 + yu**2 + zu**2)
           rl = sqrt(xu**2 + yl**2 + zu**2)
           fy(i,j,k) = -dym1 * (f(i,j,k)*((rl/ru)**mr)*((yu/yl)**my) - f(i,j+1,k))
        enddo
     enddo
!
  end if
!
  if (Symmetry == NO_SYMM) then
!
     k = kmin
     zu = Z(1,1,kmin)
     zl = Z(1,1,kmin + 1)
     do i = imin, imax
        do j = jmin,  jmax
           xu = X(i,1,1)
           yu = Y(1,j,1)
           !           yu = Y(k)
           ru = sqrt(xu**2 + yu**2 + zu**2)
           rl = sqrt(xu**2 + yu**2 + zl**2)
           fz(i,j,k) = -dzm1 * (f(i,j,k)*((rl/ru)**mr)*((zu/zl)**mz) - f(i,j,k+1))
        enddo
     enddo
  end if
end subroutine dfix
