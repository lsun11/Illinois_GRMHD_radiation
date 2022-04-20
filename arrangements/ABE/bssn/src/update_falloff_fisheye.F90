!-----------------------------------------------------------------------------
!
! $Id: update_falloff.f90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!
!-----------------------------------------------------------------------------
!
! Update outer boundaries with x^n1 y^n2 z^n3 / r^m falloff condition
!
!-----------------------------------------------------------------------------
subroutine update_falloff_fisheye(ex, X, Y, Z, Rp, f, f_old, &
     n1, n2, n3, m , Symmetry)
implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: f, f_old, Rp
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
  real*8                     :: xpu,ypu,zpu,rpu,rpl,xpl,ypl,zpl
  real*8                     :: HALF, ONE, TWO
!  real*8, dimension(ex(1),ex(2),ex(3))     :: f_old
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
  f_old = f
!
! x = max
!
  i = imax
  il = i - 1
  xu = X(imax,1,1)
  xl = X(imax - 1,1,1)
  do k = kmin,  kmax
     do j = jmin, jmax
	yu = Y(1,j,1)
	zu = Z(1,1,k)
        ru = sqrt(xu**2 + yu**2 + zu**2)
	rl = sqrt(xl**2 + yu**2 + zu**2)
	rpu = Rp(imax,j,k)
	rpl = Rp(imax-1,j,k)
        xpu = xu/ru * rpu
	xpl = xl/rl * rpl
        f(i,j,k) = f_old(il,j,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1)
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
     do k = kmin,  kmax
        do i = imin, imax
	   xu = X(i,1,1)
	   zu = Z(1,1,k)
           ru = sqrt(xu**2 + yu**2 + zu**2)
           rl = sqrt(xu**2 + yl**2 + zu**2)
	   rpu = Rp(i,jmax,k)
           rpl = Rp(i,jmax-1,k)
	   ypu = yu/ru * rpu
	   ypl = yl/rl * rpl
           f(i,j,k) = f_old(i,jl,k)*(rpl**m / rpu**m) * (ypu**n2 / ypl**n2) 
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
  do j = jmin, jmax - 1
     do i = imin, imax - 1
        xu = X(i,1,1)
        yu = Y(1,j,1)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xu**2 + yu**2 + zl**2)
	rpu = Rp(i,j,kmax)
  	rpl = Rp(i,j,kmax-1)
	zpu = zu/ru * rpu
	zpl = zl/rl * rpl
        f(i,j,k) = f_old(i,j,kl)*(rpl**m / rpu**m) * (zpu**n3 / zpl**n3)
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
	rpu = Rp(imax,jmax,k)
        rpl = Rp(imax-1,jmax-1,k)
	xpu = xu/ru * rpu
	xpl = xl/rl * rpl
	ypu = yu/ru * rpu
	ypl = yl/rl * rpl
        f(i,j,k) = f_old(il,jl,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1) * &
           (ypu**n2 / ypl**n2)
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
     do k = kmin, kmax
        do j = jmin, jmax 
           yu = Y(1,j,1)
           zu = Z(1,1,k)
           ru = sqrt(xu**2 + yu**2 + zu**2)
           rl = sqrt(xl**2 + yu**2 + zu**2) 
	   rpu = Rp(imin,j,k)
	   rpl = Rp(imin+1,j,k)
	   xpu = xu/ru * rpu
	   xpl = xl/rl * rpl
           f(i,j,k) = f_old(il,j,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1)
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
  do k = kmin, kmax
     do i = imin, imax
        xu = X(i,1,1)
        zu = Z(1,1,k)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xu**2 + yl**2 + zu**2) 
	rpu = Rp(i,jmin,k)
	rpl = Rp(i,jmin+1,k)
	ypu = yu/ru * rpu
	ypl = yl/rl * rpl
        f(i,j,k) = f_old(i,jl,k)*(rpl**m / rpu**m) * (ypu**n2 / ypl**n2) 
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
     rpu = Rp(imin,jmax,k)
     rpl = Rp(imin+1,jmax-1,k)
     xpu = xu/ru * rpu
     xpl = xl/rl * rpl
     ypu = yu/ru * rpu
     ypl = yl/rl * rpl
     f(i,j,k) = f_old(il,jl,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1) * (ypu**n2 / ypl**n2)
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
     rpu = Rp(imax,jmin,k)
     rpl = Rp(imax-1,jmin+1,k)
     xpu = xu/ru * rpu
     xpl = xl/rl * rpl
     ypu = yu/ru * rpu
     ypl = yl/rl * rpl
     f(i,j,k) = f_old(il,jl,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1) * (ypu**n2 / ypl**n2)
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
     rpu = Rp(imin,jmin,k)
     rpl = Rp(imin+1,jmin+1,k)
     xpu = xu/ru * rpu
     xpl = xl/rl * rpl
     ypu = yu/ru * rpu
     ypl = yl/rl * rpl
     f(i,j,k) = f_old(il,jl,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1) * (ypu**n2 / ypl**n2)
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
  do j = jmin + 1, jmax
     do i = imin + 1, imax
        xu = X(i,1,1)
        yu = Y(1,j,1)
        ru = sqrt(xu**2 + yu**2 + zu**2)
        rl = sqrt(xu**2 + yu**2 + zl**2)
	rpu = Rp(i,j,kmin)
	rpl = Rp(i,j,kmin+1)
	zpu = zu/ru * rpu
	zpl = zl/rl * rpl
        f(i,j,k) = f_old(i,j,kl)*(rpl**m / rpu**m) * (zpu**n3 / zpl**n3)
     enddo
  enddo
end if
!

end subroutine update_falloff_fisheye
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!
! Update outer boundaries with x^n1 y^n2 z^n3 / r^m falloff condition
!  with specified imin,imax,jmin,jmax,kmin and kmax
!
!-----------------------------------------------------------------------------
subroutine update_falloff_fisheye2(ex, X, Y, Z, Rp, f, &
     n1, n2, n3, m , Symmetry,imin,imax,jmin,jmax,kmin,kmax)
implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: f, Rp
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
  real*8                     :: xpu,ypu,zpu,rpu,rpl,xpl,ypl,zpl
  real*8                     :: HALF, ONE, TWO
  real*8, dimension(ex(1),ex(2),ex(3))     :: f_old
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2)
  parameter(PI_SYMM = 3, AXISYM = 4)
!
! where the action is...
!
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
!
  f_old = f
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
	rpu = Rp(imax,j,k)
	rpl = Rp(imax-1,j,k)
        xpu = xu/ru * rpu
	xpl = xl/rl * rpl
        f(i,j,k) = f_old(il,j,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1)
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
	   rpu = Rp(i,jmax,k)
           rpl = Rp(i,jmax-1,k)
	   ypu = yu/ru * rpu
	   ypl = yl/rl * rpl
           f(i,j,k) = f_old(i,jl,k)*(rpl**m / rpu**m) * (ypu**n2 / ypl**n2) 
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
	rpu = Rp(i,j,kmax)
  	rpl = Rp(i,j,kmax-1)
	zpu = zu/ru * rpu
	zpl = zl/rl * rpl
        f(i,j,k) = f_old(i,j,kl)*(rpl**m / rpu**m) * (zpu**n3 / zpl**n3)
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
	rpu = Rp(imax,jmax,k)
        rpl = Rp(imax-1,jmax-1,k)
	xpu = xu/ru * rpu
	xpl = xl/rl * rpl
	ypu = yu/ru * rpu
	ypl = yl/rl * rpl
        f(i,j,k) = f_old(il,jl,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1) * &
           (ypu**n2 / ypl**n2)
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
	   rpu = Rp(imin,j,k)
	   rpl = Rp(imin+1,j,k)
	   xpu = xu/ru * rpu
	   xpl = xl/rl * rpl
           f(i,j,k) = f_old(il,j,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1)
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
	rpu = Rp(i,jmin,k)
	rpl = Rp(i,jmin+1,k)
	ypu = yu/ru * rpu
	ypl = yl/rl * rpl
        f(i,j,k) = f_old(i,jl,k)*(rpl**m / rpu**m) * (ypu**n2 / ypl**n2) 
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
     rpu = Rp(imin,jmax,k)
     rpl = Rp(imin+1,jmax-1,k)
     xpu = xu/ru * rpu
     xpl = xl/rl * rpl
     ypu = yu/ru * rpu
     ypl = yl/rl * rpl
     f(i,j,k) = f_old(il,jl,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1) * (ypu**n2 / ypl**n2)
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
     rpu = Rp(imax,jmin,k)
     rpl = Rp(imax-1,jmin+1,k)
     xpu = xu/ru * rpu
     xpl = xl/rl * rpl
     ypu = yu/ru * rpu
     ypl = yl/rl * rpl
     f(i,j,k) = f_old(il,jl,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1) * (ypu**n2 / ypl**n2)
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
     rpu = Rp(imin,jmin,k)
     rpl = Rp(imin+1,jmin+1,k)
     xpu = xu/ru * rpu
     xpl = xl/rl * rpl
     ypu = yu/ru * rpu
     ypl = yl/rl * rpl
     f(i,j,k) = f_old(il,jl,k)*(rpl**m / rpu**m) * (xpu**n1 / xpl**n1) * (ypu**n2 / ypl**n2)
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
	rpu = Rp(i,j,kmin)
	rpl = Rp(i,j,kmin+1)
	zpu = zu/ru * rpu
	zpl = zl/rl * rpl
        f(i,j,k) = f_old(i,j,kl)*(rpl**m / rpu**m) * (zpu**n3 / zpl**n3)
     enddo
  enddo
end if
!

end subroutine update_falloff_fisheye2
!-----------------------------------------------------------------------------
