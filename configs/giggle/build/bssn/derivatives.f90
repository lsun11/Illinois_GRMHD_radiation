!-----------------------------------------------------------------------------
!
! $Id: new_derivatives.f90
!
!-----------------------------------------------------------------------------
!
!  This file contains all kinds of derivative operators, all of them
!  *global* derivatives!
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! First derivatives
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gderivs(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm1,dym1,dzm1
  real*8, parameter    :: HALF = 0.5D0, TWO = 2.D0
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
! Derivative with respect to x
!
! TEMPORARY: THIS LINE INVALIDATES THE USE OF octant SYMMETRY
!  fx(imin,:,:) = TWO * dxm1 * (f(imin+1,:,:) - f(imin,:,:))
  fx(imin,:,:) = dxm1 * (f(imin+1,:,:) - SYM1*f(imin,:,:))
  do i = imin + 1, imax - 1
     fx(i,:,:) = dxm1 * (f(i+1,:,:) - f(i-1,:,:))
  enddo
  fx(imax,:,:) = TWO * dxm1 * (f(imax,:,:) - f(imax-1,:,:))
!
! Derivative with respect to y
!
!
! TEMPORARY: THIS LINE INVALIDATES THE USE OF octant SYMMETRY
!  fy(:,jmin,:) = TWO * dym1 * (f(:,jmin+1,:) - f(:,jmin,:))
  fy(:,jmin,:) = dym1 * (f(:,jmin+1,:) - SYM2*f(:,jmin,:))
  do j = jmin + 1, jmax - 1
     fy(:,j,:) = dym1 * (f(:,j+1,:) - f(:,j-1,:))
  enddo
  fy(:,jmax,:) = TWO * dym1 * (f(:,jmax,:) - f(:,jmax-1,:))
!
! Derivative with respect to z
!
  fz(:,:,kmin) = dzm1 * (f(:,:,kmin+1) - SYM3*f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     fz(:,:,k) = dzm1 * (f(:,:,k+1) - f(:,:,k-1))
  enddo
  fz(:,:,kmax) = TWO * dzm1 * (f(:,:,kmax) - f(:,:,kmax-1))
end subroutine gderivs
subroutine gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm1,dym1,dzm1
  real*8, parameter    :: HALF = 0.5D0, TWO = 2.D0
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
! Derivative with respect to x
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  fx(imin,:,:) = dxm1 * (f(imin+1,:,:) - SYM1*f(imin,:,:))
  do i = imin + 1, imax - 1
     fx(i,:,:) = dxm1 * (f(i+1,:,:) - f(i-1,:,:))
  enddo
  fx(imax,:,:) = TWO * dxm1 * (f(imax,:,:) - f(imax-1,:,:))
!
! Derivative with respect to y
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  fy(:,jmin,:) = dym1 * (f(:,jmin+1,:) - SYM2*f(:,jmin,:))
  do j = jmin + 1, jmax - 1
     fy(:,j,:) = dym1 * (f(:,j+1,:) - f(:,j-1,:))
  enddo
  fy(:,jmax,:) = TWO * dym1 * (f(:,jmax,:) - f(:,jmax-1,:))
!
! Derivative with respect to z
!
  fz(:,:,kmin) = dzm1 * (f(:,:,kmin+1) - SYM3*f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     fz(:,:,k) = dzm1 * (f(:,:,k+1) - f(:,:,k-1))
  enddo
  fz(:,:,kmax) = TWO * dzm1 * (f(:,:,kmax) - f(:,:,kmax-1))
end subroutine gderivs_oct
subroutine gderivs_oct_fix(ex,X,Y,Z,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm1,dym1,dzm1
  real*8, parameter    :: HALF = 0.5D0, TWO = 2.D0
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
  do i=imin,imax
     if(X(i,1,1) .lt. 0.D0) then
        imin = imin + 1
     end if
  end do
  do j=jmin,jmax
     if(Y(1,j,1) .lt. 0.D0) then
        jmin = jmin + 1
     end if
  end do
  do k=kmin,kmax
     if(Z(1,1,k) .lt. 0.D0) then
        kmin = kmin + 1
     end if
  end do
!
! Derivative with respect to x
!
! TEMPORARY: THIS LINE INVALIDATES THE USE OF octant SYMMETRY
!  fx(imin,:,:) = TWO * dxm1 * (f(imin+1,:,:) - f(imin,:,:))
  fx(imin,:,:) = dxm1 * (f(imin+1,:,:) - SYM1*f(imin,:,:))
  do i = imin + 1, imax - 1
     fx(i,:,:) = dxm1 * (f(i+1,:,:) - f(i-1,:,:))
  enddo
  fx(imax,:,:) = TWO * dxm1 * (f(imax,:,:) - f(imax-1,:,:))
!
! Derivative with respect to y
!
!
! TEMPORARY: THIS LINE INVALIDATES THE USE OF octant SYMMETRY
!  fy(:,jmin,:) = TWO * dym1 * (f(:,jmin+1,:) - f(:,jmin,:))
  fy(:,jmin,:) = dym1 * (f(:,jmin+1,:) - SYM2*f(:,jmin,:))
  do j = jmin + 1, jmax - 1
     fy(:,j,:) = dym1 * (f(:,j+1,:) - f(:,j-1,:))
  enddo
  fy(:,jmax,:) = TWO * dym1 * (f(:,jmax,:) - f(:,jmax-1,:))
!
! Derivative with respect to z
!
  fz(:,:,kmin) = dzm1 * (f(:,:,kmin+1) - SYM3*f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     fz(:,:,k) = dzm1 * (f(:,:,k+1) - f(:,:,k-1))
  enddo
  fz(:,:,kmax) = TWO * dzm1 * (f(:,:,kmax) - f(:,:,kmax-1))
end subroutine gderivs_oct_fix
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! First derivatives in equatorial symmetry
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gderivs_eq(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm1,dym1,dzm1
  real*8, parameter    :: HALF = 0.5D0, TWO = 2.D0
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
! Derivative with respect to x
!
  fx(imin,:,:) = TWO * dxm1 * (f(imin+1,:,:) - f(imin,:,:))
  do i = imin + 1, imax - 1
     fx(i,:,:) = dxm1 * (f(i+1,:,:) - f(i-1,:,:))
  enddo
  fx(imax,:,:) = TWO * dxm1 * (f(imax,:,:) - f(imax-1,:,:))
!  write(*,*)derivs_eq:,f(1,2,1),f(3,2,1),fx(2,2,1),dxm1
!
! Derivative with respect to y
!
  fy(:,jmin,:) = TWO * dym1 * (f(:,jmin+1,:) - f(:,jmin,:))
  do j = jmin + 1, jmax - 1
     fy(:,j,:) = dym1 * (f(:,j+1,:) - f(:,j-1,:))
  enddo
  fy(:,jmax,:) = TWO * dym1 * (f(:,jmax,:) - f(:,jmax-1,:))
!
! Derivative with respect to z
!
  fz(:,:,kmin) = dzm1 * (f(:,:,kmin+1) - SYM3*f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     fz(:,:,k) = dzm1 * (f(:,:,k+1) - f(:,:,k-1))
  enddo
  fz(:,:,kmax) = TWO * dzm1 * (f(:,:,kmax) - f(:,:,kmax-1))
end subroutine gderivs_eq
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! First derivatives in Pi symmetry
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gderivs_pi(ex,f,fx,fy,fz,dX,dY,dZ,SYM3,  &
                        glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
                        PI_SYMM_f)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8                                   :: dX,dY,dZ,SYM3
  integer                                  :: Nx,Nz
  integer                                  :: glob_imin,glob_jmin,glob_kmin
  real*8, dimension(Nx+1,Nz+1)                   :: PI_SYMM_f
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm1,dym1,dzm1
  real*8, parameter    :: HALF = 0.5D0, TWO = 2.D0
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
! Derivative with respect to x
!
  fx(imin,:,:) = TWO * dxm1 * (f(imin+1,:,:) - f(imin,:,:))
  do i = imin + 1, imax - 1
     fx(i,:,:) = dxm1 * (f(i+1,:,:) - f(i-1,:,:))
  enddo
  fx(imax,:,:) = TWO * dxm1 * (f(imax,:,:) - f(imax-1,:,:))
!
! Derivative with respect to y
!
  do k = kmin, kmax
    do i = imin, imax
      fy(i,jmin,k) = dym1 * (f(i,jmin+1,k) - PI_SYMM_f(i+glob_imin,k+glob_kmin))
    end do
  end do
  do j = jmin + 1, jmax - 1
     fy(:,j,:) = dym1 * (f(:,j+1,:) - f(:,j-1,:))
  enddo
  fy(:,jmax,:) = TWO * dym1 * (f(:,jmax,:) - f(:,jmax-1,:))
!
! Derivative with respect to z
!
  fz(:,:,kmin) = dzm1 * (f(:,:,kmin+1) - SYM3*f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     fz(:,:,k) = dzm1 * (f(:,:,k+1) - f(:,:,k-1))
  enddo
  fz(:,:,kmax) = TWO * dzm1 * (f(:,:,kmax) - f(:,:,kmax-1))
end subroutine gderivs_pi
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! First derivatives in Axisymmetry
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gderivs_axi(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm1,dym1,dzm1
  real*8, parameter    :: HALF = 0.5D0, TWO = 2.D0
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
! Derivative with respect to x
!
  fx(imin,:,:) = dxm1 * (f(imin+1,:,:) - SYM1*f(imin,:,:))
  do i = imin + 1, imax - 1
     fx(i,:,:) = dxm1 * (f(i+1,:,:) - f(i-1,:,:))
  enddo
  fx(imax,:,:) = TWO * dxm1 * (f(imax,:,:) - f(imax-1,:,:))
!
! Derivative with respect to y
!
  fy(:,jmin,:) = dym1 * (f(:,jmin+1,:) - SYM2*f(:,jmin,:))
  do j = jmin + 1, jmax - 1
     fy(:,j,:) = dym1 * (f(:,j+1,:) - f(:,j-1,:))
  enddo
  fy(:,jmax,:) = TWO * dym1 * (f(:,jmax,:) - f(:,jmax-1,:))
!
! Derivative with respect to z
!
  fz(:,:,kmin) = dzm1 * (f(:,:,kmin+1) - SYM3*f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     fz(:,:,k) = dzm1 * (f(:,:,k+1) - f(:,:,k-1))
  enddo
  fz(:,:,kmax) = TWO * dzm1 * (f(:,:,kmax) - f(:,:,kmax-1))
end subroutine gderivs_axi
subroutine gderivs_scal(ex,f,fx,fy,fz,dX,dY,dZ,Symmetry)
  implicit none
  interface
     subroutine gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3),              intent(in ):: ex
       real*8,dimension(ex(1),ex(2),ex(3)),intent(in ):: f
       real*8,dimension(ex(1),ex(2),ex(3)),intent(out):: fx,fy,fz
       real*8,                             intent(in ):: dX,dy,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_oct
     subroutine gderivs_eq(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_eq
!     subroutine gderivs_pi(ex,f,fx,fy,fz,dX,dY,dZ,SYM3,   &
!                     glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
!                     PI_SYMM_f)
!       implicit none
!       integer, dimension(3)                    :: ex
!       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
!       real*8                                   :: dX,dY,dZ,SYM3
!       integer                                  :: Nx,Nz
!       integer                                  :: glob_imin,glob_jmin,glob_kmin
!       real*8, dimension(Nx+1,Nz+1)             :: PI_SYMM_f
!     end subroutine gderivs_pi
     subroutine gderivs_axi(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_axi
     subroutine gderivs_gen(ex,f,fx,fy,fz,dX,dY,dZ)
       implicit none
       integer, dimension(3),              intent(in ):: ex
       real*8,dimension(ex(1),ex(2),ex(3)),intent(in ):: f
       real*8,dimension(ex(1),ex(2),ex(3)),intent(out):: fx,fy,fz
       real*8,                             intent(in ):: dX,dY,dZ
     end subroutine gderivs_gen
  end interface
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8                                   :: dX,dY,dZ
  integer                                  :: Symmetry
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  integer           :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter ( NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  !
  if(Symmetry==OCTANT) then
     call gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM, SYM, SYM)
  elseif(Symmetry==EQUATORIAL) then
     call gderivs_eq(ex,f,fx,fy,fz,dX,dY,dZ,SYM,SYM, SYM)
  elseif(Symmetry==AXISYM) then
     call gderivs_axi(ex,f,fx,fy,fz,dX,dY,dZ,SYM,SYM, SYM)
!  elseif (Symmetry==PI_SYMM) then
!     call gderivs_pi(ex,fx,fxx,fxy,fxz,dX,dY,dZ,SYM,   &
!          glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
!          PI_SYMM_fx)
  elseif(Symmetry==NO_SYMM) then
     call gderivs_gen(ex,f,fx,fy,fz,dX,dY,dZ)
  end if
end subroutine gderivs_scal
!
!****************************************************888
!
subroutine gderivs_vec(ex,fx,fy,fz,fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz, &
     dX,dY,dZ,Symmetry)
  implicit none
  interface
     subroutine gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3),              intent(in ):: ex
       real*8,dimension(ex(1),ex(2),ex(3)),intent(in ):: f
       real*8,dimension(ex(1),ex(2),ex(3)),intent(out):: fx,fy,fz
       real*8,                             intent(in ):: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_oct
     subroutine gderivs_eq(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_eq
 !    subroutine gderivs_pi(ex,f,fx,fy,fz,dX,dY,dZ,SYM3,   &
!                     glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
!                     PI_SYMM_f)
!       implicit none
!       integer, dimension(3)                    :: ex
!       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
!       real*8                                   :: dX,dY,dZ,SYM3
!       integer                                  :: Nx,Nz
!       integer                                  :: glob_imin,glob_jmin,glob_kmin
!       real*8, dimension(Nx+1,Nz+1)             :: PI_SYMM_f
!     end subroutine gderivs_pi
     subroutine gderivs_axi(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_axi
     subroutine gderivs_gen(ex,f,fx,fy,fz,dX,dY,dZ)
       implicit none
       integer, dimension(3),              intent(in ):: ex
       real*8,dimension(ex(1),ex(2),ex(3)),intent(in ):: f
       real*8,dimension(ex(1),ex(2),ex(3)),intent(out):: fx,fy,fz
       real*8,                             intent(in ):: dX,dY,dZ
     end subroutine gderivs_gen
  end interface
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz
  real*8, dimension(ex(1),ex(2),ex(3))     :: fyx,fyy,fyz
  real*8, dimension(ex(1),ex(2),ex(3))     :: fzx,fzy,fzz
  real*8                                   :: dX,dY,dZ
  integer                                  :: Symmetry
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  integer           :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter ( NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  !
  if(Symmetry==OCTANT) then
     call gderivs_oct(ex,fx,fxx,fxy,fxz,dX,dY,dZ,ANTI, SYM, SYM)
     call gderivs_oct(ex,fy,fyx,fyy,fyz,dX,dY,dZ, SYM,ANTI, SYM)
     call gderivs_oct(ex,fz,fzx,fzy,fzz,dX,dY,dZ, SYM, SYM,ANTI)
  elseif(Symmetry==EQUATORIAL) then
     call gderivs_eq(ex,fx,fxx,fxy,fxz,dX,dY,dZ,ANTI,SYM, SYM)
     call gderivs_eq(ex,fy,fyx,fyy,fyz,dX,dY,dZ,SYM, ANTI,SYM)
     call gderivs_eq(ex,fz,fzx,fzy,fzz,dX,dY,dZ,SYM, SYM, ANTI)
  elseif(Symmetry==AXISYM) then
     call gderivs_axi(ex,fx,fxx,fxy,fxz,dX,dY,dZ,ANTI,SYM, SYM)
     call gderivs_axi(ex,fy,fyx,fyy,fyz,dX,dY,dZ,ANTI, ANTI,SYM)
     call gderivs_axi(ex,fz,fzx,fzy,fzz,dX,dY,dZ,SYM, SYM, ANTI)
!  elseif (Symmetry==PI_SYMM) then
!     call gderivs_pi(ex,fx,fxx,fxy,fxz,dX,dY,dZ,SYM,   &
!          glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
!          PI_SYMM_fx)
!     call gderivs_pi(ex,fy,fyx,fyy,fyz,dX,dY,dZ,SYM,   &
!          glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
!          PI_SYMM_fy)
!     call gderivs_pi(ex,fz,fzx,fzy,fzz,dX,dY,dZ,ANTI,   &
!          glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
!          PI_SYMM_fz)
  elseif(Symmetry==NO_SYMM) then
     call gderivs_gen(ex,fx,fxx,fxy,fxz,dX,dY,dZ)
     call gderivs_gen(ex,fy,fyx,fyy,fyz,dX,dY,dZ)
     call gderivs_gen(ex,fz,fzx,fzy,fzz,dX,dY,dZ)
  end if
end subroutine gderivs_vec
!
!*************************************************************
!
subroutine gderivs_gen(ex,f,fx,fy,fz,dX,dY,dZ)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8                                   :: dX,dY,dZ
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm1,dym1,dzm1
  real*8, parameter    :: HALF = 0.5D0, TWO = 2.D0
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
! Derivative with respect to x
!
  fx(imin,:,:) = TWO * dxm1 * (f(imin+1,:,:) - f(imin,:,:))
  do i = imin + 1, imax - 1
     fx(i,:,:) = dxm1 * (f(i+1,:,:) - f(i-1,:,:))
  enddo
  fx(imax,:,:) = TWO * dxm1 * (f(imax,:,:) - f(imax-1,:,:))
!
! Derivative with respect to y
!
  fy(:,jmin,:) = TWO * dym1 * (f(:,jmin+1,:) - f(:,jmin,:))
  do j = jmin + 1, jmax - 1
     fy(:,j,:) = dym1 * (f(:,j+1,:) - f(:,j-1,:))
  enddo
  fy(:,jmax,:) = TWO * dym1 * (f(:,jmax,:) - f(:,jmax-1,:))
!
! Derivative with respect to z
!
  fz(:,:,kmin) = TWO * dzm1 * (f(:,:,kmin+1) - f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     fz(:,:,k) = dzm1 * (f(:,:,k+1) - f(:,:,k-1))
  enddo
  fz(:,:,kmax) = TWO * dzm1 * (f(:,:,kmax) - f(:,:,kmax-1))
end subroutine gderivs_gen
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Second derivatives
!
! NOTE:  We are not evaluating second derivatives at upper limits of grid!
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gdderivs(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
     dX,dY,dZ,SYM1,SYM2,SYM3)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm2,dym2,dzm2,dxm1,dym1,dzm1,dxdy,dxdz,dydz
  real*8, parameter    :: HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, F1o4 = 0.25D0
!
! Input translation
!
  imin = lbound(f,1)
  imax = ubound(f,1)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
  dxm2 = ONE/(dX*dX)
  dym2 = ONE/(dY*dY)
  dzm2 = ONE/(dZ*dZ)
  dxm1 = HALF/dX
  dym1 = HALF/dY
  dzm1 = HALF/dZ
  dxdy = F1o4/(dX*dY)
  dxdz = F1o4/(dX*dZ)
  dydz = F1o4/(dY*dZ)
!
! fxx
! 
  fxx(imin,:,:) = dxm2 * ((SYM1 - TWO)*f(imin,:,:) + f(imin+1,:,:))
  do i = imin + 1, imax - 1
     fxx(i,:,:) = dxm2 * (f(i-1,:,:) - TWO * f(i,:,:) + f(i+1,:,:))
  enddo
  fxx(imax,:,:) = 0.D0
!
! fyy
!
  fyy(:,jmin,:) = dym2 * ((SYM2 - TWO)*f(:,jmin,:) + f(:,jmin+1,:))
  do j = jmin + 1, jmax - 1
     fyy(:,j,:) = dym2 * (f(:,j-1,:) - TWO * f(:,j,:) + f(:,j+1,:))
  enddo
  fyy(:,jmax,:) = 0.D0
!
! fzz
!
  fzz(:,:,kmin) = dzm2 * ((SYM3 - TWO)*f(:,:,kmin) + f(:,:,kmin+1))
  do k = kmin + 1, kmax - 1
     fzz(:,:,k) = dzm2 * (f(:,:,k-1) - TWO * f(:,:,k) + f(:,:,k+1))
  enddo
  fzz(:,:,kmax) = 0.D0
!
! fxy
!
  fxy(imin,:,:) = dxm1 * (fy(imin+1,:,:) - SYM1*fy(imin,:,:))
  do i = imin + 1, imax - 1
     fxy(i,:,:) = dxm1 * (fy(i+1,:,:) - fy(i-1,:,:))
  end do
  fxy(imax,:,:) = 0.D0
!
! fxz
!
  fxz(imin,:,:) = dxm1 * (fz(imin+1,:,:) - SYM1*fz(imin,:,:))
  do i = imin + 1, imax - 1
     fxz(i,:,:) = dxm1 * (fz(i+1,:,:) - fz(i-1,:,:))
  end do
  fxz(imax,:,:) = 0.D0
!
! fyz
!
  fyz(:,jmin,:) = dym1 * (fz(:,jmin+1,:) - SYM2*fz(:,jmin,:))
  do j = jmin + 1, jmax - 1
     fyz(:,j,:) = dym1 * (fz(:,j+1,:) - fz(:,j-1,:))
  end do
  fyz(:,jmax,:) = 0.D0
end subroutine gdderivs
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Second derivatives
!
! NOTE:  We are not evaluating second derivatives at upper limits of grid!
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gdderivs_oct(ex,X,Y,Z,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
     dX,dY,dZ,SYM1,SYM2,SYM3)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm2,dym2,dzm2,dxm1,dym1,dzm1,dxdy,dxdz,dydz
  real*8, parameter    :: HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, F1o4 = 0.25D0
!
! Input translation
!
  imin = lbound(f,1)
  imax = ubound(f,1)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
  dxm2 = ONE/(dX*dX)
  dym2 = ONE/(dY*dY)
  dzm2 = ONE/(dZ*dZ)
  dxm1 = HALF/dX
  dym1 = HALF/dY
  dzm1 = HALF/dZ
  dxdy = F1o4/(dX*dY)
  dxdz = F1o4/(dX*dZ)
  dydz = F1o4/(dY*dZ)
  do i=imin,imax
     if(X(i,1,1) .lt. 0.D0) then
        imin = imin + 1
     end if
  end do
  do j=jmin,jmax
     if(Y(1,j,1) .lt. 0.D0) then
        jmin = jmin + 1
     end if
  end do
  do k=kmin,kmax
     if(Z(1,1,k) .lt. 0.D0) then
        kmin = kmin + 1
     end if
  end do
!
! fxx
! 
  fxx(imin,:,:) = dxm2 * ((SYM1 - TWO)*f(imin,:,:) + f(imin+1,:,:))
  do i = imin + 1, imax - 1
     fxx(i,:,:) = dxm2 * (f(i-1,:,:) - TWO * f(i,:,:) + f(i+1,:,:))
  enddo
  fxx(imax,:,:) = 0.D0
!
! fyy
!
  fyy(:,jmin,:) = dym2 * ((SYM2 - TWO)*f(:,jmin,:) + f(:,jmin+1,:))
  do j = jmin + 1, jmax - 1
     fyy(:,j,:) = dym2 * (f(:,j-1,:) - TWO * f(:,j,:) + f(:,j+1,:))
  enddo
  fyy(:,jmax,:) = 0.D0
!
! fzz
!
  fzz(:,:,kmin) = dzm2 * ((SYM3 - TWO)*f(:,:,kmin) + f(:,:,kmin+1))
  do k = kmin + 1, kmax - 1
     fzz(:,:,k) = dzm2 * (f(:,:,k-1) - TWO * f(:,:,k) + f(:,:,k+1))
  enddo
  fzz(:,:,kmax) = 0.D0
!
! fxy
!
  fxy(imin,:,:) = dxm1 * (fy(imin+1,:,:) - SYM1*fy(imin,:,:))
  do i = imin + 1, imax - 1
     fxy(i,:,:) = dxm1 * (fy(i+1,:,:) - fy(i-1,:,:))
  end do
  fxy(imax,:,:) = 0.D0
!
! fxz
!
  fxz(imin,:,:) = dxm1 * (fz(imin+1,:,:) - SYM1*fz(imin,:,:))
  do i = imin + 1, imax - 1
     fxz(i,:,:) = dxm1 * (fz(i+1,:,:) - fz(i-1,:,:))
  end do
  fxz(imax,:,:) = 0.D0
!
! fyz
!
  fyz(:,jmin,:) = dym1 * (fz(:,jmin+1,:) - SYM2*fz(:,jmin,:))
  do j = jmin + 1, jmax - 1
     fyz(:,j,:) = dym1 * (fz(:,j+1,:) - fz(:,j-1,:))
  end do
  fyz(:,jmax,:) = 0.D0
end subroutine gdderivs_oct
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Second derivatives
!
! NOTE:  Version of gdderivs for axisymmetry
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gdderivs_axi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
     dX,dY,dZ,SYM1,SYM2,SYM3)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm2,dym2,dzm2,dxm1,dym1,dzm1,dxdy,dxdz,dydz
  real*8, parameter    :: HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, F1o4 = 0.25D0
!
! Input translation
!
  imin = lbound(f,1)
  imax = ubound(f,1)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
  dxm2 = ONE/(dX*dX)
  dym2 = ONE/(dY*dY)
  dzm2 = ONE/(dZ*dZ)
  dxm1 = HALF/dX
  dym1 = HALF/dY
  dzm1 = HALF/dZ
  dxdy = F1o4/(dX*dY)
  dxdz = F1o4/(dX*dZ)
  dydz = F1o4/(dY*dZ)
!
! fxx
!
  fxx(imin,:,:) = dxm2 * ((SYM1 - TWO)*f(imin,:,:) + f(imin+1,:,:))
  do i = imin + 1, imax - 1
     fxx(i,:,:) = dxm2 * (f(i-1,:,:) - TWO * f(i,:,:) + f(i+1,:,:))
  enddo
  fxx(imax,:,:) = 0.D0
!
! fyy
!
  fyy(:,jmin,:) = dym2 * ((SYM2 - TWO)*f(:,jmin,:) + f(:,jmin+1,:))
  do j = jmin + 1, jmax - 1
     fyy(:,j,:) = dym2 * (f(:,j-1,:) - TWO * f(:,j,:) + f(:,j+1,:))
  enddo
  fyy(:,jmax,:) = 0.D0
!
! fzz
!
  fzz(:,:,kmin) = dzm2 * ((SYM3 - TWO)*f(:,:,kmin) + f(:,:,kmin+1))
  do k = kmin + 1, kmax - 1
     fzz(:,:,k) = dzm2 * (f(:,:,k-1) - TWO * f(:,:,k) + f(:,:,k+1))
  enddo
  fzz(:,:,kmax) = 0.D0
!
! fxy
!
  fxy(imin,:,:) = dxm1 * (fy(imin+1,:,:) + SYM1*fy(imin,:,:))
  do i = imin + 1, imax - 1
     fxy(i,:,:) = dxm1 * (fy(i+1,:,:) - fy(i-1,:,:))
  end do
!  write(*,*) meow1,fxy(imin+1,2,2),fy(imin+2,2,2),-fy(imin,2,2),dxm1,fy(imin+2,2,2) - fy(imin,2,2),SYM1
  fxy(imax,:,:) = 0.D0
!
! fxz
!
  fxz(imin,:,:) = dxm1 * (fz(imin+1,:,:) - SYM1*fz(imin,:,:))
  do i = imin + 1, imax - 1
     fxz(i,:,:) = dxm1 * (fz(i+1,:,:) - fz(i-1,:,:))
  end do
!  write(*,*) meow2,fxz(imin+1,2,2),fz(imin+2,2,2),fz(imin,2,2),dxm1,fz(imin+2,2,2) - fz(imin,2,2),SYM1,fy(imin,2,2),fy(imin+1,2,2)
  fxz(imax,:,:) = 0.D0
!
! fyz
!
  fyz(:,jmin,:) = dym1 * (fz(:,jmin+1,:) - SYM2*fz(:,jmin,:))
  do j = jmin + 1, jmax - 1
     fyz(:,j,:) = dym1 * (fz(:,j+1,:) - fz(:,j-1,:))
  end do
  fyz(:,jmax,:) = 0.D0
end subroutine gdderivs_axi
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Second derivatives
!
! NOTE:  Kludge for computing second derivatives of shift in compute_rhs
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gdderivs_axi_shift_compute_rhs(ex,X,Z,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
     dX,dY,dZ,SYM1,SYM2,SYM3)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm2,dym2,dzm2,dxm1,dym1,dzm1,dxdy,dxdz,dydz
  real*8, parameter    :: HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, F1o4 = 0.25D0
!
! Input translation
!
  imin = lbound(f,1)
  imax = ubound(f,1)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
  dxm2 = ONE/(dX*dX)
  dym2 = ONE/(dY*dY)
  dzm2 = ONE/(dZ*dZ)
  dxm1 = HALF/dX
  dym1 = HALF/dY
  dzm1 = HALF/dZ
  dxdy = F1o4/(dX*dY)
  dxdz = F1o4/(dX*dZ)
  dydz = F1o4/(dY*dZ)
!  do i=imin,imax
!     if(X(1,1,1) .lt. 0.D0) then 
!        imin = imin + 1
!     end if
!  end do
!  do k=kmin,kmax
!     if(Z(1,1,1) .lt. 0.D0) then 
!        kmin = kmin + 1
!     end if
!  end do
!
! fxx
!
  fxx(imin,:,:) = dxm2 * ((SYM1 - TWO)*f(imin,:,:) + f(imin+1,:,:))
  do i = imin + 1, imax - 1
     fxx(i,:,:) = dxm2 * (f(i-1,:,:) - TWO * f(i,:,:) + f(i+1,:,:))
  enddo
  fxx(imax,:,:) = 0.D0
!
! fyy
!
  fyy(:,jmin,:) = dym2 * ((SYM2 - TWO)*f(:,jmin,:) + f(:,jmin+1,:))
  do j = jmin + 1, jmax - 1
     fyy(:,j,:) = dym2 * (f(:,j-1,:) - TWO * f(:,j,:) + f(:,j+1,:))
  enddo
  fyy(:,jmax,:) = 0.D0
!
! fzz
!
  fzz(:,:,kmin) = dzm2 * ((SYM3 - TWO)*f(:,:,kmin) + f(:,:,kmin+1))
  do k = kmin + 1, kmax - 1
     fzz(:,:,k) = dzm2 * (f(:,:,k-1) - TWO * f(:,:,k) + f(:,:,k+1))
  enddo
  fzz(:,:,kmax) = 0.D0
!  write(*,*) should be the same magnitude: ,f(1,2,1),f(2,2,1),f(1,2,2),f(2,2,2)
!
! fxy
!
  fxy(imin,:,:) = dxm1 * (fy(imin+1,:,:) + SYM1*fy(imin,:,:))
!  write(*,*) meow,fxy(imin,2,2),fy(imin+1,2,2),fy(imin,2,2)
  do i = imin + 1, imax - 1
     fxy(i,:,:) = dxm1 * (fy(i+1,:,:) - fy(i-1,:,:))
  end do
!  write(*,*) meow1,fxy(imin+1,2,2),fy(imin+2,2,2),-fy(imin,2,2),dxm1,fy(imin+2,2,2) - fy(imin,2,2),SYM1,fy(imin,2,2),fy(imin+1,2,2)
  fxy(imax,:,:) = 0.D0
!
! fxz
!
  fxz(imin,:,:) = dxm1 * (fz(imin+1,:,:) - SYM1*fz(imin,:,:))
  do i = imin + 1, imax - 1
     fxz(i,:,:) = dxm1 * (fz(i+1,:,:) - fz(i-1,:,:))
  end do
!  write(*,*) meow2,fxz(imin+1,2,2),fz(imin+2,2,2),fz(imin,2,2),dxm1,fz(imin+2,2,2) - fz(imin,2,2),SYM1
  fxz(imax,:,:) = 0.D0
!
! fyz
!
  fyz(:,jmin,:) = dym1 * (fz(:,jmin+1,:) - SYM2*fz(:,jmin,:))
  do j = jmin + 1, jmax - 1
     fyz(:,j,:) = dym1 * (fz(:,j+1,:) - fz(:,j-1,:))
  end do
  fyz(:,jmax,:) = 0.D0
end subroutine gdderivs_axi_shift_compute_rhs
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Second derivatives
!
! NOTE:  version of gdderivs for Pi symmetry
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gdderivs_pi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
                             dX,dY,dZ,SYM1,SYM3,  &
                        glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
                        PI_SYMM_f)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8                                   :: dX,dY,dZ,SYM1,SYM3
  integer                                  :: Nx,Nz
  integer                                  :: glob_imin,glob_jmin,glob_kmin
  real*8, dimension(Nx+1,Nz+1)                   :: PI_SYMM_f
  !
  ! other variables
  !
  integer              :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: dxm2,dym2,dzm2,dxm1,dym1,dzm1,dxdy,dxdz,dydz
  real*8, parameter    :: HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, F1o4 = 0.25D0
!
! Input translation
!
  imin = lbound(f,1)
  imax = ubound(f,1)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
  dxm2 = ONE/(dX*dX)
  dym2 = ONE/(dY*dY)
  dzm2 = ONE/(dZ*dZ)
  dxm1 = HALF/dX
  dym1 = HALF/dY
  dzm1 = HALF/dZ
  dxdy = F1o4/(dX*dY)
  dxdz = F1o4/(dX*dZ)
  dydz = F1o4/(dY*dZ)
!
! fxx
! 
  fxx(imin,:,:) = 0.D0
  do i = imin + 1, imax - 1
     fxx(i,:,:) = dxm2 * (f(i-1,:,:) - TWO * f(i,:,:) + f(i+1,:,:))
  enddo
  fxx(imax,:,:) = 0.D0
!
! fyy
!
  do k = kmin, kmax
    do i = imin, imax
      fyy(i,jmin,k) = dym2 * (f(i,jmin+1,k) - TWO*f(i,jmin,k)  &
                + PI_SYMM_f(i+glob_imin,k+glob_kmin))
    end do
  end do
  do j = jmin + 1, jmax - 1
     fyy(:,j,:) = dym2 * (f(:,j-1,:) - TWO * f(:,j,:) + f(:,j+1,:))
  enddo
  fyy(:,jmax,:) = 0.D0
!
! fzz
!
  fzz(:,:,kmin) = dzm2 * ((SYM3 - TWO)*f(:,:,kmin) + f(:,:,kmin+1))
  do k = kmin + 1, kmax - 1
     fzz(:,:,k) = dzm2 * (f(:,:,k-1) - TWO * f(:,:,k) + f(:,:,k+1))
  enddo
  fzz(:,:,kmax) = 0.D0
!
! fxy
!
  fxy(imin,:,:) = TWO * dxm1 * (fy(imin+1,:,:) - fy(imin,:,:))
  do i = imin + 1, imax - 1
     fxy(i,:,:) = dxm1 * (fy(i+1,:,:) - fy(i-1,:,:))
  end do
  fxy(imax,:,:) = 0.D0
!
! fxz
!
  fxz(imin,:,:) = TWO * dxm1 * (fz(imin+1,:,:) - fz(imin,:,:))
  do i = imin + 1, imax - 1
     fxz(i,:,:) = dxm1 * (fz(i+1,:,:) - fz(i-1,:,:))
  end do
  fxz(imax,:,:) = 0.D0
!
! fyz
!
  do i = imin, imax
      fyz(i,jmin,kmin) = dym1 * dzm1 *  &
        ( f(i,jmin+1,kmin+1)-SYM3*f(i,jmin+1,kmin) &
         -PI_SYMM_f(i+glob_imin,kmin+1+glob_kmin)     &
         +SYM3*PI_SYMM_f(i+glob_imin,kmin+glob_kmin) )
  end do
  do k = kmin+1, kmax-1
    do i = imin, imax
      fyz(i,jmin,k) = dym1 * dzm1 *  ( f(i,jmin+1,k+1)-f(i,jmin+1,k-1) &
         -PI_SYMM_f(i+glob_imin,k+1+glob_kmin) &
         +PI_SYMM_f(i+glob_imin,k-1+glob_kmin) )
    end do
  end do
  do i = imin, imax
     fyz(i,jmin,kmax) = dym1 * dzm1 * TWO* (f(i,jmin+1,k)-f(i,jmin+1,k-1) &
         -PI_SYMM_f(i+glob_imin,k  +glob_kmin) &
         +PI_SYMM_f(i+glob_imin,k-1+glob_kmin) )
  end do
!    write(*,*) gdderivs_pi fyz(6,1,10) 
!    write(*,11) fyz(6,1,10),(f(6,2,10)-f(6,2,9))*dzm1,  &
!        (TWO*PI_SYMM_f(6+glob_imin,10  +glob_kmin) &
!         +TWO*PI_SYMM_f(6+glob_imin,9 +glob_kmin))*dzm1
!11 format(6f20.15)
  do j = jmin + 1, jmax - 1
     fyz(:,j,:) = dym1 * (fz(:,j+1,:) - fz(:,j-1,:))
  end do
  fyz(:,jmax,:) = 0.D0
end subroutine gdderivs_pi
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Laplace operator:       nablaf = g^{ij} f_{,ij}
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine glaplace(ex, f, fx, fy, fz, nablaf, dX, dY, dZ, SYM1, SYM2, SYM3, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
  implicit none
  interface
     subroutine gdderivs(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
          dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gdderivs
  end interface
  !
  ! Input parameters:
  !
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: f, fx, fy, fz, nablaf
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx, gupxy, gupxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy, gupyz, gupzz
  integer                                     :: i, j, k
  real*8                                      :: dX, dY, dZ
  real*8                                      :: SYM1, SYM2, SYM3
  !
  ! other variables
  !
  real*8, dimension(ex(1),ex(2),ex(3))        :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8, parameter                           :: TWO = 2.D0
  !
  ! compute second derivatives of f...
  !
  call gdderivs(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ,SYM1,SYM2,SYM3)
  !
  ! and contract with g:
  !
  nablaf = gupxx * fxx + gupyy * fyy + gupzz * fzz + &
       TWO * ( gupxy * fxy + gupxz * fxz + gupyz * fyz )
end subroutine glaplace
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Laplace operator:       nablaf = g^{ij} f_{,ij}
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine glaplace_oct(ex,X,Y,Z, f, fx, fy, fz, nablaf, dX, dY, dZ, SYM1, SYM2, SYM3, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
  implicit none
  interface
     subroutine gdderivs_oct(ex,X,Y,Z,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
          dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gdderivs_oct
  end interface
  !
  ! Input parameters:
  !
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: f, fx, fy, fz, nablaf
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx, gupxy, gupxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy, gupyz, gupzz
  integer                                     :: i, j, k
  real*8                                      :: dX, dY, dZ
  real*8                                      :: SYM1, SYM2, SYM3
  !
  ! other variables
  !
  real*8, dimension(ex(1),ex(2),ex(3))        :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8, parameter                           :: TWO = 2.D0
  !
  ! compute second derivatives of f...
  !
  call gdderivs_oct(ex,X,Y,Z,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ,SYM1,SYM2,SYM3)
  !
  ! and contract with g:
  !
  nablaf = gupxx * fxx + gupyy * fyy + gupzz * fzz + &
       TWO * ( gupxy * fxy + gupxz * fxz + gupyz * fyz )
end subroutine glaplace_oct
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Laplace operator for Pi symmetry:       nablaf = g^{ij} f_{,ij}
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine glaplace_pi(ex, f, fx, fy, fz, nablaf, dX, dY, dZ, SYM1, SYM3, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,  &
     glob_imin,glob_jmin,glob_kmin,Nx,Nz, PI_SYMM_f)
  implicit none
  interface
     subroutine gdderivs_pi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
          dX,dY,dZ,SYM1,SYM3,glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
          PI_SYMM_f)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
       real*8                                   :: dX,dY,dZ,SYM1,SYM3
       integer                                  :: Nx,Nz
       integer                                  :: glob_imin,glob_jmin,glob_kmin
       real*8, dimension(Nx+1,Nz+1)             :: PI_SYMM_f
     end subroutine gdderivs_pi
  end interface
  !
  ! Input parameters:
  !
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: f, fx, fy, fz, nablaf
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx, gupxy, gupxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy, gupyz, gupzz
  integer                                     :: i, j, k
  real*8                                      :: dX, dY, dZ
  real*8                                      :: SYM1, SYM3
  integer                                     :: Nx,Nz
  integer                                     :: glob_imin,glob_jmin,glob_kmin
  real*8, dimension(Nx+1,Nz+1)                :: PI_SYMM_f
  !
  ! other variables
  !
  real*8, dimension(ex(1),ex(2),ex(3))        :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8, parameter                           :: TWO = 2.D0
  !
  ! compute second derivatives of f...
  !
  call gdderivs_pi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ,SYM1,SYM3,  &
        glob_imin,glob_jmin,glob_kmin,Nx,Nz,PI_SYMM_f)
  !
  ! and contract with g:
  !
  nablaf = gupxx * fxx + gupyy * fyy + gupzz * fzz + &
       TWO * ( gupxy * fxy + gupxz * fxz + gupyz * fyz )
end subroutine glaplace_pi
!-----------------------------------------------------------------------------
!
! Laplace operator in axisymmetry
!
!-----------------------------------------------------------------------------
subroutine glaplace_axi(ex,X,Z,f,fx,fy,fz,nablaf,dX,dY,dZ,SYM1,SYM2,SYM3, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
  implicit none
  interface
     subroutine gdderivs_axi_shift_compute_rhs(ex,X,Z,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
          dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: X,Z
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gdderivs_axi_shift_compute_rhs
  end interface
  !
  ! Input parameters:
  !
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: f, fx, fy, fz, nablaf
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx, gupxy, gupxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy, gupyz, gupzz
  integer                                     :: i, j, k
  real*8                                      :: dX, dY, dZ
  real*8                                      :: SYM1, SYM2, SYM3
  !
  ! other variables
  !
  real*8, dimension(ex(1),ex(2),ex(3))        :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8, parameter                           :: TWO = 2.D0
  if(X(1,1,1).lt.0.D0) then
!     fy(1,:,:) = fy(2,:,:)
  end if
  if(Z(1,1,1).lt.0.D0) then
!     fy(:,:,1) = fy(:,:,2)
  end if
  !
  ! compute second derivatives of f...
  !
  call gdderivs_axi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ,SYM1, &
       SYM2,SYM3)
!  call gdderivs_axi_shift_compute_rhs(ex,X,Z,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ,SYM1, &
!       SYM2,SYM3)
  !
  ! and contract with g:
  !
  nablaf = gupxx * fxx + gupyy * fyy + gupzz * fzz + &
       TWO * ( gupxy * fxy + gupxz * fxz + gupyz * fyz )
!  write(*,*) hi. glaplace_axi:,fy(1,2,15),fy(2,2,15),f(2,1,15),f(2,3,15),fxx(2,2,15),fyy(2,2,15),fzz(2,2,15),fxy(2,2,15),fxz(2,2,15
end subroutine glaplace_axi
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! covariant derivatives w.r.t *tilded metric* of *symmetric* variable
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gcov_deriv(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ, &
     Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
     Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
     Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz)
  implicit none
  interface
     subroutine gdderivs(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
          dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gdderivs
  end interface
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: f
  real*8, dimension(ex(1),ex(2),ex(3))        :: fx, fy, fz
  real*8, dimension(ex(1),ex(2),ex(3))        :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamzyy, Gamzyz, Gamzzz
  real*8                                      :: dX, dY, dZ
!
! Other variables:
!
  real*8, parameter                           :: SYM = 1.0
!
! Compute second derivatives of f...
!
  call gdderivs(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ,SYM,SYM,SYM)
!
! Add Connection terms
!
  fxx = fxx - Gamxxx*fx - Gamyxx*fy - Gamzxx*fz
  fxy = fxy - Gamxxy*fx - Gamyxy*fy - Gamzxy*fz
  fxz = fxz - Gamxxz*fx - Gamyxz*fy - Gamzxz*fz
  fyy = fyy - Gamxyy*fx - Gamyyy*fy - Gamzyy*fz
  fyz = fyz - Gamxyz*fx - Gamyyz*fy - Gamzyz*fz
  fzz = fzz - Gamxzz*fx - Gamyzz*fy - Gamzzz*fz
end subroutine gcov_deriv
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! covariant derivatives w.r.t *tilded metric* of *symmetric* variable
! for Pi symmetry
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gcov_deriv_pi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ, &
     Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
     Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
     Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
     glob_imin,glob_jmin,glob_kmin,Nx,Nz,PI_SYMM_f)
  implicit none
  interface
     subroutine gdderivs_pi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
          dX,dY,dZ,SYM1,SYM3,glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
          PI_SYMM_f)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
       real*8                                   :: dX,dY,dZ,SYM1,SYM3
       integer                                  :: Nx,Nz
       integer                                  :: glob_imin,glob_jmin,glob_kmin
       real*8, dimension(Nx+1,Nz+1)             :: PI_SYMM_f
     end subroutine gdderivs_pi
  end interface
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: f
  real*8, dimension(ex(1),ex(2),ex(3))        :: fx, fy, fz
  real*8, dimension(ex(1),ex(2),ex(3))        :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamzyy, Gamzyz, Gamzzz
  real*8                                      :: dX, dY, dZ
  integer                                     :: Nx,Nz
  integer                                     :: glob_imin,glob_jmin,glob_kmin
  real*8, dimension(Nx+1,Nz+1)                      :: PI_SYMM_f
!
! Other variables:
!
  real*8, parameter                           :: SYM = 1.0
!
! Compute second derivatives of f...
!
  call gdderivs_pi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,  &
        dX,dY,dZ,SYM,SYM,glob_imin,glob_jmin,glob_kmin,Nx,Nz,PI_SYMM_f)
!
! Add Connection terms
!
  fxx = fxx - Gamxxx*fx - Gamyxx*fy - Gamzxx*fz
  fxy = fxy - Gamxxy*fx - Gamyxy*fy - Gamzxy*fz
  fxz = fxz - Gamxxz*fx - Gamyxz*fy - Gamzxz*fz
  fyy = fyy - Gamxyy*fx - Gamyyy*fy - Gamzyy*fz
  fyz = fyz - Gamxyz*fx - Gamyyz*fy - Gamzyz*fz
  fzz = fzz - Gamxzz*fx - Gamyzz*fy - Gamzzz*fz
end subroutine gcov_deriv_pi
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! covariant derivatives w.r.t. *untilded* (i.e. physical) metric 
! of *symmetric* variable
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gnt_cov_deriv(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ, &
     phix,phiy,phiz,gxx,gxy,gxz,gyy,gyz,gzz, &
     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
     Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
     Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
     Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz)
  implicit none
  interface
     subroutine gdderivs(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
          dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gdderivs
  end interface
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3))        :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: phix,phiy,phiz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy,gupyz,gupzz
  real*8                                      :: dX, dY, dZ
!
! Other variables:
!
  real*8, dimension(ex(1),ex(2),ex(3))        :: phiupx,phiupy,phiupz
  real*8, parameter                           :: SYM = 1.0, TWO = 2.D0
!
! Compute second derivatives of f...
!
  call gdderivs(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ,SYM,SYM,SYM)
!
! Make untilded Gammas out of tilded ones - first raise index on phi_i...
!
  phiupx = gupxx*phix + gupxy*phiy + gupxz*phiz
  phiupy = gupxy*phix + gupyy*phiy + gupyz*phiz
  phiupz = gupxz*phix + gupyz*phiy + gupzz*phiz
!
! ... and then add reconstructed *untilded* Christofels...
!
  fxx = fxx - &
       (Gamxxx + TWO * (phix + phix - gxx * phiupx) ) * fx - &
       (Gamyxx + TWO * (            - gxx * phiupy) ) * fy - &
       (Gamzxx + TWO * (            - gxx * phiupz) ) * fz
  fxy = fxy - &
       (Gamxxy + TWO * (phiy        - gxy * phiupx) ) * fx - &
       (Gamyxy + TWO * (       phix - gxy * phiupy) ) * fy - &
       (Gamzxy + TWO * (            - gxy * phiupz) ) * fz
  fxz = fxz - &
       (Gamxxz + TWO * (phiz        - gxz * phiupx) ) * fx - &
       (Gamyxz + TWO * (            - gxz * phiupy) ) * fy - &
       (Gamzxz + TWO * (       phix - gxz * phiupz) ) * fz
  fyy = fyy - &
       (Gamxyy + TWO * (            - gyy * phiupx) ) * fx - &
       (Gamyyy + TWO * (phiy + phiy - gyy * phiupy) ) * fy - &
       (Gamzyy + TWO * (            - gyy * phiupz) ) * fz
  fyz = fyz - &
       (Gamxyz + TWO * (            - gyz * phiupx) ) * fx - &
       (Gamyyz + TWO * (phiz        - gyz * phiupy) ) * fy - &
       (Gamzyz + TWO * (       phiy - gyz * phiupz) ) * fz
  fzz = fzz - &
       (Gamxzz + TWO * (            - gzz * phiupx) ) * fx - &
       (Gamyzz + TWO * (            - gzz * phiupy) ) * fy - &
       (Gamzzz + TWO * (phiz + phiz - gzz * phiupz) ) * fz
end subroutine gnt_cov_deriv
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! covariant derivatives w.r.t. *untilded* (i.e. physical) metric 
! of *symmetric* variable
! for Pi symmetry
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gnt_cov_deriv_pi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,dX,dY,dZ, &
     phix,phiy,phiz,gxx,gxy,gxz,gyy,gyz,gzz, &
     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
     Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
     Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
     Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
     glob_imin,glob_jmin,glob_kmin,Nx,Nz,PI_SYMM_f)
  implicit none
  interface
     subroutine gdderivs_pi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz, &
          dX,dY,dZ,SYM1,SYM3,glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
          PI_SYMM_f)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8, dimension(ex(1),ex(2),ex(3))     :: fxx,fxy,fxz,fyy,fyz,fzz
       real*8                                   :: dX,dY,dZ,SYM1,SYM3
       integer                                  :: Nx,Nz
       integer                                  :: glob_imin,glob_jmin,glob_kmin
       real*8, dimension(Nx+1,Nz+1)             :: PI_SYMM_f
     end subroutine gdderivs_pi
  end interface
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3))        :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: phix,phiy,phiz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy,gupyz,gupzz
  real*8                                      :: dX, dY, dZ
  integer                                     :: Nx,Nz
  integer                                     :: glob_imin,glob_jmin,glob_kmin
  real*8, dimension(Nx+1,Nz+1)                      :: PI_SYMM_f
!
! Other variables:
!
  real*8, dimension(ex(1),ex(2),ex(3))        :: phiupx,phiupy,phiupz
  real*8, parameter                           :: SYM = 1.0, TWO = 2.D0
!
! Compute second derivatives of f...
!
  call gdderivs_pi(ex,f,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,  &
        dX,dY,dZ,SYM,SYM,glob_imin,glob_jmin,glob_kmin,Nx,Nz,PI_SYMM_f)
!
! Make untilded Gammas out of tilded ones - first raise index on phi_i...
!
  phiupx = gupxx*phix + gupxy*phiy + gupxz*phiz
  phiupy = gupxy*phix + gupyy*phiy + gupyz*phiz
  phiupz = gupxz*phix + gupyz*phiy + gupzz*phiz
!
! ... and then add reconstructed *untilded* Christofels...
!
  fxx = fxx - &
       (Gamxxx + TWO * (phix + phix - gxx * phiupx) ) * fx - &
       (Gamyxx + TWO * (            - gxx * phiupy) ) * fy - &
       (Gamzxx + TWO * (            - gxx * phiupz) ) * fz
  fxy = fxy - &
       (Gamxxy + TWO * (phiy        - gxy * phiupx) ) * fx - &
       (Gamyxy + TWO * (       phix - gxy * phiupy) ) * fy - &
       (Gamzxy + TWO * (            - gxy * phiupz) ) * fz
  fxz = fxz - &
       (Gamxxz + TWO * (phiz        - gxz * phiupx) ) * fx - &
       (Gamyxz + TWO * (            - gxz * phiupy) ) * fy - &
       (Gamzxz + TWO * (       phix - gxz * phiupz) ) * fz
  fyy = fyy - &
       (Gamxyy + TWO * (            - gyy * phiupx) ) * fx - &
       (Gamyyy + TWO * (phiy + phiy - gyy * phiupy) ) * fy - &
       (Gamzyy + TWO * (            - gyy * phiupz) ) * fz
  fyz = fyz - &
       (Gamxyz + TWO * (            - gyz * phiupx) ) * fx - &
       (Gamyyz + TWO * (phiz        - gyz * phiupy) ) * fy - &
       (Gamzyz + TWO * (       phiy - gyz * phiupz) ) * fz
  fzz = fzz - &
       (Gamxzz + TWO * (            - gzz * phiupx) ) * fx - &
       (Gamyzz + TWO * (            - gzz * phiupy) ) * fy - &
       (Gamzzz + TWO * (phiz + phiz - gzz * phiupz) ) * fz
end subroutine gnt_cov_deriv_pi
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!  Interface to compute the Covariant Laplace operator of f 
!
!       nablaf = g^{ij} f_{,ij} - Gamma^i f_{,i}
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gdd(ex,f,ddf, dX, dY, dZ, &
     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Gamx,Gamy,Gamz)
  implicit none
  interface
     subroutine glaplace(ex,f,fx,fy,fz,nablaf,dX,dY,dZ,SYM1,SYM2,SYM3, &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
       integer, dimension(3)                       :: ex
       real*8, dimension(ex(1),ex(2),ex(3))        :: f, fx, fy, fz, nablaf
       real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx, gupxy, gupxz
       real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy, gupyz, gupzz
       integer                                     :: i, j, k
       real*8                                      :: dX, dY, dZ
       real*8                                      :: SYM1, SYM2, SYM3
     end subroutine glaplace
     subroutine gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_oct
  end interface
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f
  real*8                                   :: dX, dY, dZ
  real*8, dimension(ex(1),ex(2),ex(3))     :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3))     :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3))     :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3))     :: ddf
!
! Other variables
!
  real*8, dimension(ex(1),ex(2),ex(3))     :: fx,fy,fz
  real*8, parameter                        :: SYM = 1.D0
!
! Compute first derivatives of f:
!
  call gderivs_oct(ex, f, fx, fy, fz, dX, dY, dZ, SYM, SYM, SYM)
!
! Compute Laplace operator:
!
  call glaplace(ex,f,fx,fy,fz,ddf,dX,dY,dZ,SYM,SYM,SYM, &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
!
! Add connection terms:
!
  ddf = ddf - Gamx * fx - Gamy * fy - Gamz * fz
end subroutine gdd
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!  Interface to compute the Covariant Laplace operator of f 
!
!       nablaf = g^{ij} f_{,ij} - Gamma^i f_{,i}
!
!  for Pi symmetry
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gdd_pi(ex,f,ddf, dX, dY, dZ, &
     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Gamx,Gamy,Gamz,  &
     glob_imin,glob_jmin,glob_kmin,Nx,Nz,PI_SYMM_f)
  implicit none
  interface
     subroutine glaplace_pi(ex,f,fx,fy,fz,nablaf,dX,dY,dZ,SYM1,SYM3, &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,  &
               glob_imin,glob_jmin,glob_kmin,Nx,Nz,PI_SYMM_f)
       integer, dimension(3)                       :: ex
       real*8, dimension(ex(1),ex(2),ex(3))        :: f, fx, fy, fz, nablaf
       real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx, gupxy, gupxz
       real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy, gupyz, gupzz
       integer                                     :: i, j, k
       real*8                                      :: dX, dY, dZ
       real*8                                      :: SYM1, SYM3
       integer                                     :: Nx,Nz
       integer                                     :: glob_imin,glob_jmin,glob_kmin
       real*8, dimension(Nx+1,Nz+1)                :: PI_SYMM_f
     end subroutine glaplace_pi
     subroutine gderivs_pi(ex,f,fx,fy,fz,dX,dY,dZ,SYM3,  &
                        glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
                        PI_SYMM_f)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM3
       integer                                  :: Nx,Nz
       integer                                  :: glob_imin,glob_jmin,glob_kmin
       real*8, dimension(Nx+1,Nz+1)                :: PI_SYMM_f
     end subroutine gderivs_pi
  end interface
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f
  real*8                                   :: dX, dY, dZ
  real*8, dimension(ex(1),ex(2),ex(3))     :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3))     :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3))     :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3))     :: ddf
  integer                                  :: Nx,Nz
  integer                                  :: glob_imin,glob_jmin,glob_kmin
  real*8, dimension(Nx+1,Nz+1)                   :: PI_SYMM_f
!
! Other variables
!
  real*8, dimension(ex(1),ex(2),ex(3))     :: fx,fy,fz
  real*8, parameter                        :: SYM = 1.D0
!
! Compute first derivatives of f:
!
  call gderivs_pi(ex, f, fx, fy, fz, dX, dY, dZ, SYM,  &
          glob_imin,glob_jmin,glob_kmin,Nx,Nz, PI_SYMM_f)
!
!
! Compute Laplace operator:
!
  call glaplace_pi(ex,f,fx,fy,fz,ddf,dX,dY,dZ,SYM,SYM, &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,  &
          glob_imin,glob_jmin,glob_kmin,Nx,Nz, PI_SYMM_f)
!
!
! Add connection terms:
!
  ddf = ddf - Gamx * fx - Gamy * fy - Gamz * fz
end subroutine gdd_pi
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!  Interface to compute the non-tilded Covariant Laplace operator of f 
!
!       nablaf = e^{-4 phi} ( g^{ij} f_{,ij} - Gamma^i f_{,i} +
!                             2 g^{ij} phi_{,i} f_{,j} )
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gntdd(ex,f,ddf, dX, dY, dZ, phi, phix, phiy, phiz, &
     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Gamx,Gamy,Gamz)
  implicit none
  interface
     subroutine glaplace(ex,f,fx,fy,fz,nablaf,dX,dY,dZ,SYM1,SYM2,SYM3, &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
       integer, dimension(3)                       :: ex
       real*8, dimension(ex(1),ex(2),ex(3))        :: f, fx, fy, fz, nablaf
       real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx, gupxy, gupxz
       real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy, gupyz, gupzz
       integer                                     :: i, j, k
       real*8                                      :: dX, dY, dZ
       real*8                                      :: SYM1, SYM2, SYM3
     end subroutine glaplace
     subroutine gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_oct
  end interface
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f
  real*8                                   :: dX, dY, dZ
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi,phix,phiy,phiz
  real*8, dimension(ex(1),ex(2),ex(3))     :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3))     :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3))     :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3))     :: ddf
!
! Other variables
!
  real*8, dimension(ex(1),ex(2),ex(3))     :: fx,fy,fz
  real*8, parameter                        :: SYM  = 1.D0
  real*8, parameter                        :: TWO  = 2.D0
  real*8, parameter                        :: FOUR = 4.D0
!
! Compute first derivatives of f:
!
  call gderivs_oct(ex, f, fx, fy, fz, dX, dY, dZ, SYM, SYM, SYM)
!
! Compute Laplace operator:
!
  call glaplace(ex,f,fx,fy,fz,ddf,dX,dY,dZ,SYM,SYM,SYM, &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
!
! Add connection terms:
!
  ddf = exp(-FOUR*phi) * ( ddf - Gamx * fx - Gamy * fy - Gamz * fz &
       + TWO * (gupxx * fx * phix + gupyy * fy * phiy + gupzz * fz * phiz &
       + gupxy * (fx * phiy + fy * phix ) &
       + gupxz * (fx * phiz + fz * phix ) &
       + gupyz * (fy * phiz + fz * phiy ) ) )
!write(*,*) higntdd:,phi(2,33,2),ddf(2,33,2),Gamx(2,33,2),fx(2,33,2),Gamy(2,33,2),fy(2,33,2),Gamz(2,33,2),fz(2,33,2)
end subroutine gntdd
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!  Interface to compute the non-tilded Covariant Laplace operator of f 
!
!       nablaf = e^{-4 phi} ( g^{ij} f_{,ij} - Gamma^i f_{,i} +
!                             2 g^{ij} phi_{,i} f_{,j} )
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine gntdd_pi(ex,f,ddf, dX, dY, dZ, phi, phix, phiy, phiz, &
     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Gamx,Gamy,Gamz,  &
     glob_imin,glob_jmin,glob_kmin,Nx,Nz,PI_SYMM_f)
  implicit none
  interface
     subroutine glaplace_pi(ex,f,fx,fy,fz,nablaf,dX,dY,dZ,SYM1,SYM3, &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,  &
               glob_imin,glob_jmin,glob_kmin,Nx,Nz,PI_SYMM_f)
       integer, dimension(3)                       :: ex
       real*8, dimension(ex(1),ex(2),ex(3))        :: f, fx, fy, fz, nablaf
       real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx, gupxy, gupxz
       real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy, gupyz, gupzz
       integer                                     :: i, j, k
       real*8                                      :: dX, dY, dZ
       real*8                                      :: SYM1, SYM3
       integer                                     :: Nx,Nz
       integer                                     :: glob_imin,glob_jmin,glob_kmin
       real*8, dimension(Nx+1,Nz+1)                :: PI_SYMM_f
     end subroutine glaplace_pi
     subroutine gderivs_pi(ex,f,fx,fy,fz,dX,dY,dZ,SYM3,  &
                        glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
                        PI_SYMM_f)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM3
       integer                                  :: Nx,Nz
       integer                                  :: glob_imin,glob_jmin,glob_kmin
       real*8, dimension(Nx+1,Nz+1)                :: PI_SYMM_f
     end subroutine gderivs_pi
  end interface
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f
  real*8                                   :: dX, dY, dZ
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi,phix,phiy,phiz
  real*8, dimension(ex(1),ex(2),ex(3))     :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3))     :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3))     :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3))     :: ddf
  integer                                  :: Nx,Nz
  integer                                  :: glob_imin,glob_jmin,glob_kmin
  real*8, dimension(Nx+1,Nz+1)                   :: PI_SYMM_f
!
! Other variables
!
  real*8, dimension(ex(1),ex(2),ex(3))     :: fx,fy,fz
  real*8, parameter                        :: SYM  = 1.D0
  real*8, parameter                        :: TWO  = 2.D0
  real*8, parameter                        :: FOUR = 4.D0
!
! Compute first derivatives of f:
!
  call gderivs_pi(ex, f, fx, fy, fz, dX, dY, dZ, SYM,  &
          glob_imin,glob_jmin,glob_kmin,Nx,Nz, PI_SYMM_f)
!
! Compute Laplace operator:
!
  call glaplace_pi(ex,f,fx,fy,fz,ddf,dX,dY,dZ,SYM,SYM, &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,  &
          glob_imin,glob_jmin,glob_kmin,Nx,Nz, PI_SYMM_f)
!
! Add connection terms:
!
  ddf = exp(-FOUR*phi) * ( ddf - Gamx * fx - Gamy * fy - Gamz * fz &
       + TWO * (gupxx * fx * phix + gupyy * fy * phiy + gupzz * fz * phiz &
       + gupxy * (fx * phiy + fy * phix ) &
       + gupxz * (fx * phiz + fz * phix ) &
       + gupyz * (fy * phiz + fz * phiy ) ) )
end subroutine gntdd_pi
!-----------------------------------------------------------------------------
!
! Computes partial derivative of function f with respect to x
!
!-----------------------------------------------------------------------------
subroutine ddx(ext,f,df,dX,sym)
!
! Input parameters:
!
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
  real*8                                   :: sym,dX
!
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
!
  integer                    :: i
  integer                    :: imin, imax
  real*8                     :: dxm1
  real*8                     :: HALF
  parameter(HALF = 0.5D0)
!
! Input translation
!
  imin = lbound(f,1)
  imax = ubound(f,1)
!
! Derivative with respect to x
!
  dxm1 = HALF/dX
  df(imin,:,:) = dxm1 * (f(imin+1,:,:) - sym*f(imin,:,:))
  do i = imin + 1, imax - 1
     df(i,:,:) = dxm1 * (f(i+1,:,:) - f(i-1,:,:))
  enddo
!  df(imax,:,:) = dxm1 * (sym*f(imax,:,:) - f(imax-1,:,:))
  df(imax,:,:) = 2.d0*dxm1 * (f(imax,:,:) - f(imax-1,:,:))
end subroutine ddx
!-----------------------------------------------------------------------------
!
! Computes partial derivative of function f with respect to y
!
!-----------------------------------------------------------------------------
subroutine ddy(ext,f,df,dY,sym)
!
! Input parameters:
!
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
  real*8                                   :: sym,dY
!
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
!
  integer                    :: j
  integer                    :: jmin, jmax
  real*8                     :: dym1
  real*8                     :: HALF
  parameter(HALF = 0.5D0)
!
! Input translation
!
  jmin = lbound(f,2)
  jmax = ubound(f,2)
!
! Derivative with respect to y
!
  dym1 = HALF/dY
  df(:,jmin,:) = dym1 * (f(:,jmin+1,:) - sym*f(:,jmin,:))
  do j = jmin + 1, jmax - 1
     df(:,j,:) = dym1 * (f(:,j+1,:) - f(:,j-1,:))
  enddo
!  df(:,jmax,:) = dym1 * (sym*f(:,jmax,:) - f(:,jmax-1,:))
  df(:,jmax,:) = 2.d0 * dym1 * (f(:,jmax,:) - f(:,jmax-1,:))
end subroutine ddy
!-----------------------------------------------------------------------------
!
! Computes partial derivative of function f with respect to y
!
!-----------------------------------------------------------------------------
subroutine ddy_pi(ext,f,df,dY,sym,  &
                glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
                PI_SYMM_f)
!
! Input parameters:
!
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
  real*8                                   :: dY,sym
  integer                                  :: Nx,Nz
  integer                                  :: glob_imin,glob_jmin,glob_kmin
  real*8, dimension(Nx+1,Nz+1)                   :: PI_SYMM_f
  !
  ! other variables
  !
  integer                    :: i,j,k
  integer                    :: imin,jmin,kmin,imax,jmax,kmax
  real*8                     :: dym1
  real*8                     :: HALF
  parameter(HALF = 0.5D0)
!
! Input translation
!
  imin = lbound(f,1)
  imax = ubound(f,1)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
!
! Derivative with respect to y
!
  dym1 = HALF/dY
  do k = kmin, kmax
    do i = imin, imax
      df(i,jmin,k) = dym1 * (f(i,jmin+1,k) - PI_SYMM_f(i+glob_imin,k+glob_kmin))
    end do
  end do
  do j = jmin + 1, jmax - 1
     df(:,j,:) = dym1 * (f(:,j+1,:) - f(:,j-1,:))
  enddo
  df(:,jmax,:) = dym1 * (sym*f(:,jmax,:) - f(:,jmax-1,:))
end subroutine ddy_pi
!-----------------------------------------------------------------------------
!
! Computes partial derivative of function f with respect to z
!
!-----------------------------------------------------------------------------
subroutine ddz(ext,f,df,dZ,sym)
!
! Input parameters:
!
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
  real*8                                   :: sym,dZ
!
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
!
  integer                    :: k
  integer                    :: kmin, kmax
  real*8                     :: dzm1
  real*8                     :: HALF
  parameter(HALF = 0.5D0)
!
! Input translation
!
  kmin = lbound(f,3)
  kmax = ubound(f,3)
!
! Derivative with respect to z
!
  dzm1 = HALF/dZ
  df(:,:,kmin) = dzm1 * (f(:,:,kmin+1) - sym*f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     df(:,:,k) = dzm1 * (f(:,:,k+1) - f(:,:,k-1))
  enddo
!  df(:,:,kmax) = dzm1 * (sym*f(:,:,kmax) - f(:,:,kmax-1))
  df(:,:,kmax) = 2.d0 * dzm1 * (f(:,:,kmax) - f(:,:,kmax-1))
end subroutine ddz
!-----------------------------------------------------------------------------
!
! Computes second partial derivative of function f with respect to x
!
!-----------------------------------------------------------------------------
subroutine ddx2(ext,f,d2f,dX,sym)
!
! Input parameters:
!
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,d2f
  real*8                                   :: sym,dX
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
  integer                    :: i
  integer                    :: imin, imax
  real*8                     :: dxm2
  real*8                     :: ONE,TWO
  parameter(ONE = 1.D0, TWO = 2.D0)
  imin = lbound(f,1)
  imax = ubound(f,1)
! Derivative with respect to x
  dxm2 = ONE/(dX*dX)
  d2f(imin,:,:) = dxm2 * (f(imin+1,:,:) + (sym - TWO)*f(imin,:,:))
  do i = imin + 1, imax - 1
     d2f(i,:,:) = dxm2 * (f(i+1,:,:) + f(i-1,:,:) - TWO*f(i,:,:))
  enddo
  d2f(imax,:,:) = d2f(imax-1,:,:)
end subroutine ddx2
!-----------------------------------------------------------------------------
!
! Computes partial derivative of function f with respect to y
!
!-----------------------------------------------------------------------------
subroutine ddy2(ext,f,d2f,dY,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,d2f
  real*8                                   :: sym,dY
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
  integer                    :: j
  integer                    :: jmin, jmax
  real*8                     :: dym2
  real*8                     :: ONE, TWO
  parameter(ONE = 1.0D0, TWO = 2.D0)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
! Derivative with respect to y
  dym2 = ONE/(dY*dY)
  d2f(:,jmin,:) = dym2 * (f(:,jmin+1,:) + (sym - TWO)*f(:,jmin,:))
  do j = jmin + 1, jmax - 1
     d2f(:,j,:) = dym2 * (f(:,j+1,:) + f(:,j-1,:) - TWO*f(:,j,:))
  enddo
  d2f(:,jmax,:) = d2f(:,jmax-1,:)
end subroutine ddy2
!-----------------------------------------------------------------------------
!
! Computes partial derivative of function f with respect to z
!
!-----------------------------------------------------------------------------
subroutine ddz2(ext,f,d2f,dZ,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,d2f
  real*8                                   :: sym,dZ
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
  integer                    :: k
  integer                    :: kmin, kmax
  real*8                     :: dzm1, dzm2
  real*8                     :: ONE, TWO
  parameter(ONE = 1.0D0, TWO = 2.D0)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
!
! Derivative with respect to z
!
  dzm2 = ONE/(dZ*dZ)
  d2f(:,:,kmin) = dzm2 * (f(:,:,kmin+1) + (sym - TWO)*f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     d2f(:,:,k) = dzm2 * (f(:,:,k+1) + f(:,:,k-1) - TWO*f(:,:,k))
  enddo
  d2f(:,:,kmax) = d2f(:,:,kmax-1)
end subroutine ddz2
!-----------------------------------------------------------------------------
!
! Computes second partial derivative of function f with respect to x. The 
! values of f are stored at the cell interfaces.
!
!-----------------------------------------------------------------------------
subroutine ddx2b(ext,f,d2f,dX,sym)
!
! Input parameters:
!
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,d2f
  real*8                                   :: sym,dX
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
  integer                    :: i
  integer                    :: imin, imax
  real*8                     :: dxm2
  real*8                     :: HALF, ONE, TWO
  parameter(ONE = 1.0D0, TWO = 2.D0)
  imin = lbound(f,1)
  imax = ubound(f,1)
! Derivative with respect to x
  dxm2 = ONE/(dX*dX)
  d2f(imin,:,:) = dxm2 * ((ONE + sym)*f(imin+1,:,:) - TWO*f(imin,:,:))
  do i = imin + 1, imax - 1
     d2f(i,:,:) = dxm2 * (f(i+1,:,:) + f(i-1,:,:) - TWO*f(i,:,:))
  enddo
  d2f(imax,:,:) = d2f(imax-1,:,:)
end subroutine ddx2b
!-----------------------------------------------------------------------------
!
! Computes partial derivative of function f with respect to y. The values of 
! f are stored at the cell interfaces.
!
!-----------------------------------------------------------------------------
subroutine ddy2b(ext,f,d2f,dY,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,d2f
  real*8                                   :: sym,dY
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
  integer                    :: j
  integer                    :: jmin, jmax
  real*8                     :: dym2
  real*8                     :: ONE, TWO
  parameter(ONE = 1.0D0, TWO = 2.D0)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
! Derivative with respect to y
  dym2 = ONE/(dY*dY)
  d2f(:,jmin,:) = dym2 * ((ONE+sym)*f(:,jmin+1,:) - TWO*f(:,jmin,:))
  do j = jmin + 1, jmax - 1
     d2f(:,j,:) = dym2 * (f(:,j+1,:) + f(:,j-1,:) - TWO*f(:,j,:))
  enddo
  d2f(:,jmax,:) = d2f(:,jmax-1,:)
end subroutine ddy2b
!-----------------------------------------------------------------------------
!
! Computes partial derivative of function f with respect to z. The values 
! of f are stored at the cell interfaces.
!
!-----------------------------------------------------------------------------
subroutine ddz2b(ext,f,d2f,dZ,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,d2f
  real*8                                   :: sym,dZ
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
  integer                    :: k
  integer                    :: kmin, kmax
  real*8                     :: dzm2
  real*8                     :: ONE, TWO
  parameter(ONE = 1.0D0, TWO = 2.D0)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
!
! Derivative with respect to z
!
  dzm2 = ONE/(dZ*dZ)
  d2f(:,:,kmin) = dzm2 * ((ONE+sym)*f(:,:,kmin+1) - TWO*f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     d2f(:,:,k) = dzm2 * (f(:,:,k+1) + f(:,:,k-1) - TWO*f(:,:,k))
  enddo
  d2f(:,:,kmax) = d2f(:,:,kmax-1)
end subroutine ddz2b
!
! Compute x-derivative of f to third-order,
! knowing f_f, a third-order approximate to f_{i-1/2}
subroutine ddx_o3(ext,f_f,df,dX,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f_f,df
  real*8                                   :: sym,dX
  integer                    :: i
  integer                    :: imin, imax
  real*8                     :: dxm1
  real*8                     :: HALF, ONE, TWO, GAM1
  real*8, dimension(ext(1),ext(2),ext(3))  :: f_temp
  parameter(ONE = 1.0D0, TWO = 2.D0, GAM1 = 1.D0/24.D0)
  imin = lbound(f_f,1)
  imax = ubound(f_f,1)
! Derivative with respect to x
  dxm1 = ONE/dX
!!$  call ddx2b(ext,f_f,f_temp,ONE,sym)
!!$  f_temp = f_f - GAM1 * f_temp
  f_temp = f_f
  do i=imin, imax-1
     df(i,:,:) = dxm1 * (f_temp(i+1,:,:) - f_temp(i,:,:))
  end do
  df(imax,:,:) = df(imax-1,:,:)
end subroutine ddx_o3
!
! Compute y-derivative of f to third-order,
! knowing f_f, a third-order approximate to f_{i-1/2}
subroutine ddy_o3(ext,f_f,df,dY,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f_f,df
  real*8                                   :: sym,dY
  integer                    :: j
  integer                    :: jmin, jmax
  real*8                     :: dym1
  real*8                     :: HALF, ONE, TWO, GAM1
  real*8, dimension(ext(1),ext(2),ext(3))  :: f_temp
  parameter(ONE = 1.0D0, TWO = 2.D0, GAM1 = 1.D0/24.D0)
  jmin = lbound(f_f,2)
  jmax = ubound(f_f,2)
! Derivative with respect to y
  dym1 = ONE/dY
!!$  call ddy2b(ext,f_f,f_temp,ONE,sym)
!!$  f_temp = f_f - GAM1 * f_temp
  f_temp = f_f
  do j=jmin, jmax-1
     df(:,j,:) = dym1 * (f_temp(:,j+1,:) - f_temp(:,j,:))
  end do
  df(:,jmax,:) = df(:,jmax-1,:)
end subroutine ddy_o3
!
! Compute z-derivative of f to third-order,
! knowing f_f, a third-order approximate to f_{i-1/2}
subroutine ddz_o3(ext,f_f,df,dZ,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f_f,df
  real*8                                   :: sym,dZ
  integer                    :: k
  integer                    :: kmin, kmax
  real*8                     :: dzm1
  real*8                     :: HALF, ONE, TWO, GAM1
  real*8, dimension(ext(1),ext(2),ext(3))  :: f_temp
  parameter(ONE = 1.0D0, TWO = 2.D0, GAM1 = 1.D0/24.D0)
  kmin = lbound(f_f,3)
  kmax = ubound(f_f,3)
! Derivative with respect to z
  dzm1 = ONE/dZ
!!$  call ddz2b(ext,f_f,f_temp,ONE,sym)
!!$  f_temp = f_f - GAM1 * f_temp
  f_temp = f_f
  do k=kmin, kmax-1
     df(:,:,k) = dzm1 * (f_temp(:,:,k+1) - f_temp(:,:,k))
  end do
  df(:,:,kmax) = df(:,:,kmax-1)
end subroutine ddz_o3
!-----------------------------------------------------------------------------
! Computes for x-advection
!-----------------------------------------------------------------------------
subroutine ddx3(ext,f,ddf,s,dX,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,ddf,s
  real*8                                   :: dX,sym
  integer                    :: i,j,k
  integer                    :: imin, imax, jmin, jmax, kmin, kmax
  real*8                     :: dxm1,dxm2,tdxm1
  real*8                     :: HALF, ONE, TWO, ZERO, EPSILON
!  parameter(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0, EPSILON = 1.0e-8)
  parameter(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0, EPSILON = 1.0D-8)
!
  imin = lbound(f,1)
  imax = ubound(f,1)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
!
  tdxm1 = ONE/dX
  dxm1 = tdxm1*HALF
  dxm2 = ONE/(dX*dX)
!
  do k = kmin, kmax
     do j = jmin, jmax
        ddf(imin,j,k) = dxm2 * (f(imin+1,j,k) + (sym - TWO)*f(imin,j,k))
        if((f(imin+1,j,k)-f(imin,j,k))*f(imin,j,k)*(ONE-sym) < ZERO) then
           s(imin,j,k) = ZERO
        else
           s(imin,j,k) = (TWO*(f(imin+1,j,k)-f(imin,j,k)) &
                * f(imin,j,k)*(ONE-sym) + EPSILON) &
          / ((f(imin+1,j,k)-f(imin,j,k))**2 + (f(imin,j,k)*(ONE-sym))**2 + EPSILON)
        end if
        do i = imin + 1, imax - 1
           ddf(i,j,k) = dxm2 * (f(i+1,j,k) - TWO*f(i,j,k) + f(i-1,j,k))
           if((f(i+1,j,k)-f(i,j,k))*(f(i,j,k)-f(i-1,j,k)) < ZERO) then
              s(i,j,k) = ZERO
           else
              s(i,j,k) = (TWO*(f(i+1,j,k)-f(i,j,k))*(f(i,j,k)-f(i-1,j,k)) + EPSILON) &
                   / ((f(i+1,j,k)-f(i,j,k))**2 + (f(i,j,k)-f(i-1,j,k))**2 + EPSILON)
           end if
        enddo
        ddf(imax,j,k) = ZERO
        s(imax,j,k) = ZERO
     end do
  end do
end subroutine ddx3
!-----------------------------------------------------------------------------
! Computes for y-advection
!-----------------------------------------------------------------------------
subroutine ddy3(ext,f,ddf,s,dY,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,ddf,s
  real*8                                   :: dY,sym
  integer                    :: i,j,k
  integer                    :: imin, imax, jmin, jmax, kmin, kmax
  real*8                     :: dym1,dym2,tdym1
  real*8                     :: HALF, ONE, TWO, ZERO, EPSILON
!  parameter(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0, EPSILON = 1.0e-8)
  parameter(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0, EPSILON = 1.0D-8)
!
  imin = lbound(f,1)
  imax = ubound(f,1)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
!
  tdym1 = ONE/dY
  dym1 = HALF*tdym1
  dym2 = ONE/(dY*dY)
!
  do k = kmin, kmax
     do i = imin, imax
        ddf(i,jmin,k) = dym2 * (f(i,jmin+1,k) + (sym - TWO)*f(i,jmin,k))
        if((f(i,jmin+1,k)-f(i,jmin,k))*f(i,jmin,k)*(ONE-sym) < ZERO) then
           s(i,jmin,k) = ZERO
        else
           s(i,jmin,k) = (TWO*(f(i,jmin+1,k)-f(i,jmin,k)) &
                *f(i,jmin,k)*(ONE-sym) + EPSILON) &
                / ((f(i,jmin+1,k)-f(i,jmin,k))**2 + (f(i,jmin,k)*(ONE-sym))**2 + EPSILON)
        end if
        do j = jmin + 1, jmax - 1
           ddf(i,j,k) = dym2 * (f(i,j+1,k) - TWO*f(i,j,k) + f(i,j-1,k))
           if((f(i,j+1,k)-f(i,j,k))*(f(i,j,k)-f(i,j-1,k)) < ZERO) then
              s(i,j,k) = ZERO
           else
              s(i,j,k) = (TWO*(f(i,j+1,k)-f(i,j,k))*(f(i,j,k)-f(i,j-1,k)) + EPSILON) &
                   / ((f(i,j+1,k)-f(i,j,k))**2 + (f(i,j,k)-f(i,j-1,k))**2 + EPSILON)
           end if
        enddo
        ddf(i,jmax,k) = ZERO
        s(i,jmax,k) = ZERO
     enddo
  end do
end subroutine ddy3
!-----------------------------------------------------------------------------
! Computes for z-advection
!-----------------------------------------------------------------------------
subroutine ddz3(ext,f,ddf,s,dZ,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,ddf,s
  real*8                                   :: dZ,sym
  integer                    :: i,j,k
  integer                    :: imin, imax, jmin, jmax, kmin, kmax
  real*8                     :: dzm1,dzm2,tdzm1
  real*8                     :: HALF, ONE, TWO, ZERO, EPSILON
!  parameter(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0, EPSILON = 1.0e-8)
  parameter(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0, EPSILON = 1.0D-8)
!
  imin = lbound(f,1)
  imax = ubound(f,1)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
!
  tdzm1 = ONE/dZ
  dzm1 = HALF*tdzm1
  dzm2 = ONE/(dZ*dZ)
!
  do j = jmin, jmax
     do i = imin, imax
        ddf(i,j,kmin) = dzm2 * (f(i,j,kmin+1) + (sym - TWO)*f(i,j,kmin))
        if((f(i,j,kmin+1)-f(i,j,kmin))*f(i,j,kmin)*(ONE-sym) < ZERO) then
           s(i,j,kmin) = ZERO
        else
           s(i,j,kmin) = (TWO*(f(i,j,kmin+1)-f(i,j,kmin)) &
                * f(i,j,kmin)*(ONE-sym) + EPSILON) &
                / ((f(i,j,kmin+1)-f(i,j,kmin))**2 + (f(i,j,kmin)*(ONE-sym))**2 + EPSILON)
        end if
        do k = kmin + 1, kmax - 1
           ddf(i,j,k) = dzm2 * (f(i,j,k+1) - TWO*f(i,j,k) + f(i,j,k-1))
           if((f(i,j,k+1)-f(i,j,k))*(f(i,j,k)-f(i,j,k-1)) < ZERO) then
              s(i,j,k) = ZERO
           else
              s(i,j,k) = (TWO*(f(i,j,k+1)-f(i,j,k))*(f(i,j,k)-f(i,j,k-1)) + EPSILON) &
                   / ((f(i,j,k+1)-f(i,j,k))**2 + (f(i,j,k)-f(i,j,k-1))**2 + EPSILON)
           end if
        enddo
        ddf(i,j,kmax) = ZERO
        s(i,j,kmax) = ZERO
     enddo
  end do
end subroutine ddz3
