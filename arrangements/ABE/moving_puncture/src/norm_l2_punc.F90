!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine norm2_dv_punc(ex,f1,X,Y,Z,norm,adjimin,adjimax,adjjmin,adjjmax, &
     adjkmin,adjkmax,Symmetry,xbh1,ybh1,zbh1)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                 :: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: f1
  real*8, dimension(ex(1),ex(2),ex(3))  :: X,Y,Z
  real*8                                :: norm,xbh1,ybh1,zbh1
  integer                               :: adjimin,adjimax
  integer                               :: adjjmin,adjjmax
  integer                               :: adjkmin,adjkmax
  integer                               :: Symmetry
!
  integer         :: imin, jmin, kmin, imax, jmax, kmax
  integer         :: i,j,k
  real*8          :: dX,dY,dZ,dV,xi,yj,zk,rr
  real*8          :: ZERO, TWO,PI
  integer         :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter ( NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3)
  parameter ( AXISYM = 4 )
  parameter ( ZERO = 0.D0, TWO = 2.D0 )
!
! Input translation
!
  PI = acos(-1.D0)
  norm = ZERO
  imin = lbound(f1,1) - adjimin
  jmin = lbound(f1,2) - adjjmin
  kmin = lbound(f1,3) - adjkmin
  imax = ubound(f1,1) - adjimax
  jmax = ubound(f1,2) - adjjmax
  kmax = ubound(f1,3) - adjkmax
  dX = X(imin + 1,1,1) - X(imin,1,1)
  dY = Y(1,jmin + 1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin + 1) - Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
  end if

  do k = kmin, kmax
     zk=z(1,1,k)-zbh1
     do j = jmin, jmax
        yj=y(1,j,1)-ybh1
        do i = imin, imax
           xi=x(i,1,1)-xbh1
           rr=sqrt(xi**2+yj**2+zk**2)
           if(rr > 0.5) then
              if (Symmetry==AXISYM) then
                 norm = norm + f1(i,j,k)*f1(i,j,k)*X(i,1,1)
              else
                 norm = norm + f1(i,j,k)*f1(i,j,k)
              end if
           endif
        enddo
     enddo
  enddo

  if (Symmetry.eq.PI_SYMM) norm =  TWO * norm

  norm = norm * dV * dV

  if (Symmetry==AXISYM) norm=norm*4.d0*PI/dY

  return

end subroutine norm2_dv_punc
