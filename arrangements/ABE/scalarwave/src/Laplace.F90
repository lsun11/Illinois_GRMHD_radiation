!-----------------------------------------------------------------------------
!
! $Id: Laplace.F90,v 1.1.1.1 2006/02/17 21:31:26 zetienne Exp $
!
!-----------------------------------------------------------------------------
!
! Laplace Operator of Function f (for interior, and inner boundaries)
!
!-----------------------------------------------------------------------------
subroutine Laplace(f,DDf,ext,X,Y,Z,dX,dY,dZ,xsym,ysym,zsym) 
!
! Input parameters:
!
  integer, dimension(3)                        :: ext
  real*8                                       :: dX, dY, dZ
  real*8, dimension(ext(1),ext(2),ext(3))      :: f,DDf,X,Y,Z
  integer                                      :: i,j,k
  integer                                      :: xsym,ysym,zsym
!
! sym give symmetries across coordinate planes: 
!     1 :   symmetric
!    -1 :   antisymmetric
! 
! Other variables
!
  integer                    :: imin, jmin, kmin, imax, jmax, kmax  
  real*8                     :: dxm2, dym2, dzm2
  real*8                     :: HALF, ONE, TWO
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
!
! sanity check:
!
  if (abs(xsym) /= ONE) then
     write(*,*) "OOOPS! ", xsym, "is no valid symmetry!"
     stop
  endif
  if (abs(ysym) /= ONE) then
     write(*,*) "OOOPS! ", ysym, "is no valid symmetry!"
     stop
  endif
  if (abs(zsym) /= ONE) then
     write(*,*) "OOOPS! ", zsym, "is no valid symmetry!"
     stop
  endif
!
! Input translation
!


!  if(x(lbound(f,1),1,1) < 0) then
!     imin = lbound(f,1)+1
!  else 
!     imin = lbound(f,1)
!  endif
!  if(y(1,lbound(f,2),1) < 0) then
!     jmin = lbound(f,2)+1
!  else 
!     jmin = lbound(f,2)
!  endif
!  if(z(1,1,lbound(f,3)) < 0) then
!     kmin = lbound(f,3)+1
!  else
!     kmin = lbound(f,3)
!  endif
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  imax = ubound(f,1)
  jmax = ubound(f,2)
  kmax = ubound(f,3)
  dxm2 = ONE/(dX*dX)
  dym2 = ONE/(dY*dY)
  dzm2 = ONE/(dZ*dZ)
!
! x = 0, y = 0, z = 0
!
  i = imin
  j = jmin
  k = kmin
 ! write(*,*) imin, jmin, kmin, i,j,k,xsym,f(i,j,k)
  DDf(i,j,k) = & 
       dxm2*(xsym*f(i,j,k) - TWO*f(i,j,k) + f(i+1,j,k)) + &
       dym2*(ysym*f(i,j,k) - TWO*f(i,j,k) + f(i,j+1,k)) + &
       dzm2*(zsym*f(i,j,k) - TWO*f(i,j,k) + f(i,j,k+1))

!
! x = 0, y = 0, z = int  
!
  i = imin
  j = jmin
  do k = kmin + 1, kmax - 1
     DDf(i,j,k) = &
          dxm2*(xsym*f(i,j,k) - TWO*f(i,j,k) + f(i+1,j,k)) + &
          dym2*(ysym*f(i,j,k) - TWO*f(i,j,k) + f(i,j+1,k)) + &
          dzm2*(f(i,j,k-1) - TWO*f(i,j,k) + f(i,j,k+1))
  enddo
!
! x = 0, y = int, z = 0
!
  i = imin 
  k = kmin
  do j = jmin + 1, jmax - 1
     DDf(i,j,k) = &
          dxm2*(xsym*f(i,j,k) - TWO*f(i,j,k) + f(i+1,j,k)) + &
          dym2*(f(i,j-1,k) - TWO*f(i,j,k) + f(i,j+1,k)) + &
          dzm2*(zsym*f(i,j,k) - TWO*f(i,j,k) + f(i,j,k+1))
  enddo
!
! x = int, y = 0, z = 0
!
  j = jmin
  k = kmin
  do i = imin + 1, imax - 1 
     DDf(i,j,k) = &
          dxm2*(f(i-1,j,k) - TWO*f(i,j,k) + f(i+1,j,k)) + &
          dym2*(ysym*f(i,j,k) - TWO*f(i,j,k) + f(i,j+1,k)) + &
          dzm2*(zsym*f(i,j,k) - TWO*f(i,j,k) + f(i,j,k+1))
  enddo
!
! x = 0, y = int, z = int
!
  i = imin
  do k = kmin + 1, kmax - 1
     do j = jmin + 1, jmax - 1
        DDf(i,j,k) = &
             dxm2*(xsym*f(i,j,k) - TWO*f(i,j,k) + f(i+1,j,k)) + &
             dym2*(f(i,j-1,k) - TWO*f(i,j,k) + f(i,j+1,k)) + &
             dzm2*(f(i,j,k-1) - TWO*f(i,j,k) + f(i,j,k+1))
     enddo
  enddo
!
! x = int, y = 0, z = int
!
  j = jmin
  do k = kmin + 1, kmax - 1
     do i = imin + 1, imax - 1
        DDf(i,j,k) = &
             dxm2*(f(i-1,j,k) - TWO*f(i,j,k) + f(i+1,j,k)) + &
             dym2*(ysym*f(i,j,k) - TWO*f(i,j,k) + f(i,j+1,k)) + &
             dzm2*(f(i,j,k-1) - TWO*f(i,j,k) + f(i,j,k+1))
     enddo
  enddo
!
! x = int, y = int, z = 0
!
  k = kmin
  do j = jmin + 1, jmax - 1
     do i = imin + 1, imax - 1
        DDf(i,j,k) = &
             dxm2*(f(i-1,j,k) - TWO*f(i,j,k) + f(i+1,j,k)) + &
             dym2*(f(i,j-1,k) - TWO*f(i,j,k) + f(i,j+1,k)) + &
             dzm2*(zsym*f(i,j,k) - TWO*f(i,j,k) + f(i,j,k+1))
     enddo
  enddo
!
! x = int, y = int, z = int
!
do k = kmin + 1, kmax - 1
  do j = jmin + 1, jmax - 1
     do i = imin + 1, imax - 1
           DDf(i,j,k) = &
                dxm2*(f(i-1,j,k) - TWO*f(i,j,k) + f(i+1,j,k)) + &
                dym2*(f(i,j-1,k) - TWO*f(i,j,k) + f(i,j+1,k)) + &
                dzm2*(f(i,j,k-1) - TWO*f(i,j,k) + f(i,j,k+1))
        enddo
     enddo
  enddo
  return 
end subroutine Laplace

