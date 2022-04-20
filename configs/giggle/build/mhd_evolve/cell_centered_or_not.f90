!-----------------------------------------------------------------------------
! Determines whether our local grid is cell-centered or not, for symmetry
!  ghostzone filling purposes.
! 
!  Note that this test will yield cell_centering_enabled==1 if the local
!  processor grid does not touch a X=0, Y=0, or Z=0 boundary.  This is not
!  a problem since there will be no symmetry ghostzones on such grids!
!-----------------------------------------------------------------------------
subroutine cell_centered_or_not(ext,X,Y,Z,cell_centering_enabled)
  implicit none
  !INPUT PARAMETERS: 
  integer,dimension(3)                     :: ext
  real*8,dimension(ext(1),ext(2),ext(3))   :: X,Y,Z
  !OUTPUT PARAMETERS: 
  integer                                  :: cell_centering_enabled
  !OTHER PARAMETERS: 
  integer                                  :: i,j,k
  real*8                                   :: TOL
  TOL = (Z(1,1,2) - Z(1,1,1)) * 1.D-3
  cell_centering_enabled = 1
  do k=1,ext(3)
     if(abs(Z(1,1,k)).lt.TOL) cell_centering_enabled=0
  end do
  do j=1,ext(2)
     if(abs(Y(1,j,1)).lt.TOL) cell_centering_enabled=0
  end do
  do i=1,ext(1)
     if(abs(X(i,1,1)).lt.TOL) cell_centering_enabled=0
  end do
end subroutine cell_centered_or_not
