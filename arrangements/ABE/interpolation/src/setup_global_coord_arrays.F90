#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine setup_global_coord_arrays(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3)                    :: ext
  real*8                                   :: dT,dX,dY,dZ,xoffset,yoffset,zoffset
  real*8                                   :: xmin,ymin,zmin,xmax,ymax,zmax
  integer                                  :: i,j,k
  integer                                  :: handle,index,ierr
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymax,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax,1,index)
  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmin,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymin,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmin,1,index)
  
  do i=1,cctk_gsh(1)
     Xglobal(i) = xmin + (i-1.D0)*dX
  end do
  do j=1,cctk_gsh(2)
     Yglobal(j) = ymin + (j-1.D0)*dY
  end do
  do k=1,cctk_gsh(3)
     Zglobal(k) = zmin + (k-1.D0)*dZ
  end do
  
  do i=1,cctk_lsh(1)
     Xlocal1d(i) = X(i,1,1)
  end do
  do j=1,cctk_lsh(2)
     Ylocal1d(j) = Y(1,j,1)
  end do
  do k=1,cctk_lsh(3)
     Zlocal1d(k) = Z(1,1,k)
  end do


!  write(*,*) "blah1",Xlocal1d(12),cctk_lsh(1)
!  write(*,*) "blah2",Ylocal1d(2),cctk_lsh(2)
!  write(*,*) "blah3",Zlocal1d(12),cctk_lsh(3)

!  if(Symmetry==AXISYM) then
!     Yglobal(1) = -dY
!     Yglobal(2) = 0.D0
!     Yglobal(3) = dY
!  end if

!  do i=1,cctk_gsh(1)
!     write(*,*) "X:",Xglobal(i)     
!  end do
!  do i=1,cctk_gsh(2)
!     write(*,*) "Y:",Yglobal(i)     
!  end do
!  do i=1,cctk_gsh(3)
!     write(*,*) "Z:",Zglobal(i)     
!  end do

end subroutine setup_global_coord_arrays
