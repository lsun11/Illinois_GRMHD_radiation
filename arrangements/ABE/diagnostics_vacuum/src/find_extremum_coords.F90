#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!-------------------------------------------------------------------------------------------
! Find the global spatial coordinates corresponding to a value for specified gridfunction
!
! Don't expect this function to work if the exact value occurs in > 1 place on the grid!
!-------------------------------------------------------------------------------------------

subroutine find_gfvalue_coords(cctkGH,ext,variable_index,variable_gf,X,Y,Z,gfvalue,x_coord,y_coord,z_coord)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_POINTER                             :: cctkGH
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: variable_gf,X,Y,Z
  integer                                  :: variable_index
  real*8                                   :: x_coord,y_coord,z_coord,gfvalue

  real*8                                   :: gfvalue_temp, TOL
  integer                                  :: i,j,k
  integer                                  :: sum_handle,ierr

  ! TOL = tolerance.
  TOL = 1.D-8

  gfvalue_temp = 1.D-200
  x_coord = 0.D0
  y_coord = 0.D0
  z_coord = 0.D0
  
  !----------------------------------------------------------------
  !First figure out which point contains the extremum value.
  !We search for the extremum point everywhere, even the ghostzones
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           if(abs(variable_gf(i,j,k)-gfvalue) < TOL) then
              gfvalue_temp = variable_gf(i,j,k)
              x_coord = X(i,j,k)
              y_coord = Y(i,j,k)
              z_coord = Z(i,j,k)
           end if
        end do
     end do
  end do

  call CCTK_Reductionhandle(sum_handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,sum_handle,gfvalue_temp,gfvalue_temp,CCTK_VARIABLE_REAL)
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,sum_handle,x_coord,x_coord,CCTK_VARIABLE_REAL)
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,sum_handle,y_coord,y_coord,CCTK_VARIABLE_REAL)
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,sum_handle,z_coord,z_coord,CCTK_VARIABLE_REAL)

  !Now correct for any ghostzone overlap
  x_coord = x_coord * gfvalue/gfvalue_temp
  y_coord = y_coord * gfvalue/gfvalue_temp
  z_coord = z_coord * gfvalue/gfvalue_temp
  
  write(*,*) "COORDINATES CORRESPONDING TO FUNCTION MAXIMUM:",x_coord,y_coord,z_coord,gfvalue,gfvalue_temp

end subroutine find_gfvalue_coords
