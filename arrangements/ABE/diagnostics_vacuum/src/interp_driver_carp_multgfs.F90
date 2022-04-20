#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!------------------------------------
! Do stationary puncture diagnostics
!------------------------------------
subroutine interp_driver_carp_multgfs(cctkGH,num_points,num_gfs,pointcoords,gridfunc_varindices,output)
  implicit none
  DECLARE_CCTK_FUNCTIONS

  CCTK_POINTER                             :: cctkGH
  CCTK_POINTER, dimension(3)               :: interp_coords

  integer                                  :: num_points,num_gfs
  real*8, dimension(num_points,3)          :: pointcoords

  character(60)                            :: options_string,interpolation_operator
  character(128)                           :: operator
  integer                                  :: interpolation_order,i
  integer                                  :: nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle
  integer                                  :: N_interp_points
  integer                                  :: N_input_arrays,N_output_arrays,input_array_type_codes,output_array_type_codes
  real*8, dimension(num_points)            :: xinterp,yinterp,zinterp
  integer,dimension(num_gfs)               :: gridfunc_varindices,input_array_indices
  CCTK_POINTER,dimension(num_gfs)          :: output_array_pointers
  real*8,dimension(num_points)             :: output_array1,output_array2,output_array3,output_array4,output_array5
  real*8,dimension(num_gfs,num_points)     :: output

  if(num_gfs.gt.5 .or. num_gfs.lt.1) then
     write(*,*) "Error: interp_driver_carp_multgfs() only supports 1 <= N <= 5 gridfunctions at once."
     stop
  end if

  N_dims = 3
  interpolation_order = 4
!  interpolation_order = 3
  N_interp_points = num_points
  N_input_arrays = num_gfs
  N_output_arrays = num_gfs

  input_array_indices = gridfunc_varindices

  param_table_handle = -1
  interp_handle = -1
  coord_system_handle = -1

!  call CCTK_InterpHandle (interp_handle, "generalized polynomial interpolation")
  call CCTK_InterpHandle (interp_handle,  "Lagrange polynomial interpolation")
!  call CCTK_InterpHandle (interp_handle, "uniform cartesian")
  if (interp_handle .lt. 0) then
     call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
  endif
  
  options_string = "order = " // char(ichar('0') + interpolation_order)
  call Util_TableCreateFromString (param_table_handle, options_string)
  if (param_table_handle .lt. 0) then
     call CCTK_WARN(0,"Cannot create parameter table for interpolator")
  endif

  call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
  if (coord_system_handle .lt. 0) then
     call CCTK_WARN(0,"Cannot get handle for cart3d coordinate system ! Forgot to activate an implementation providing coordinates ??")
  endif

  input_array_type_codes = CCTK_VARIABLE_REAL
  output_array_type_codes = CCTK_VARIABLE_REAL

  xinterp = pointcoords(:,1)
  yinterp = pointcoords(:,2)
  zinterp = pointcoords(:,3)

  !Zach says: Sorry about the ugliness below.  
  ! Feel free to fix it if you can figure out how to reduce the number of lines!
  if(num_gfs==1) then
     output_array_pointers(1) = CCTK_PointerTo(output_array1)
  else if(num_gfs==2) then
     output_array_pointers(1) = CCTK_PointerTo(output_array1)
     output_array_pointers(2) = CCTK_PointerTo(output_array2)
  else if(num_gfs==3) then
     output_array_pointers(1) = CCTK_PointerTo(output_array1)
     output_array_pointers(2) = CCTK_PointerTo(output_array2)
     output_array_pointers(3) = CCTK_PointerTo(output_array3)
  else if(num_gfs==4) then
     output_array_pointers(1) = CCTK_PointerTo(output_array1)
     output_array_pointers(2) = CCTK_PointerTo(output_array2)
     output_array_pointers(3) = CCTK_PointerTo(output_array3)
     output_array_pointers(4) = CCTK_PointerTo(output_array4)
  else if(num_gfs==5) then
     output_array_pointers(1) = CCTK_PointerTo(output_array1)
     output_array_pointers(2) = CCTK_PointerTo(output_array2)
     output_array_pointers(3) = CCTK_PointerTo(output_array3)
     output_array_pointers(4) = CCTK_PointerTo(output_array4)
     output_array_pointers(5) = CCTK_PointerTo(output_array5)
  end if

  interp_coords(1) = CCTK_PointerTo(xinterp)
  interp_coords(2) = CCTK_PointerTo(yinterp)
  interp_coords(3) = CCTK_PointerTo(zinterp)

  call CCTK_InterpGridArrays(ierr,cctkGH,N_dims,interp_handle, &
       param_table_handle,coord_system_handle, &
       N_interp_points,input_array_type_codes,interp_coords, &
       N_input_arrays,input_array_indices, &
       N_output_arrays, output_array_type_codes, output_array_pointers)

  if(num_gfs==1) then
     output(1,:) = output_array1
  else if(num_gfs==2) then
     output(1,:) = output_array1
     output(2,:) = output_array2
  else if(num_gfs==3) then
     output(1,:) = output_array1
     output(2,:) = output_array2
     output(3,:) = output_array3
  else if(num_gfs==4) then
     output(1,:) = output_array1
     output(2,:) = output_array2
     output(3,:) = output_array3
     output(4,:) = output_array4
  else if(num_gfs==5) then
     output(1,:) = output_array1
     output(2,:) = output_array2
     output(3,:) = output_array3
     output(4,:) = output_array4
     output(5,:) = output_array5
  end if

end subroutine interp_driver_carp_multgfs
