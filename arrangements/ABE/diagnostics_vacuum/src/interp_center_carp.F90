#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!------------------------------------
! Do stationary puncture diagnostics
!------------------------------------
subroutine interp_center_carp(cctkGH,gridfunc_varindex,output)
  implicit none
  DECLARE_CCTK_FUNCTIONS

  CCTK_POINTER                             :: cctkGH
  CCTK_POINTER, dimension(3)               :: interp_coords
  character(60)                            :: options_string,interpolation_operator
  character(128)                           :: operator
  integer                                  :: interpolation_order
  integer                                  :: nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle
  integer                                  :: N_interp_points
  integer                                  :: N_input_arrays,N_output_arrays,input_array_type_codes,output_array_type_codes
  real*8, dimension(1)                     :: xinterp,yinterp,zinterp
  integer                                  :: gridfunc_varindex
  integer,dimension(1)                     :: input_array_indices
  CCTK_POINTER,dimension(1)                :: output_array_pointers
  real*8,dimension(1)                      :: output_arrays
  real*8                                   :: output

  N_dims = 3
  interpolation_order = 3
  N_interp_points = 1
  N_input_arrays = 1
  N_output_arrays = 1

  input_array_indices(1) = gridfunc_varindex

  param_table_handle = -1
  interp_handle = -1
  coord_system_handle = -1

!  call CCTK_InterpHandle (interp_handle, "generalized polynomial interpolation")
  call CCTK_InterpHandle (interp_handle, "Lagrange polynomial interpolation")
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

  xinterp = 0.D0
  yinterp = 0.D0
  zinterp = 0.D0

  output_array_pointers(1) = CCTK_PointerTo(output_arrays)

  interp_coords(1) = CCTK_PointerTo(xinterp)
  interp_coords(2) = CCTK_PointerTo(yinterp)
  interp_coords(3) = CCTK_PointerTo(zinterp)

  call CCTK_InterpGridArrays(ierr,cctkGH,N_dims,interp_handle, &
       param_table_handle,coord_system_handle, &
       N_interp_points,input_array_type_codes,interp_coords, &
       N_input_arrays,input_array_indices, &
       N_output_arrays, output_array_type_codes, output_array_pointers)

  output=output_arrays(1)

end subroutine interp_center_carp
