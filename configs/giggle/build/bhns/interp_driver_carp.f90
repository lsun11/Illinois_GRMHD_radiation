!------------------------------------
! Do stationary puncture diagnostics
!------------------------------------
subroutine bhns_interp_driver_carp_arborder(cctkGH,num_points,pointcoords,gridfunc_varindex,output,interpolation_order)
  implicit none
  external     CCTK_PointerTo
  integer*8 CCTK_PointerTo
  interface
  integer function CCTK_Equals (arg1, arg2)
  implicit none
  integer*8 arg1
  character(*) arg2
  end function CCTK_Equals
  integer function CCTK_MyProc (cctkGH)
  implicit none
  integer*8 cctkGH
  end function CCTK_MyProc
  integer function CCTK_nProcs (cctkGH)
  implicit none
  integer*8 cctkGH
  end function CCTK_nProcs
  integer function CCTK_IsThornActive (name)
  implicit none
  character(*) name
  end function CCTK_IsThornActive
  integer*8 function CCTK_NullPointer ()
  implicit none
  end function CCTK_NullPointer
  end interface
  interface
  INTEGER*4 function HorizonRadiusInDirection (horizon_number, N_points, x, y, z, radius)
  implicit none
  INTEGER*4 horizon_number
  INTEGER*4 N_points
  REAL*8 x(*)
  REAL*8 y(*)
  REAL*8 z(*)
  REAL*8 radius(*)
  end function HorizonRadiusInDirection
  end interface
  
  integer*8                             :: cctkGH
  integer*8, dimension(3)               :: interp_coords
  integer                                  :: num_points
  real*8, dimension(num_points,3)          :: pointcoords
  character*9                            :: options_string
!  character*128                           :: operator
  integer                                  :: interpolation_order
  integer                                  :: nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle
  integer                                  :: N_interp_points
  integer                                  :: N_input_arrays,N_output_arrays,input_array_type_codes,output_array_type_codes
  real*8, dimension(num_points)            :: xinterp,yinterp,zinterp
  integer                                  :: gridfunc_varindex
  integer,dimension(1)                     :: input_array_indices
  integer*8,dimension(1)                :: output_array_pointers
! real*8,dimension(num_points)             :: output_arrays
  real*8,dimension(num_points)             :: output
  N_dims = 3
!  interpolation_order = 4
!  interpolation_order = 3
  N_interp_points = num_points
  N_input_arrays = 1
  N_output_arrays = 1
  input_array_indices(1) = gridfunc_varindex
  param_table_handle = -1
  interp_handle = -1
  coord_system_handle = -1
!  call CCTK_InterpHandle (interp_handle, generalized polynomial interpolation)
  call CCTK_InterpHandle (interp_handle,  "Lagrange polynomial interpolation")
!  call CCTK_InterpHandle (interp_handle, uniform cartesian)
  if (interp_handle .lt. 0) then
     call CCTK_Warn(0,50,"interp_driver_carp.F90","bhns","Cannot get handle for interpolation ! Forgot to activate an implementatio&
  &n providing interpolation operators ??")
  endif
  options_string = "order = " // char(ichar('0') + interpolation_order)
  call Util_TableCreateFromString (param_table_handle, options_string)
  if (param_table_handle .lt. 0) then
     call CCTK_Warn(0,56,"interp_driver_carp.F90","bhns","Cannot create parameter table for interpolator")
  endif
  call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
  if (coord_system_handle .lt. 0) then
     call CCTK_Warn(0,61,"interp_driver_carp.F90","bhns","Cannot get handle for cart3d coordinate system ! Forgot to activate an im&
  &plementation providing coordinates ??")
  endif
  input_array_type_codes = 107
  output_array_type_codes = 107
  xinterp = pointcoords(:,1)
  yinterp = pointcoords(:,2)
  zinterp = pointcoords(:,3)
!  output_array_pointers(1) = CCTK_PointerTo(output_arrays)
  output_array_pointers(1) = CCTK_PointerTo(output)
  interp_coords(1) = CCTK_PointerTo(xinterp)
  interp_coords(2) = CCTK_PointerTo(yinterp)
  interp_coords(3) = CCTK_PointerTo(zinterp)
  call CCTK_InterpGridArrays(ierr,cctkGH,N_dims,interp_handle, &
       param_table_handle,coord_system_handle, &
       N_interp_points,input_array_type_codes,interp_coords, &
       N_input_arrays,input_array_indices, &
       N_output_arrays, output_array_type_codes, output_array_pointers)
!  output=output_arrays
end subroutine bhns_interp_driver_carp_arborder
