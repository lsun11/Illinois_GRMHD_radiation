!------------------------------------
! Do stationary puncture diagnostics
!------------------------------------
subroutine interp_driver_carp_multgfs(cctkGH,num_points,num_gfs,pointcoords,gridfunc_varindices,output)
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
  INTEGER*4 function HorizonCentroid (horizon_number, centroid_x, centroid_y, centroid_z)
  implicit none
  INTEGER*4 horizon_number
  REAL*8 centroid_x
  REAL*8 centroid_y
  REAL*8 centroid_z
  end function HorizonCentroid
  end interface
  interface
  INTEGER*4 function HorizonLocalCoordinateOrigin (horizon_number, origin_x, origin_y, origin_z)
  implicit none
  INTEGER*4 horizon_number
  REAL*8 origin_x
  REAL*8 origin_y
  REAL*8 origin_z
  end function HorizonLocalCoordinateOrigin
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
  interface
  INTEGER*4 function HorizonWasFound (horizon_number)
  implicit none
  INTEGER*4 horizon_number
  end function HorizonWasFound
  end interface
  
  integer*8                             :: cctkGH
  integer*8, dimension(3)               :: interp_coords
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
  integer*8,dimension(num_gfs)          :: output_array_pointers
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
!  call CCTK_InterpHandle (interp_handle, generalized polynomial interpolation)
  call CCTK_InterpHandle (interp_handle,  "Lagrange polynomial interpolation")
!  call CCTK_InterpHandle (interp_handle, uniform cartesian)
  if (interp_handle .lt. 0) then
     call CCTK_Warn(0,54,"interp_driver_carp_multgfs.F90","diagnostics_vacuum","Cannot get handle for interpolation ! Forgot to act&
  &ivate an implementation providing interpolation operators ??")
  endif
  options_string = "order = " // char(ichar('0') + interpolation_order)
  call Util_TableCreateFromString (param_table_handle, options_string)
  if (param_table_handle .lt. 0) then
     call CCTK_Warn(0,60,"interp_driver_carp_multgfs.F90","diagnostics_vacuum","Cannot create parameter table for interpolator")
  endif
  call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
  if (coord_system_handle .lt. 0) then
     call CCTK_Warn(0,65,"interp_driver_carp_multgfs.F90","diagnostics_vacuum","Cannot get handle for cart3d coordinate system ! Fo&
  &rgot to activate an implementation providing coordinates ??")
  endif
  input_array_type_codes = 107
  output_array_type_codes = 107
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
