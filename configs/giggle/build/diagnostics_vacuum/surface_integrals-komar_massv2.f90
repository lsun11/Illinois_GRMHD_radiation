!-----------------------------------------------------------------------------
! Compute Komar mass via surface integral
!-----------------------------------------------------------------------------
subroutine Komar_surf_integralv2(cctkGH,output_integral)
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
  
  REAL*8  excision_radius
  REAL*8  run_time
  INTEGER*4 Symmetry
  INTEGER*4 bssn_enable
  INTEGER*4 cowling_enable
  INTEGER*4 excision_enable
  INTEGER*4 fisheye_enable
  INTEGER*4 iter_count
  INTEGER*4 number_of_mol_ministeps
  INTEGER*4 rot_metric
  INTEGER*4 trA_detg_enforce
  COMMON /cctk_params_global/excision_radius,run_time,Symmetry,bssn_enable,cowling_enable,excision_enable,fisheye_enable,iter_count&
  &,number_of_mol_ministeps,rot_metric,trA_detg_enforce
  REAL*8  BH_Vol_Excise_Radius
  REAL*8  dcostheta
  REAL*8  ddrbddr
  REAL*8  dphi
  REAL*8  drbdr
  REAL*8  inner_volInt_radius
  REAL*8  rbr
  REAL*8  rsurf2
  REAL*8  surf_radius
  INTEGER*4 Compute_VolIntegrands_Every
  INTEGER*4 N_phi
  INTEGER*4 N_theta
  INTEGER*4 WhichIntegral
  INTEGER*4 arithsurf
  INTEGER*4 enable_Jz_constraint
  INTEGER*4 enable_M_constraint
  INTEGER*4 enable_P_constraint
  INTEGER*4 nsurf
  INTEGER*4 ntot
  INTEGER*4 num_BHs
  INTEGER*4 numphi
  INTEGER*4 numtheta
  INTEGER*4 scaledsurf
  INTEGER*4 sym_factor
  COMMON /diagnostics_vacuumrest/BH_Vol_Excise_Radius,dcostheta,ddrbddr,dphi,drbdr,inner_volInt_radius,rbr,rsurf2,surf_radius,Compu&
  &te_VolIntegrands_Every,N_phi,N_theta,WhichIntegral,arithsurf,enable_Jz_constraint,enable_M_constraint,enable_P_constraint,nsurf,&
  &ntot,num_BHs,numphi,numtheta,scaledsurf,sym_factor
  REAL*8  CCTKH0
  REAL*8  CCTKH1
  REAL*8  CCTKH2
  REAL*8  CCTKH3
  REAL*8  CCTKH13
  REAL*8  CCTKH14
  REAL*8  CCTKH15
  REAL*8  CCTKH16
  REAL*8  CCTKH17
  REAL*8  CCTKH18
  REAL*8  CCTKH19
  REAL*8  CCTKH20
  integer*8  bitant_plane
  integer*8  domain
  integer*8  CCTKH4
  integer*8  CCTKH5
  integer*8  CCTKH12
  INTEGER*4 CCTKH6
  INTEGER*4 CCTKH7
  INTEGER*4 CCTKH8
  INTEGER*4 CCTKH9
  INTEGER*4 CCTKH10
  INTEGER*4 CCTKH11
  COMMON /GRIDrest/CCTKH0,CCTKH1,CCTKH2,CCTKH3,CCTKH13,CCTKH14,CCTKH15,CCTKH16,CCTKH17,CCTKH18,CCTKH19,CCTKH20,bitant_plane,domain,&
  &CCTKH4,CCTKH5,CCTKH12,CCTKH6,CCTKH7,CCTKH8,CCTKH9,CCTKH10,CCTKH11
  REAL*8  CCTKH23
  REAL*8  CCTKH24
  REAL*8  CCTKH25
  integer*8  slicing_type
  INTEGER*4 CCTKH21
  INTEGER*4 CCTKH22
  INTEGER*4 CCTKH26
  INTEGER*4 CCTKH27
  INTEGER*4 CCTKH28
  INTEGER*4 CCTKH29
  INTEGER*4 CCTKH30
  COMMON /LAPSErest/CCTKH23,CCTKH24,CCTKH25,slicing_type,CCTKH21,CCTKH22,CCTKH26,CCTKH27,CCTKH28,CCTKH29,CCTKH30
  INTEGER*4 Spatial_Gauge
  INTEGER*4 CCTKH31
  INTEGER*4 CCTKH32
  COMMON /SHIFTrest/Spatial_Gauge,CCTKH31,CCTKH32
  
  ! Variables in the function call
  integer*8                             :: cctkGH
  real*8                                   :: output_integral
  ! Variables needed for interpolation
  integer*8, dimension(3)               :: interp_coords
  character(60)                            :: options_string
  integer                                  :: interpolation_order,nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle,N_interp_points
  integer                                  :: N_input_arrays,N_output_arrays
  real*8, dimension(N_theta*N_phi)         :: xinterp,yinterp,zinterp
  integer,dimension(26)                    :: input_array_type_codes,input_array_varindices,output_array_type_codes
  integer*8,dimension(26)               :: output_array_pointers
  ! Output arrays, dummy indices, parameters
  real*8,dimension(N_theta*N_phi)          :: Psi4int,gupxxint,gupxyint,gupxzint,gupyyint,gupyzint,gupzzint
  real*8,dimension(N_theta*N_phi)          :: gxxint,gxyint,gxzint,gyyint,gyzint,gzzint
  real*8,dimension(N_theta*N_phi)          :: Axxint,Axyint,Axzint,Ayyint,Ayzint,Azzint,trKint
  real*8,dimension(N_theta*N_phi)          :: lapsexint,lapseyint,lapsezint
  real*8,dimension(N_theta*N_phi)          :: shiftxint,shiftyint,shiftzint
  real*8                                   :: phiangle,costheta,sintheta,PI,nxp,nyp,nzp
  integer                                  :: i,j,n
  PI = 3.14159265358979323846D0
  N_dims = 3
  interpolation_order = 2
  N_interp_points = N_theta*N_phi
  interp_handle = -1
  call CCTK_InterpHandle (interp_handle, "Lagrange polynomial interpolation")
!  call CCTK_InterpHandle (interp_handle, uniform cartesian)
  if (interp_handle .lt. 0) then
     call CCTK_Warn(0,45,"surface_integrals-komar_massv2.F90","diagnostics_vacuum","Cannot get handle for interpolation ! Forgot to&
  & activate an implementation providing interpolation operators ??")
  endif
  param_table_handle = -1
  options_string = "order = " // char(ichar('0') + interpolation_order)
  call Util_TableCreateFromString (param_table_handle, options_string)
  if (param_table_handle .lt. 0) then
     call CCTK_Warn(0,52,"surface_integrals-komar_massv2.F90","diagnostics_vacuum","Cannot create parameter table for interpolator")
  endif
  coord_system_handle = -1
  call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
  if (coord_system_handle .lt. 0) then
     call CCTK_Warn(0,58,"surface_integrals-komar_massv2.F90","diagnostics_vacuum","Cannot get handle for cart3d coordinate system &
  &! Forgot to activate an implementation providing coordinates ??")
  endif
  input_array_type_codes = 107
  output_array_type_codes = 107
  ! Specify interpolation input arrays, output arrays:
  N_input_arrays = 26
  N_output_arrays = 26
  call CCTK_VarIndex (input_array_varindices(1), "bssn::phi")
  call CCTK_VarIndex (input_array_varindices(2), "bssn::gupxx")
  call CCTK_VarIndex (input_array_varindices(3), "bssn::gupxy")
  call CCTK_VarIndex (input_array_varindices(4), "bssn::gupxz")
  call CCTK_VarIndex (input_array_varindices(5), "bssn::gupyy")
  call CCTK_VarIndex (input_array_varindices(6), "bssn::gupyz")
  call CCTK_VarIndex (input_array_varindices(7), "bssn::gupzz")
  call CCTK_VarIndex (input_array_varindices(8), "bssn::gxx")
  call CCTK_VarIndex (input_array_varindices(9), "bssn::gxy")
  call CCTK_VarIndex (input_array_varindices(10), "bssn::gxz")
  call CCTK_VarIndex (input_array_varindices(11), "bssn::gyy")
  call CCTK_VarIndex (input_array_varindices(12), "bssn::gyz")
  call CCTK_VarIndex (input_array_varindices(13), "bssn::gzz")
  call CCTK_VarIndex (input_array_varindices(14), "bssn::Axx")
  call CCTK_VarIndex (input_array_varindices(15), "bssn::Axy")
  call CCTK_VarIndex (input_array_varindices(16), "bssn::Axz")
  call CCTK_VarIndex (input_array_varindices(17), "bssn::Ayy")
  call CCTK_VarIndex (input_array_varindices(18), "bssn::Ayz")
  call CCTK_VarIndex (input_array_varindices(19), "bssn::Azz")
  call CCTK_VarIndex (input_array_varindices(20), "bssn::trK")
  call CCTK_VarIndex (input_array_varindices(21), "lapse::lapsex")
  call CCTK_VarIndex (input_array_varindices(22), "lapse::lapsey")
  call CCTK_VarIndex (input_array_varindices(23), "lapse::lapsez")
  call CCTK_VarIndex (input_array_varindices(24), "shift::shiftx")
  call CCTK_VarIndex (input_array_varindices(25), "shift::shifty")
  call CCTK_VarIndex (input_array_varindices(26), "shift::shiftz")
! FIXME?!
!  call CCTK_SyncGroup(i,cctkGH,lapse::lapse_derivatives)
  output_array_pointers(1) = CCTK_PointerTo(Psi4int)
  output_array_pointers(2) = CCTK_PointerTo(gupxxint)
  output_array_pointers(3) = CCTK_PointerTo(gupxyint)
  output_array_pointers(4) = CCTK_PointerTo(gupxzint)
  output_array_pointers(5) = CCTK_PointerTo(gupyyint)
  output_array_pointers(6) = CCTK_PointerTo(gupyzint)
  output_array_pointers(7) = CCTK_PointerTo(gupzzint)
  output_array_pointers(8) = CCTK_PointerTo(gxxint)
  output_array_pointers(9) = CCTK_PointerTo(gxyint)
  output_array_pointers(10) = CCTK_PointerTo(gxzint)
  output_array_pointers(11) = CCTK_PointerTo(gyyint)
  output_array_pointers(12) = CCTK_PointerTo(gyzint)
  output_array_pointers(13) = CCTK_PointerTo(gzzint)
  output_array_pointers(14) = CCTK_PointerTo(Axxint)
  output_array_pointers(15) = CCTK_PointerTo(Axyint)
  output_array_pointers(16) = CCTK_PointerTo(Axzint)
  output_array_pointers(17) = CCTK_PointerTo(Ayyint)
  output_array_pointers(18) = CCTK_PointerTo(Ayzint)
  output_array_pointers(19) = CCTK_PointerTo(Azzint)
  output_array_pointers(20) = CCTK_PointerTo(trKint)
  output_array_pointers(21) = CCTK_PointerTo(lapsexint)
  output_array_pointers(22) = CCTK_PointerTo(lapseyint)
  output_array_pointers(23) = CCTK_PointerTo(lapsezint)
  output_array_pointers(24) = CCTK_PointerTo(shiftxint)
  output_array_pointers(25) = CCTK_PointerTo(shiftyint)
  output_array_pointers(26) = CCTK_PointerTo(shiftzint)
  ! Set up interpolation coordinate arrays:
  n = 1
  do i=1,N_theta
     costheta = 1.D0 - (i - 0.5D0)*dcostheta
     sintheta = sqrt(1.D0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5D0)*dphi
        if(N_phi==1) phiangle = 0.D0
        xinterp(n) = surf_radius*sintheta*cos(phiangle)
        yinterp(n) = surf_radius*sintheta*sin(phiangle)
        zinterp(n) = surf_radius*costheta
        n = n + 1
     end do
  end do
  interp_coords(1) = CCTK_PointerTo(xinterp)
  interp_coords(2) = CCTK_PointerTo(yinterp)
  interp_coords(3) = CCTK_PointerTo(zinterp)
  ! Perform interpolation:
  call CCTK_InterpGridArrays(ierr,cctkGH,N_dims,interp_handle, &
       param_table_handle,coord_system_handle, &
       N_interp_points,input_array_type_codes,interp_coords, &
       N_input_arrays,input_array_varindices, &
       N_output_arrays, output_array_type_codes, output_array_pointers)
  Psi4int = exp(4.D0*Psi4int)
  !Convert Aijint to Kijint!
  Axxint = (Axxint + (1.D0/3.D0)*trKint*gxxint)*Psi4int
  Axyint = (Axyint + (1.D0/3.D0)*trKint*gxyint)*Psi4int
  Axzint = (Axzint + (1.D0/3.D0)*trKint*gxzint)*Psi4int
  Ayyint = (Ayyint + (1.D0/3.D0)*trKint*gyyint)*Psi4int
  Ayzint = (Ayzint + (1.D0/3.D0)*trKint*gyzint)*Psi4int
  Azzint = (Azzint + (1.D0/3.D0)*trKint*gzzint)*Psi4int
  if(ierr.ge.0) then
     !|~~~~~> Integrate the dot product of Dphi with the surface normal:
     output_integral = 0.D0
     do i=1,ntot
        output_integral = output_integral + sqrt(Psi4int(i))*( &
             (lapsexint(i) - Axxint(i)*shiftxint(i) - Axyint(i)*shiftyint(i) - Axzint(i)*shiftzint(i)) &
             *( xinterp(i)*gupxxint(i) + yinterp(i)*gupxyint(i) + zinterp(i)*gupxzint(i) ) + &
             (lapseyint(i) - Axyint(i)*shiftxint(i) - Ayyint(i)*shiftyint(i) - Ayzint(i)*shiftzint(i)) &
             *( xinterp(i)*gupxyint(i) + yinterp(i)*gupyyint(i) + zinterp(i)*gupyzint(i) ) + &
             (lapsezint(i) - Axzint(i)*shiftxint(i) - Ayzint(i)*shiftyint(i) - Azzint(i)*shiftzint(i)) &
             *( xinterp(i)*gupxzint(i) + yinterp(i)*gupyzint(i) + zinterp(i)*gupzzint(i) ))/surf_radius
     end do
     !  Multiply by -1/2\pi * surface element.
     output_integral = output_integral * surf_radius*surf_radius * dphi * dcostheta * (1.D0/(4.D0*PI))*sym_factor
  else
     output_integral = -9999.D0
  end if
end subroutine Komar_surf_integralv2
