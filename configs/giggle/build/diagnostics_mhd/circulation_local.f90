!-----------------------------------------------------------------------------
! Compute c= \int_0^2pi h u^m dS_m.
!  Note that we work around a surface centered at xcenter,ycenter,zcenter
!-----------------------------------------------------------------------------
subroutine circulation(cctkGH,x_cm,y_cm,z_cm,r_local,output_integral)
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
  REAL*8  CoMExcision_Radius
  REAL*8  M_ADM
  REAL*8  initial_monopole_value
  REAL*8  inner_lum_rad_ratio
  REAL*8  lum_outer_rad
  REAL*8  radius_esc1
  REAL*8  radius_esc2
  REAL*8  radius_esc3
  REAL*8  radius_esc4
  REAL*8  rhob_cutoff
  REAL*8  rhosurf
  REAL*8  rhosurf_rmax
  REAL*8  rhosurf_rmin
  REAL*8  u0sch
  INTEGER*4 N_rad_ray
  INTEGER*4 Nphi_points
  INTEGER*4 Ntheta_points
  INTEGER*4 const_rad_surf_diagnostics
  INTEGER*4 drho_dtau_calc_enable
  INTEGER*4 escape_mass_diag
  INTEGER*4 luminosity_diagnostics
  COMMON /diagnostics_mhdrest/CoMExcision_Radius,M_ADM,initial_monopole_value,inner_lum_rad_ratio,lum_outer_rad,radius_esc1,radius_&
  &esc2,radius_esc3,radius_esc4,rhob_cutoff,rhosurf,rhosurf_rmax,rhosurf_rmin,u0sch,N_rad_ray,Nphi_points,Ntheta_points,const_rad_s&
  &urf_diagnostics,drho_dtau_calc_enable,escape_mass_diag,luminosity_diagnostics
  REAL*8  CCTKH5
  REAL*8  CCTKH6
  REAL*8  CCTKH7
  REAL*8  CCTKH8
  REAL*8  binary_orb_freq
  REAL*8  CCTKH11
  REAL*8  CCTKH12
  REAL*8  CCTKH13
  REAL*8  CCTKH14
  REAL*8  CCTKH15
  REAL*8  CCTKH16
  REAL*8  rah1_initial
  REAL*8  rah2_initial
  REAL*8  CCTKH17
  REAL*8  xbh1_initial
  REAL*8  xbh2_initial
  REAL*8  CCTKH19
  INTEGER*4 CCTKH0
  INTEGER*4 CCTKH1
  INTEGER*4 CCTKH2
  INTEGER*4 CCTKH3
  INTEGER*4 CCTKH4
  INTEGER*4 CCTKH9
  INTEGER*4 CCTKH10
  INTEGER*4 CCTKH18
  INTEGER*4 CCTKH20
  COMMON /BSSNrest/CCTKH5,CCTKH6,CCTKH7,CCTKH8,binary_orb_freq,CCTKH11,CCTKH12,CCTKH13,CCTKH14,CCTKH15,CCTKH16,rah1_initial,rah2_in&
  &itial,CCTKH17,xbh1_initial,xbh2_initial,CCTKH19,CCTKH0,CCTKH1,CCTKH2,CCTKH3,CCTKH4,CCTKH9,CCTKH10,CCTKH18,CCTKH20
  REAL*8  CCTKH21
  REAL*8  dcostheta
  REAL*8  ddrbddr
  REAL*8  dphi
  REAL*8  drbdr
  REAL*8  CCTKH26
  REAL*8  rbr
  REAL*8  CCTKH30
  REAL*8  surf_radius
  INTEGER*4 Compute_VolIntegrands_Every
  INTEGER*4 N_phi
  INTEGER*4 N_theta
  INTEGER*4 WhichIntegral
  INTEGER*4 CCTKH22
  INTEGER*4 CCTKH23
  INTEGER*4 CCTKH24
  INTEGER*4 CCTKH25
  INTEGER*4 CCTKH27
  INTEGER*4 ntot
  INTEGER*4 num_BHs
  INTEGER*4 CCTKH28
  INTEGER*4 CCTKH29
  INTEGER*4 CCTKH31
  INTEGER*4 sym_factor
  COMMON /DIAGNOSTICS_VACUUMrest/CCTKH21,dcostheta,ddrbddr,dphi,drbdr,CCTKH26,rbr,CCTKH30,surf_radius,Compute_VolIntegrands_Every,N&
  &_phi,N_theta,WhichIntegral,CCTKH22,CCTKH23,CCTKH24,CCTKH25,CCTKH27,ntot,num_BHs,CCTKH28,CCTKH29,CCTKH31,sym_factor
  REAL*8  CCTKH32
  REAL*8  CCTKH33
  REAL*8  CCTKH34
  REAL*8  CCTKH35
  REAL*8  CCTKH45
  REAL*8  CCTKH46
  REAL*8  CCTKH47
  REAL*8  CCTKH48
  REAL*8  CCTKH49
  REAL*8  CCTKH50
  REAL*8  CCTKH51
  REAL*8  CCTKH52
  integer*8  bitant_plane
  integer*8  domain
  integer*8  CCTKH36
  integer*8  CCTKH37
  integer*8  CCTKH44
  INTEGER*4 CCTKH38
  INTEGER*4 CCTKH39
  INTEGER*4 CCTKH40
  INTEGER*4 CCTKH41
  INTEGER*4 CCTKH42
  INTEGER*4 CCTKH43
  COMMON /GRIDrest/CCTKH32,CCTKH33,CCTKH34,CCTKH35,CCTKH45,CCTKH46,CCTKH47,CCTKH48,CCTKH49,CCTKH50,CCTKH51,CCTKH52,bitant_plane,dom&
  &ain,CCTKH36,CCTKH37,CCTKH44,CCTKH38,CCTKH39,CCTKH40,CCTKH41,CCTKH42,CCTKH43
  REAL*8  CCTKH55
  REAL*8  CCTKH56
  REAL*8  CCTKH57
  integer*8  slicing_type
  INTEGER*4 CCTKH53
  INTEGER*4 CCTKH54
  INTEGER*4 CCTKH58
  INTEGER*4 CCTKH59
  INTEGER*4 CCTKH60
  INTEGER*4 CCTKH61
  INTEGER*4 CCTKH62
  COMMON /LAPSErest/CCTKH55,CCTKH56,CCTKH57,slicing_type,CCTKH53,CCTKH54,CCTKH58,CCTKH59,CCTKH60,CCTKH61,CCTKH62
  REAL*8  CCTKH65
  REAL*8  CCTKH66
  REAL*8  K_poly
  REAL*8  CCTKH69
  REAL*8  CCTKH71
  REAL*8  P_fact
  REAL*8  CCTKH72
  REAL*8  CCTKH73
  REAL*8  Sym_Bz
  REAL*8  CCTKH75
  REAL*8  CCTKH79
  REAL*8  CCTKH80
  REAL*8  CCTKH83
  REAL*8  CCTKH90
  REAL*8  ergo_sigma
  REAL*8  CCTKH92
  REAL*8  gamma_th
  REAL*8  CCTKH94
  REAL*8  CCTKH96
  REAL*8  CCTKH97
  REAL*8  CCTKH98
  REAL*8  CCTKH99
  REAL*8  CCTKH101
  REAL*8  pfloor
  REAL*8  CCTKH105
  REAL*8  CCTKH106
  REAL*8  CCTKH107
  REAL*8  CCTKH108
  REAL*8  CCTKH109
  REAL*8  CCTKH111
  REAL*8  CCTKH112
  REAL*8  CCTKH115
  REAL*8  CCTKH116
  REAL*8  rho_b_atm
  REAL*8  rho_b_max
  REAL*8  rho_fact
  REAL*8  CCTKH119
  REAL*8  CCTKH120
  REAL*8  CCTKH121
  REAL*8  tau_atm
  REAL*8  tau_fact
  INTEGER*4 CCTKH63
  INTEGER*4 CCTKH64
  INTEGER*4 CCTKH67
  INTEGER*4 CCTKH68
  INTEGER*4 CCTKH70
  INTEGER*4 CCTKH74
  INTEGER*4 CCTKH76
  INTEGER*4 CCTKH77
  INTEGER*4 CCTKH78
  INTEGER*4 CCTKH81
  INTEGER*4 constrained_transport_scheme
  INTEGER*4 CCTKH82
  INTEGER*4 em_evolve_enable
  INTEGER*4 CCTKH84
  INTEGER*4 CCTKH85
  INTEGER*4 enable_OS_collapse
  INTEGER*4 CCTKH86
  INTEGER*4 CCTKH87
  INTEGER*4 CCTKH88
  INTEGER*4 CCTKH89
  INTEGER*4 eps_flag
  INTEGER*4 ergo_star
  INTEGER*4 CCTKH91
  INTEGER*4 CCTKH93
  INTEGER*4 hyperbolic_divergence_cleaning_centered_differencing
  INTEGER*4 hyperbolic_divergence_cleaning_enable
  INTEGER*4 CCTKH95
  INTEGER*4 CCTKH100
  INTEGER*4 neos
  INTEGER*4 CCTKH102
  INTEGER*4 CCTKH103
  INTEGER*4 CCTKH104
  INTEGER*4 CCTKH110
  INTEGER*4 rad_evolve_enable
  INTEGER*4 CCTKH113
  INTEGER*4 CCTKH114
  INTEGER*4 CCTKH117
  INTEGER*4 CCTKH118
  INTEGER*4 CCTKH122
  INTEGER*4 CCTKH123
  INTEGER*4 CCTKH124
  INTEGER*4 CCTKH125
  COMMON /MHD_EVOLVErest/CCTKH65,CCTKH66,K_poly,CCTKH69,CCTKH71,P_fact,CCTKH72,CCTKH73,Sym_Bz,CCTKH75,CCTKH79,CCTKH80,CCTKH83,CCTKH&
  &90,ergo_sigma,CCTKH92,gamma_th,CCTKH94,CCTKH96,CCTKH97,CCTKH98,CCTKH99,CCTKH101,pfloor,CCTKH105,CCTKH106,CCTKH107,CCTKH108,CCTKH&
  &109,CCTKH111,CCTKH112,CCTKH115,CCTKH116,rho_b_atm,rho_b_max,rho_fact,CCTKH119,CCTKH120,CCTKH121,tau_atm,tau_fact,CCTKH63,CCTKH64&
  &,CCTKH67,CCTKH68,CCTKH70,CCTKH74,CCTKH76,CCTKH77,CCTKH78,CCTKH81,constrained_transport_scheme,CCTKH82,em_evolve_enable,CCTKH84,C&
  &CTKH85,enable_OS_collapse,CCTKH86,CCTKH87,CCTKH88,CCTKH89,eps_flag,ergo_star,CCTKH91,CCTKH93,hyperbolic_divergence_cleaning_cent&
  &ered_differencing,hyperbolic_divergence_cleaning_enable,CCTKH95,CCTKH100,neos,CCTKH102,CCTKH103,CCTKH104,CCTKH110,rad_evolve_ena&
  &ble,CCTKH113,CCTKH114,CCTKH117,CCTKH118,CCTKH122,CCTKH123,CCTKH124,CCTKH125
  INTEGER*4 Spatial_Gauge
  INTEGER*4 CCTKH126
  INTEGER*4 CCTKH127
  COMMON /SHIFTrest/Spatial_Gauge,CCTKH126,CCTKH127
  
  ! Variables in the function call
  integer*8                             :: cctkGH
  real*8                                   :: x_cm,y_cm,z_cm,r_local
  real*8                                   :: output_integral
  ! Variables needed for interpolation
  integer*8, dimension(3)               :: interp_coords
  character(60)                            :: options_string
  integer                                  :: interpolation_order,nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle,N_interp_points
  integer                                  :: N_input_arrays,N_output_arrays
  real*8, dimension(Nphi_points)                 :: xinterp,yinterp,zinterp
  integer,dimension(15)                    :: input_array_type_codes,input_array_varindices,output_array_type_codes
  integer*8,dimension(15)               :: output_array_pointers
  ! Output arrays, dummy indices, parameters
  real*8,dimension(Nphi_points)            :: Psi4i,gxxi,gxyi,gxzi,gyyi,gyzi,gzzi
  real*8,dimension(Nphi_points)            :: shiftxi,shiftyi,shiftzi,lamp1i,vxi,vyi,vzi,hi
  real*8,dimension(Nphi_points)            :: v2
  real*8                                   :: phiangle,PI,nxp,nyp,nzp
  real*8                                   :: alp_u02,u02,u0,u_x,u_y
  integer                                  :: i,j,n
  PI = 3.14159265358979323846D0
  N_dims = 3
  interpolation_order = 2
  N_interp_points = Nphi_points
  interp_handle = -1
  call CCTK_InterpHandle (interp_handle, "Lagrange polynomial interpolation")
!  call CCTK_InterpHandle (interp_handle, uniform cartesian)
  if (interp_handle .lt. 0) then
     call CCTK_Warn(0,46,"circulation_local.F90","diagnostics_mhd","Cannot get handle for interpolation ! Forgot to activate an imp&
  &lementation providing interpolation operators ??")
  endif
  param_table_handle = -1
  options_string = "order = " // char(ichar('0') + interpolation_order)
  call Util_TableCreateFromString (param_table_handle, options_string)
  if (param_table_handle .lt. 0) then
     call CCTK_Warn(0,53,"circulation_local.F90","diagnostics_mhd","Cannot create parameter table for interpolator")
  endif
  coord_system_handle = -1
  call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
  if (coord_system_handle .lt. 0) then
     call CCTK_Warn(0,59,"circulation_local.F90","diagnostics_mhd","Cannot get handle for cart3d coordinate system ! Forgot to acti&
  &vate an implementation providing coordinates ??")
  endif
  input_array_type_codes = 107
  output_array_type_codes = 107
  ! Specify interpolation input arrays, output arrays:
  N_input_arrays = 15
  N_output_arrays = 15
!  print *,before call arrays
  call CCTK_VarIndex (input_array_varindices(1), "bssn::phi")
  call CCTK_VarIndex (input_array_varindices(2), "bssn::gxx")
  call CCTK_VarIndex (input_array_varindices(3), "bssn::gxy")
  call CCTK_VarIndex (input_array_varindices(4), "bssn::gxz")
  call CCTK_VarIndex (input_array_varindices(5), "bssn::gyy")
  call CCTK_VarIndex (input_array_varindices(6), "bssn::gyz")
  call CCTK_VarIndex (input_array_varindices(7), "bssn::gzz")
  call CCTK_VarIndex (input_array_varindices(8),  "shift::shiftx")
  call CCTK_VarIndex (input_array_varindices(9),  "shift::shifty")
  call CCTK_VarIndex (input_array_varindices(10), "shift::shiftz")
  call CCTK_VarIndex (input_array_varindices(11), "lapse::lapm1")
  call CCTK_VarIndex (input_array_varindices(12), "mhd_evolve::vx")
  call CCTK_VarIndex (input_array_varindices(13), "mhd_evolve::vy")
  call CCTK_VarIndex (input_array_varindices(14), "mhd_evolve::vz")
!  call CCTK_VarIndex(vindex,mhd_evolve::vx)
  call CCTK_VarIndex (input_array_varindices(15), "mhd_evolve::h")
!  print *,after call arrays
  output_array_pointers(1) = CCTK_PointerTo(Psi4i)
  output_array_pointers(2) = CCTK_PointerTo(gxxi)
  output_array_pointers(3) = CCTK_PointerTo(gxyi)
  output_array_pointers(4) = CCTK_PointerTo(gxzi)
  output_array_pointers(5) = CCTK_PointerTo(gyyi)
  output_array_pointers(6) = CCTK_PointerTo(gyzi)
  output_array_pointers(7) = CCTK_PointerTo(gzzi)
  output_array_pointers(8) = CCTK_PointerTo(shiftxi)
  output_array_pointers(9) = CCTK_PointerTo(shiftyi)
  output_array_pointers(10) = CCTK_PointerTo(shiftzi)
  output_array_pointers(11) = CCTK_PointerTo(lamp1i)
  output_array_pointers(12) = CCTK_PointerTo(vxi)
  output_array_pointers(13) = CCTK_PointerTo(vyi)
  output_array_pointers(14) = CCTK_PointerTo(vzi)
  output_array_pointers(15) = CCTK_PointerTo(hi)
!  print *,output  array ponters
  ! Set up interpolation coordinate arrays:
  n = 1
  dphi = 2.0d0*PI/Nphi_points
  do j=1,Nphi_points
!     phiangle = (j - 0.5D0)*dphi
     phiangle = j*dphi
!     if(N_phi==1) phiangle = 0.D0
     xinterp(n) = r_local*cos(phiangle)+x_cm
     yinterp(n) = r_local*sin(phiangle)+y_cm
     zinterp(n) = 0.0
     n = n + 1
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
!  print *,after interpolation
! define new arrays
  psi4i  = exp(4.d0*psi4i)
  lamp1i = 1.0+lamp1i
  gxxi   = psi4i*gxxi
  gxyi   = psi4i*gxyi
  gxzi   = psi4i*gxzi
  gyyi   = psi4i*gyyi
  gyzi   = psi4i*gyzi
  gzzi   = psi4i*gzzi
  v2     = (    gxxi*(vxi + shiftxi)*(vxi + shiftxi)  &
        + 2.0d0*gxyi*(vxi + shiftxi)*(vyi + shiftyi)  &
        + 2.0d0*gxzi*(vxi + shiftxi)*(vzi + shiftzi)  &
        +       gyyi*(vyi + shiftyi)*(vyi + shiftyi)  &
        + 2.0d0*gyzi*(vyi + shiftyi)*(vzi + shiftzi)  &
        +       gzzi*(vzi + shiftzi)*(vzi + shiftzi) ) /  &
        (lamp1i*lamp1i)
  if(ierr.ge.0) then
     !|~~~~~> Integrate the dot product of Dphi with the surface normal:
     output_integral = 0.D0
     do i=1,N_interp_points
        alp_u02 = 1.0/(1.0-v2(i))
        u02 = alp_u02 /(lamp1i(i)*lamp1i(i))
        u0 = sqrt(u02)
        !compute u_i
        u_x = u0*(gxxi(i)*(vxi(i)+shiftxi(i)) + gxyi(i)*(vyi(i)+shiftyi(i)) + &
              gxzi(i)*(vzi(i)+shiftzi(i)) )
        u_y = u0*(gxyi(i)*(vxi(i)+shiftxi(i)) + gyyi(i)*(vyi(i)+shiftyi(i)) + &
              gyzi(i)*(vzi(i)+shiftzi(i)) )
        !tangent vector
        nxp= -(yinterp(i)-y_cm)
        nyp=  (xinterp(i)-x_cm)
        output_integral = output_integral + hi(i)*(nxp*u_x +nyp*u_y)
     end do
     output_integral = output_integral * dphi
  else
     output_integral = -9999.D0
  end if
end subroutine circulation
