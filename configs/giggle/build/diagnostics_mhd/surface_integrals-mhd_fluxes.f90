!|----------------------------------------------------+
!| Function to compute fluxes of fluid and EM fields 
!!   across a spherical surface of radius surf_radius.
!!   and centered at the origin
!|----------------------------------------------------+
subroutine mhd_surface_fluxes(cctkGH,F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em)
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
  real*8                                   :: F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em
  ! Variables needed for interpolation
  real*8, dimension(N_theta*N_phi,3)       :: pointcoords
  integer                                  :: vindex
  ! Output arrays, dummy indices, parameters
  real*8, dimension(N_theta*N_phi)         :: Psi6int,vxint,vyint,vzint,Pint,u0int
  real*8, dimension(N_theta*N_phi)         :: sbxint,sbyint,sbzint,hint,rho_starint,sbtint
  real*8, dimension(N_theta*N_phi)         :: Psi4int,gxxint,gxyint,gxzint,gyyint,gyzint,gzzint
  real*8, dimension(N_theta*N_phi)         :: lapm1int,lapseconst,shiftxint,shiftyint,shiftzint
  real*8                                   :: phiangle,costheta,sintheta,PI
  real*8                                   :: bx_i,by_i,bz_i,bt_i,vx_i,vy_i,vz_i,ut_i,lap_i
  real*8                                   :: beta_x_i,beta_y_i,beta_z_i,betax_i,betay_i,betaz_i
  real*8                                   :: betabeta,betav,Bbeta,u_t_i,b_t_i,bb,u_phi_i,b_phi_i
  real*8                                   :: gxx_i,gxy_i,gxz_i,gyy_i,gyz_i,gzz_i
  real*8                                   :: rhos_i,h_i,P_i
  real*8                                   :: F_E_emx,F_E_emy,F_E_emz
  real*8                                   :: F_J_emx,F_J_emy,F_J_emz
  real*8                                   :: F_J_fx,F_J_fy,F_J_fz
  real*8                                   :: u_x_i,u_y_i,u_z_i,b_x_i,b_y_i,b_z_i
  real*8                                   :: nx_i,ny_i,nz_i,nx_g,ny_g,nz_g,F_M0_i
  integer                                  :: i,j,n
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  PI = 3.14159265358979323846D0
  !| Psi6int <~~~ exp(6.D0*phi)
  !| Psi4int   <~~~ exp(4.D0*phi),
  !| lapseconst   <~~~ 1.D0/(sqrt(4.0*PI)*(lapm1+1.D0))
  !| lapm1int   <~~~ lapm1,
  !| shiftxint  <~~~ shiftx, shiftyint <~~~ shifty, shiftzint <~~~ shiftz 
  !| gxxint  <~~~ gxx, gxyint <~~~ gxy, gxzint <~~~ gxz, 
  !| gyyint  <~~~ gyy, gyzint <~~~ gyz, gzzint <~~~ gzz, 
  !| sbxint   <~~~ bx,  sbyint  <~~~ by,  sbzint  <~~~ bz,  sbtint  <~~~ bt,
  !| vxint  <~~~ vx, vyint <~~~ vy, vzint <~~~ vz, u0int <~~~ u0
  !| rho_starint  <~~~ rho_star, hint <~~~ h, Pint <~~~ P
  n = 1
  do i=1,N_theta
     costheta = 1.0 - (i - 0.5)*dcostheta
     sintheta = sqrt(1.0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5)*dphi
        if(N_phi==1) phiangle = 0.D0
        pointcoords(n,1) = surf_radius*sintheta*cos(phiangle)
        pointcoords(n,2) = surf_radius*sintheta*sin(phiangle)
        pointcoords(n,3) = surf_radius*costheta
        n = n + 1
     end do
  end do
  call CCTK_VarIndex(vindex,"bssn::phi")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,Psi6int)
  call CCTK_VarIndex(vindex,"lapse::lapm1")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,lapm1int)
  do i=1,N_theta*N_phi
     Psi6int(i) = exp(6.D0*Psi6int(i))
     Psi4int(i) = Psi6int(i)**(4.D0/6.D0)
     lapm1int(i) = lapm1int(i) + 1.D0
     lapseconst(i) = 1.D0/(sqrt(4.D0*PI)*(lapm1int(i)))
  end do
  call CCTK_VarIndex(vindex,"shift::shiftx")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,shiftxint)
  call CCTK_VarIndex(vindex,"shift::shifty")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,shiftyint)
  call CCTK_VarIndex(vindex,"shift::shiftz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,shiftzint)
  call CCTK_VarIndex(vindex,"bssn::gxx")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gxxint)
  call CCTK_VarIndex(vindex,"bssn::gxy")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gxyint)
  call CCTK_VarIndex(vindex,"bssn::gxz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gxzint)
  call CCTK_VarIndex(vindex,"bssn::gyy")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gyyint)
  call CCTK_VarIndex(vindex,"bssn::gyz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gyzint)
  call CCTK_VarIndex(vindex,"bssn::gzz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,gzzint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbx")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,sbxint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sby")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,sbyint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,sbzint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbt")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,sbtint)
  do i=1,N_theta*N_phi
     gxxint(i) = gxxint(i) * Psi4int(i)
     gxyint(i) = gxyint(i) * Psi4int(i)
     gxzint(i) = gxzint(i) * Psi4int(i)
     gyyint(i) = gyyint(i) * Psi4int(i)
     gyzint(i) = gyzint(i) * Psi4int(i)
     gzzint(i) = gzzint(i) * Psi4int(i)
     sbxint(i) = sbxint(i) * lapseconst(i)
     sbyint(i) = sbyint(i) * lapseconst(i)
     sbzint(i) = sbzint(i) * lapseconst(i)
     sbtint(i) = sbtint(i) * lapseconst(i)
  end do
  call CCTK_VarIndex(vindex,"mhd_evolve::vx")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,vxint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vy")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,vyint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,vzint)
  call CCTK_VarIndex(vindex,"mhd_evolve::u0")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,u0int)
  call CCTK_VarIndex(vindex,"mhd_evolve::rho_star")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,rho_starint)
  call CCTK_VarIndex(vindex,"mhd_evolve::h")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,hint)
  call CCTK_VarIndex(vindex,"mhd_evolve::P")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,Pint)
  !Zero out our flux sums first.
  F_M0      = 0.D0
  F_E_fluid = 0.D0
  F_E_em    = 0.D0
  F_J_fluid = 0.D0
  F_J_em    = 0.D0
  !|~~~~~> Integrate the dot product of Dphi with the surface normal
  do i=1,N_theta*N_phi
     bx_i = sbxint(i)
     by_i = sbyint(i)
     bz_i = sbzint(i)
     bt_i = sbtint(i)
     vx_i = vxint(i)
     vy_i = vyint(i)
     vz_i = vzint(i)
     ut_i = u0int(i)
     lap_i = lapm1int(i)
     betax_i = shiftxint(i)
     betay_i = shiftyint(i)
     betaz_i = shiftzint(i)
     gxx_i = gxxint(i)
     gxy_i = gxyint(i)
     gxz_i = gxzint(i)
     gyy_i = gyyint(i)
     gyz_i = gyzint(i)
     gzz_i = gzzint(i)
     rhos_i = rho_starint(i)
     h_i = hint(i)
     P_i = Pint(i)
     nx_g = pointcoords(i,1)/surf_radius
     ny_g = pointcoords(i,2)/surf_radius
     nz_g = pointcoords(i,3)/surf_radius
     nx_i = pointcoords(i,1)
     ny_i = pointcoords(i,2)
     nz_i = pointcoords(i,3)
     ! \beta_i
     beta_x_i = gxxint(i)*betax_i + gxyint(i)*betay_i + gxzint(i)*betaz_i
     beta_y_i = gxyint(i)*betax_i + gyyint(i)*betay_i + gyzint(i)*betaz_i
     beta_z_i = gxzint(i)*betax_i + gyzint(i)*betay_i + gzzint(i)*betaz_i
     ! \beta.\beta
     betabeta = betax_i*beta_x_i + betay_i*beta_y_i + betaz_i*beta_z_i
     ! \beta.v
     betav = vx_i*beta_x_i + vy_i*beta_y_i + vz_i*beta_z_i
     ! \beta.b
     Bbeta = bx_i*beta_x_i + by_i*beta_y_i + bz_i*beta_z_i
     ! u_t
     u_t_i = u0int(i)*(-lap_i*lap_i + betabeta + betav)
     ! b_t
     b_t_i = bt_i*(-lap_i*lap_i + betabeta) + Bbeta
     ! u_x, u_y, and u_phi
     u_x_i = ut_i*( gxxint(i)*(vx_i + betax_i) + gxy_i*(vy_i + betay_i) &
          + gxz_i*(vz_i + betaz_i) )
     u_y_i = ut_i*( gxy_i*(vx_i + betax_i) + gyy_i*(vy_i + betay_i) &
          + gyz_i*(vz_i + betaz_i) )
     u_phi_i = nx_i*u_y_i - ny_i*u_x_i
     ! b_x, b_y, and b_phi
     b_x_i = beta_x_i*bt_i + gxx_i*bx_i + gxy_i*by_i + gxz_i*bz_i
     b_y_i = beta_y_i*bt_i + gxy_i*bx_i + gyy_i*by_i + gyz_i*bz_i
     b_phi_i = nx_i*b_y_i - ny_i*b_x_i
     ! b.b
     bb = -lap_i*lap_i*bt_i*bt_i &
          + gxx_i*(bx_i+betax_i*bt_i)*(bx_i+betax_i*bt_i) &
          + gyyint(i)*(by_i+betay_i*bt_i)*(by_i+betay_i*bt_i) &
          + gzzint(i)*(bz_i+betaz_i*bt_i)*(bz_i+betaz_i*bt_i) &
          + 2.D0*( gxyint(i)*(bx_i+betax_i*bt_i)*(by_i+betay_i*bt_i) &
                  + gxzint(i)*(bx_i+betax_i*bt_i)*(bz_i+betaz_i*bt_i) &
                  + gyzint(i)*(by_i+betay_i*bt_i)*(bz_i+betaz_i*bt_i) )
     F_E_emx = bb*ut_i*vx_i*u_t_i - bx_i*b_t_i
     F_E_emy = bb*ut_i*vy_i*u_t_i - by_i*b_t_i
     F_E_emz = bb*ut_i*vz_i*u_t_i - bz_i*b_t_i
     F_J_fx = rhos_i*h_i*vx_i*u_phi_i - lap_i*Psi6int(i)*P_i*ny_i
     F_J_fy = rhos_i*h_i*vy_i*u_phi_i + lap_i*Psi6int(i)*P_i*nx_i
     F_J_fz = rhos_i*h_i*vz_i*u_phi_i
     F_J_emx = (ut_i*vx_i*u_phi_i - 0.5*ny_i)*bb - bx_i*b_phi_i
     F_J_emy = (ut_i*vy_i*u_phi_i + 0.5*nx_i)*bb - by_i*b_phi_i
     F_J_emz = ut_i*vz_i*u_phi_i*bb - bz_i*b_phi_i
     F_M0_i    = (nx_g*vx_i + ny_g*vy_i + nz_g*vz_i)*rhos_i
     F_M0      = F_M0 + F_M0_i
     F_E_fluid = F_E_fluid - F_M0_i*h_i*u_t_i
     F_E_em    = F_E_em - (nx_g*F_E_emx + ny_g*F_E_emy + nz_g*F_E_emz)*lap_i*Psi6int(i)
     F_J_fluid = F_J_fluid + (nx_g*F_J_fx + ny_g*F_J_fy + nz_g*F_J_fz)
     F_J_em    = F_J_em + (nx_g*F_J_emx + ny_g*F_J_emy + nz_g*F_J_emz)*lap_i*Psi6int(i)
  end do
  !|  Multiply by surface element.
  F_M0      = F_M0 * surf_radius * surf_radius * dphi * dcostheta*sym_factor
  F_E_fluid = F_E_fluid * surf_radius * surf_radius * dphi * dcostheta*sym_factor
  ! F_E_em is the Poynting flux through the surface at surf_radius !
  F_E_em    = F_E_em *  surf_radius * surf_radius * dphi * dcostheta*sym_factor
  F_J_fluid = F_J_fluid *  surf_radius * surf_radius * dphi * dcostheta*sym_factor
  F_J_em    = F_J_em *  surf_radius * surf_radius * dphi * dcostheta*sym_factor
  !  write(*,*) JFINAL,N_phi,N_theta,surf_radius,F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em
end subroutine mhd_surface_fluxes
