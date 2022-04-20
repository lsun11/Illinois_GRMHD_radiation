subroutine update_movingbox(cctk_dim,cctk_gsh,cctk_lsh,cctk_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbox,cctk_delta_time,cc&
  &tk_time,cctk_delta_space,cctk_origin_space,cctk_levfac,cctk_levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,cctk_convfac,cct&
  &k_nghostzones,cctk_iteration,cctkGH, X0scalarMoLanalytic,X0scalarMoLevolve,X0scalarMoLrhs,X0scalarMoLstagger,X0scalarMoLstaggerr&
  &hs,X0scawave_have_global_boundary,X1scalarMoLanalytic,X1scalarMoLevolve,X1scalarMoLrhs,X1scalarMoLstagger,X1scalarMoLstaggerrhs,&
  &X1scawave_have_global_boundary,X2scalarMoLanalytic,X2scalarMoLevolve,X2scalarMoLrhs,X2scalarMoLstagger,X2scalarMoLstaggerrhs,phi&
  &,phi_p,phi_p_p,phi_analytic,phi_analytic_minus_numeric,phi_rhs,phi_stagger,phi_stagger_p,phi_stagger_p_p,phi_stagger_analytic_mi&
  &nus_numeric,phi_stagger_rhs,phidot,phidot_p,phidot_p_p,phidot_rhs,phidot_stagger,phidot_stagger_p,phidot_stagger_p_p,phidot_stag&
  &ger_rhs,scawave_have_bdry_max,scawave_have_bdry_min, X0coordinates,X0radii,X1coordinates,X2coordinates,active_length,num_levels_&
  &length,positions_length,radii_length,active,coarse_dx,coarse_dy,coarse_dz,num_levels,position_x,position_y,position_z,r,radius,x&
  &,y,z)
  implicit none
  INTEGER cctk_dim
  INTEGER cctk_gsh(cctk_dim),cctk_lsh(cctk_dim)
  INTEGER cctk_lbnd(cctk_dim),cctk_ubnd(cctk_dim)
  INTEGER cctk_lssh(3*cctk_dim)
  INTEGER cctk_from(cctk_dim),cctk_to(cctk_dim)
  INTEGER cctk_bbox(2*cctk_dim)
  REAL*8 cctk_delta_time, cctk_time
  REAL*8 cctk_delta_space(cctk_dim)
  REAL*8 cctk_origin_space(cctk_dim)
  INTEGER cctk_levfac(cctk_dim)
  INTEGER cctk_levoff(cctk_dim)
  INTEGER cctk_levoffdenom(cctk_dim)
  INTEGER cctk_timefac
  INTEGER cctk_convlevel
  INTEGER cctk_convfac
  INTEGER cctk_nghostzones(cctk_dim)
  INTEGER cctk_iteration
  integer*8 cctkGH
  INTEGER X0scalarMoLanalytic
  INTEGER X0scalarMoLevolve
  INTEGER X0scalarMoLrhs
  INTEGER X0scalarMoLstagger
  INTEGER X0scalarMoLstaggerrhs
  INTEGER X0scawave_have_global_boundary
  INTEGER X1scalarMoLanalytic
  INTEGER X1scalarMoLevolve
  INTEGER X1scalarMoLrhs
  INTEGER X1scalarMoLstagger
  INTEGER X1scalarMoLstaggerrhs
  INTEGER X1scawave_have_global_boundary
  INTEGER X2scalarMoLanalytic
  INTEGER X2scalarMoLevolve
  INTEGER X2scalarMoLrhs
  INTEGER X2scalarMoLstagger
  INTEGER X2scalarMoLstaggerrhs
  REAL*8 phi(X0scalarMoLevolve,X1scalarMoLevolve,X2scalarMoLevolve)
  REAL*8 phi_p(X0scalarMoLevolve,X1scalarMoLevolve,X2scalarMoLevolve)
  REAL*8 phi_p_p(X0scalarMoLevolve,X1scalarMoLevolve,X2scalarMoLevolve)
  REAL*8 phi_analytic(X0scalarMoLanalytic,X1scalarMoLanalytic,X2scalarMoLanalytic)
  REAL*8 phi_analytic_minus_numeric(X0scalarMoLanalytic,X1scalarMoLanalytic,X2scalarMoLanalytic)
  REAL*8 phi_rhs(X0scalarMoLrhs,X1scalarMoLrhs,X2scalarMoLrhs)
  REAL*8 phi_stagger(X0scalarMoLstagger,X1scalarMoLstagger,X2scalarMoLstagger)
  REAL*8 phi_stagger_p(X0scalarMoLstagger,X1scalarMoLstagger,X2scalarMoLstagger)
  REAL*8 phi_stagger_p_p(X0scalarMoLstagger,X1scalarMoLstagger,X2scalarMoLstagger)
  REAL*8 phi_stagger_analytic_minus_numeric(X0scalarMoLanalytic,X1scalarMoLanalytic,X2scalarMoLanalytic)
  REAL*8 phi_stagger_rhs(X0scalarMoLstaggerrhs,X1scalarMoLstaggerrhs,X2scalarMoLstaggerrhs)
  REAL*8 phidot(X0scalarMoLevolve,X1scalarMoLevolve,X2scalarMoLevolve)
  REAL*8 phidot_p(X0scalarMoLevolve,X1scalarMoLevolve,X2scalarMoLevolve)
  REAL*8 phidot_p_p(X0scalarMoLevolve,X1scalarMoLevolve,X2scalarMoLevolve)
  REAL*8 phidot_rhs(X0scalarMoLrhs,X1scalarMoLrhs,X2scalarMoLrhs)
  REAL*8 phidot_stagger(X0scalarMoLstagger,X1scalarMoLstagger,X2scalarMoLstagger)
  REAL*8 phidot_stagger_p(X0scalarMoLstagger,X1scalarMoLstagger,X2scalarMoLstagger)
  REAL*8 phidot_stagger_p_p(X0scalarMoLstagger,X1scalarMoLstagger,X2scalarMoLstagger)
  REAL*8 phidot_stagger_rhs(X0scalarMoLstaggerrhs,X1scalarMoLstaggerrhs,X2scalarMoLstaggerrhs)
  INTEGER*4 scawave_have_bdry_max(X0scawave_have_global_boundary,X1scawave_have_global_boundary)
  INTEGER*4 scawave_have_bdry_min(X0scawave_have_global_boundary,X1scawave_have_global_boundary)
  INTEGER X0coordinates
  INTEGER X0radii
  INTEGER X1coordinates
  INTEGER X2coordinates
  INTEGER active_length
  INTEGER num_levels_length
  INTEGER positions_length
  INTEGER radii_length
  INTEGER*4 active(active_length)
  REAL*8 coarse_dx
  REAL*8 coarse_dy
  REAL*8 coarse_dz
  INTEGER*4 num_levels(num_levels_length)
  REAL*8 position_x(positions_length)
  REAL*8 position_y(positions_length)
  REAL*8 position_z(positions_length)
  REAL*8 r(X0coordinates,X1coordinates,X2coordinates)
  REAL*8 radius(X0radii,radii_length)
  REAL*8 x(X0coordinates,X1coordinates,X2coordinates)
  REAL*8 y(X0coordinates,X1coordinates,X2coordinates)
  REAL*8 z(X0coordinates,X1coordinates,X2coordinates)
  
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
  INTEGER*4 WaveMoL_MaxNumConstrainedVars
  INTEGER*4 WaveMoL_MaxNumEvolvedVars
  COMMON /scalarwaveMoLrest/WaveMoL_MaxNumConstrainedVars,WaveMoL_MaxNumEvolvedVars
  REAL*8  amplitude
  REAL*8  width
  integer*8  bound
  INTEGER*4 enable_moving_grid
  INTEGER*4 scalarwave_Symmetry
  COMMON /scalarwaveMoLpriv/amplitude,width,bound,enable_moving_grid,scalarwave_Symmetry
  REAL*8  CCTKH0
  REAL*8  CCTKH1
  REAL*8  CCTKH2
  REAL*8  CCTKH3
  REAL*8  CCTKH11
  REAL*8  CCTKH12
  REAL*8  CCTKH13
  REAL*8  CCTKH14
  REAL*8  CCTKH15
  REAL*8  CCTKH16
  REAL*8  CCTKH17
  REAL*8  CCTKH18
  integer*8  bitant_plane
  integer*8  domain
  integer*8  quadrant_direction
  integer*8  CCTKH4
  integer*8  type
  INTEGER*4 CCTKH5
  INTEGER*4 CCTKH6
  INTEGER*4 CCTKH7
  INTEGER*4 CCTKH8
  INTEGER*4 CCTKH9
  INTEGER*4 CCTKH10
  COMMON /GRIDrest/CCTKH0,CCTKH1,CCTKH2,CCTKH3,CCTKH11,CCTKH12,CCTKH13,CCTKH14,CCTKH15,CCTKH16,CCTKH17,CCTKH18,bitant_plane,domain,&
  &quadrant_direction,CCTKH4,type,CCTKH5,CCTKH6,CCTKH7,CCTKH8,CCTKH9,CCTKH10
  INTEGER*4 CCTKH19
  INTEGER*4 CCTKH20
  INTEGER*4 CCTKH21
  INTEGER*4 CCTKH22
  INTEGER*4 CCTKH23
  INTEGER*4 CCTKH24
  INTEGER*4 CCTKH25
  INTEGER*4 CCTKH26
  INTEGER*4 CCTKH27
  INTEGER*4 CCTKH28
  INTEGER*4 CCTKH29
  INTEGER*4 MoL_Num_Constrained_Vars
  INTEGER*4 MoL_Num_Evolved_Vars
  INTEGER*4 MoL_Num_SaveAndRestore_Vars
  INTEGER*4 CCTKH30
  COMMON /METHODOFLINESrest/CCTKH19,CCTKH20,CCTKH21,CCTKH22,CCTKH23,CCTKH24,CCTKH25,CCTKH26,CCTKH27,CCTKH28,CCTKH29,MoL_Num_Constra&
  &ined_Vars,MoL_Num_Evolved_Vars,MoL_Num_SaveAndRestore_Vars,CCTKH30
  
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
  INTEGER*4 function Boundary_SelectGroupForBC (GH, faces, boundary_width, table_handle, var_name, bc_name)
  implicit none
  integer*8 GH
  INTEGER*4 faces
  INTEGER*4 boundary_width
  INTEGER*4 table_handle
  character(*) var_name
  character(*) bc_name
  end function Boundary_SelectGroupForBC
  end interface
  interface
  INTEGER*4 function MoLRegisterConstrained (ConstrainedIndex)
  implicit none
  INTEGER*4 ConstrainedIndex
  end function MoLRegisterConstrained
  end interface
  interface
  INTEGER*4 function MoLRegisterEvolvedGroup (EvolvedIndex, RHSIndex)
  implicit none
  INTEGER*4 EvolvedIndex
  INTEGER*4 RHSIndex
  end function MoLRegisterEvolvedGroup
  end interface
  
  if (enable_moving_grid==1) then
     position_x(1) = 0.96d0*sin(cctk_time*0.5d0)
     position_y(1) = 1.3d0*(1.d0-cos(cctk_time*0.5d0))
     position_z(1) = 0.d0
     write(*,*) 'Moving box center: ',position_x(1),position_y(1),position_z(1)
  end if
end subroutine update_movingbox
