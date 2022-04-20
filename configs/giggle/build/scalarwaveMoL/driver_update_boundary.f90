!-----------------------------------------------------------------------------
! Update outer boundaries
!-----------------------------------------------------------------------------
subroutine scalarwaveMoL_update_boundary(cctk_dim,cctk_gsh,cctk_lsh,cctk_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbox,cctk_&
  &delta_time,cctk_time,cctk_delta_space,cctk_origin_space,cctk_levfac,cctk_levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,cct&
  &k_convfac,cctk_nghostzones,cctk_iteration,cctkGH, X0scalarMoLanalytic,X0scalarMoLevolve,X0scalarMoLrhs,X0scalarMoLstagger,X0scal&
  &arMoLstaggerrhs,X0scawave_have_global_boundary,X1scalarMoLanalytic,X1scalarMoLevolve,X1scalarMoLrhs,X1scalarMoLstagger,X1scalarM&
  &oLstaggerrhs,X1scawave_have_global_boundary,X2scalarMoLanalytic,X2scalarMoLevolve,X2scalarMoLrhs,X2scalarMoLstagger,X2scalarMoLs&
  &taggerrhs,phi,phi_p,phi_p_p,phi_analytic,phi_analytic_minus_numeric,phi_rhs,phi_stagger,phi_stagger_p,phi_stagger_p_p,phi_stagge&
  &r_analytic_minus_numeric,phi_stagger_rhs,phidot,phidot_p,phidot_p_p,phidot_rhs,phidot_stagger,phidot_stagger_p,phidot_stagger_p_&
  &p,phidot_stagger_rhs,scawave_have_bdry_max,scawave_have_bdry_min, X0coordinates,X0radii,X1coordinates,X2coordinates,active_lengt&
  &h,num_levels_length,positions_length,radii_length,active,coarse_dx,coarse_dy,coarse_dz,num_levels,position_x,position_y,position&
  &_z,r,radius,x,y,z)
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
  ;
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
  ;
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
  ;
  integer, dimension(3)                    :: ext,fake_ext,local_have_bdry_min,local_have_bdry_max, gridfunc_syms
  real*8                                   :: dT,dx,dy,dz,Xglobmin,Yglobmin,Zglobmin,Xglobmax,Yglobmax,Zglobmax,levelnumber
  integer                                  :: i,AXISYM,ierr,handle,vindex
  integer                                    :: stagger_x,stagger_y,stagger_z
  AXISYM = 4
  ext = cctk_lsh
  dT = (cctk_delta_time/cctk_timefac)
  dx = (cctk_delta_space(1)/cctk_levfac(1))
  dy = (cctk_delta_space(2)/cctk_levfac(2))
  dz = (cctk_delta_space(3)/cctk_levfac(3))
  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0
  if(CCTK_nProcs(CCTKGH)+1.gt.16384) then
     write(*,*) "Error: Total number of processors exceeds hard-coded limit of 16384 inside update_boundary!"
     write(*,*) "       Please edit scalarwaveMoL/interface.ccl so that scawave_have_global_boundary arrays have > 16384 elements i&
  &n dimension 1."
     stop
  end if
  !First we fill in symmetry ghostzones
  call fill_scalarwave_symmetry_gz(ext,X,Y,Z,scalarwave_Symmetry,phi,phidot)
  ! *** TEST ***
  stagger_x = 1
  stagger_y = 1
  stagger_z = 0
  gridfunc_syms(1) = 1
  gridfunc_syms(2) = 1
  gridfunc_syms(3) = 1
  !!call fill_scalarwave_symmetry_gz(ext,X,Y,Z,scalarwave_Symmetry,phi_stagger,phidot_stagger)
  call scalarwave_do_symmetry_gz_one_gridfunction_stagger(ext,X,Y,Z,Symmetry,phi_stagger, &
                                                gridfunc_syms,stagger_x,stagger_y,stagger_z)
  call scalarwave_do_symmetry_gz_one_gridfunction_stagger(ext,X,Y,Z,Symmetry,phidot_stagger, &
                                                gridfunc_syms,stagger_x,stagger_y,stagger_z)
  ! ************
  if(cctk_iteration.eq.0 .and. levelnumber==1) then
     !Next we find the global coordinate maxima/minima on the grid
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(vindex,"grid::X")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, 107,Xglobmin,1,vindex)
     call CCTK_VarIndex(vindex,"grid::Y")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, 107,Yglobmin,1,vindex)
     call CCTK_VarIndex(vindex,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, 107,Zglobmin,1,vindex)
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"grid::X")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, 107,Xglobmax,1,vindex)
     call CCTK_VarIndex(vindex,"grid::Y")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, 107,Yglobmax,1,vindex)
     call CCTK_VarIndex(vindex,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, 107,Zglobmax,1,vindex)
     !Then we determine whether any of the faces of our local grid corresponds to a global max/min:
     ! This way, we only update actual outer boundaries (and not refinement or symmetry boundaries)x
     call scalarwaveMol_global_edges(dX,dY,dZ, &
          Xglobmin,Yglobmin,Zglobmin, &
          X(1,1,1),Y(1,1,1),Z(1,1,1), &
          Xglobmax,Yglobmax,Zglobmax, &
          X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
          local_have_bdry_min,local_have_bdry_max,scalarwave_Symmetry)
     scawave_have_bdry_min(CCTK_MyProc(CCTKGH)+1,:) = local_have_bdry_min
     scawave_have_bdry_max(CCTK_MyProc(CCTKGH)+1,:) = local_have_bdry_max
     write(*,*) "Just set GLOBAL EDGES for the first time: levnum, blah=",levelnumber,local_have_bdry_min,local_have_bdry_max
  end if
  if(cctk_iteration.gt.0 .and. levelnumber==1) then
     local_have_bdry_min = scawave_have_bdry_min(CCTK_MyProc(CCTKGH)+1,:)
     local_have_bdry_max = scawave_have_bdry_max(CCTK_MyProc(CCTKGH)+1,:)
     write(*,*) "Update boundary: levnum, GLOBAL EDGES:",levelnumber,local_have_bdry_min,local_have_bdry_max
     !Then we update boundaries
     do i=1,cctk_nghostzones(2)
        fake_ext = cctk_lsh - cctk_nghostzones + i
        if(scalarwave_Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
        ! *** TEST ***
        call scalarwaveMoL_update_boundary_lowlevel(ext, fake_ext, X, Y, Z, &
             dX,dY,dZ,dT, &
             phi_p,phi,phidot_p,phidot, phi_stagger_p,phi_stagger, &
             phidot_stagger_p, phidot_stagger, scalarwave_Symmetry, &
             local_have_bdry_min,local_have_bdry_max)
        ! *************
     end do
  end if
end subroutine scalarwaveMoL_update_boundary
