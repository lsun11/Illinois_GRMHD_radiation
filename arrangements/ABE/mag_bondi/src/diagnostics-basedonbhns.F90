#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for bondi thorn
!-----------------------------------------------------------------------------
subroutine diagnostic_mag_bondi_accretion(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,Zmin,dt_diag
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20

  integer                                  :: vindex,N_points
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  !Needed for convenient constraint output
  real*8 :: constraint_ratio
  !End of convenient constraint output variables

  !Needed for GW extraction:
  real*8                                   :: Mass_GW,hplus,hcross,moncrief_gw_extr_radius
  real*8                                   :: psi4r,psi4i
  integer                                  :: jj,nmodes_moncrief
  real*8,dimension(2)                      :: gw_amplitudes
  !End GW extraction variables  

  ! Variables used for BH local diagnotsics:
  integer :: found_horizon, horizon_number
  integer :: N_theta_bh_diagnostics, N_phi_bh_diagnostics
  real*8  :: Mirr, M_BH, J_BH, xi_err
  ! End BH local diagnotsics variables

  real*8                                   :: old_surf_radius,SymmFactor,horizdirn_x,horizdirn_y,horizdirn_z,out_radius
  real*8                                   :: rho_tiny,M0dot,F_E_fluid,F_E_em,F_J_fluid,F_J_em,RestMass
  real*8, dimension(1,3)                   :: pointcoords
  CCTK_POINTER :: reduction_value_pointer
  integer :: header_flag,handle
  integer :: num_cols
  integer :: ierr,myproc_rank,ii

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(MOD(cctk_iteration,out_every)==0) then
     call set_SymmFactor(Symmetry,SymmFactor)

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1
     else
        header_flag = 0
     end if

     !--------------------------------------------------------------------------------!
     filename = 'bondi.mon'
     num_cols = 16
     export_data = 0.D0
     data_headers = 'unknowndude'

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Phicent'
     call CCTK_VarIndex (vindex, "bssn::phi")
     call interp_center_carp(cctkGH,vindex,export_data(2))

     data_headers(3) = 'Kcent'
     call CCTK_VarIndex (vindex, "bssn::trK")
     call interp_center_carp(cctkGH,vindex,export_data(3))

     data_headers(4) = 'lapsecent'
     call CCTK_VarIndex (vindex, "lapse::lapm1")
     call interp_center_carp(cctkGH,vindex,export_data(4))
     export_data(4) = export_data(4) + 1.D0

     data_headers(5) = 'psicent'
     !     Psi = exp(phi)
     call CCTK_VarIndex (vindex, "bssn::psi")
     call interp_center_carp(cctkGH,vindex,export_data(5))

     data_headers(6) = 'lapsemin'
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(vindex,"lapse::lapm1")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(6),1,vindex)
     export_data(6) = export_data(6) + 1.D0

     data_headers(7) = 'phimax'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"bssn::phi")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(7),1,vindex)

     data_headers(8) = 'rhocent'
     call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
     call interp_center_carp(cctkGH,vindex,export_data(8))

     data_headers(9) = 'rho_b_max'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(9),1,vindex)

     data_headers(10) = 'Jsurf_BH_1.1'
     old_surf_radius = surf_radius
     surf_radius=1.1D0
     call J_surf_integral_offcenter(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(10))
     surf_radius = old_surf_radius

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(11) = "J_1.1exVol"
     export_data(11) = J_ADM_VolInt*SymmFactor

     data_headers(12) = "M_ADMsurf_BH_1.1"
     old_surf_radius = surf_radius
     surf_radius=1.1D0
     call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(12))
     surf_radius = old_surf_radius

     ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(13) = "M_ADM1.1exVol"
     export_data(13) = M_ADM_VolInt*SymmFactor

     data_headers(14) = 'rhostar_min'
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex (vindex, "mhd_evolve::rho_star")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(14),1,vindex)

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(15) = "J_1.1exVolr"
     export_data(15) = J_ADM_VolInt_inner*SymmFactor

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(16) = "J_1.1exVolr/2"
     export_data(16) = J_ADM_VolInt_inner2*SymmFactor

     !     data_headers(8) = 'Omega_c'
     !     Pr = (vy*X - vx*Y) / (X*X + Y*Y)
     !     call CCTK_VarIndex (vindex, "mhd_evolve::Pr")
     !     call interp_center_carp(cctkGH,vindex,export_data(8))

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     filename = 'bondi.don'

     num_cols = 7
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     ! M0 volume integral is computed in scheduler of diagnostics_mhd
     data_headers(2) = "M0"
     export_data(2) = M0_VolInt*SymmFactor

     ! M0 inside AH volume integral is computed in scheduler of diagnostics_mhd
     data_headers(3) = "M0 inside AH"
     export_data(3) = M0_AH_VolInt*SymmFactor

     ! Eint (minternal_hybrid) volume integral is computed in scheduler of diagnostics_mhd
     data_headers(4) = 'Eint'
     export_data(4) = Minternal_VolInt*SymmFactor

     ! Kinetic energy (T) volume integral is computed in scheduler of diagnostics_mhd
     data_headers(5) = 'T'
     export_data(5) = T_VolInt*SymmFactor

     ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(6) = "M_ADM_Vol"
     export_data(6) = M_ADM_VolInt*SymmFactor

     data_headers(7) = "M_Komar Surf"
     call Komar_surf_integralv2(cctkGH,export_data(7))

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     filename = 'bondi.kon'
     num_cols = 8
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Xbh1'
     export_data(2) = bh_posn_x(1)

     data_headers(3) = 'Ybh1'
     export_data(3) = bh_posn_y(1)

     data_headers(4) = 'Zbh1'
     export_data(4) = bh_posn_z(1)

     data_headers(5) = 'X_CoM'
     export_data(5) = CoMx_VolInt/CoM_VolInt_denominator
     data_headers(6) = 'Y_CoM'
     export_data(6) = CoMy_VolInt/CoM_VolInt_denominator
     data_headers(7) = 'Z_CoM'
     if(Symmetry==EQUATORIAL) then
        export_data(7) = 0.D0
     else
        export_data(7) = CoMz_VolInt/CoM_VolInt_denominator
     end if

     data_headers(8) = 'rho_b_CoM'
     call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
     pointcoords(1,1) = CoMx_VolInt/CoM_VolInt_denominator
     pointcoords(1,2) = CoMy_VolInt/CoM_VolInt_denominator
     pointcoords(1,3) = 0.D0
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(8))

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     filename = 'bondi.pon'
     num_cols = 8
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'phicent'
     call CCTK_VarIndex(vindex,"bssn::phi")

     pointcoords(1,1) = bh_posn_x(1)
     pointcoords(1,2) = bh_posn_y(1)
     pointcoords(1,3) = 0.D0
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(2))

     data_headers(3) = 'phixplus'
     call CCTK_VarIndex(vindex,"bssn::phi")
     horizdirn_x = 100000.D0
     horizdirn_y = 0.D0
     horizdirn_z = 0.D0
     call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
     write(*,*) "inside bondi diagnostics: out_radius = ",out_radius
     pointcoords(1,1) = bh_posn_x(1)+out_radius
     pointcoords(1,2) = bh_posn_y(1)
     pointcoords(1,3) = 0.D0
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(3))

     data_headers(4) = 'phiyplus'
     call CCTK_VarIndex(vindex,"bssn::phi")
     horizdirn_x = 0.D0
     horizdirn_y = 100000.D0
     horizdirn_z = 0.D0
     call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
     pointcoords(1,1) = bh_posn_x(1)
     pointcoords(1,2) = bh_posn_y(1)+out_radius
     pointcoords(1,3) = 0.D0
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(4))

     data_headers(5) = 'phixminus'
     call CCTK_VarIndex(vindex,"bssn::phi")
     horizdirn_x = -100000.D0
     horizdirn_y = 0.D0
     horizdirn_z = 0.D0
     call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
     pointcoords(1,1) = bh_posn_x(1)-out_radius
     pointcoords(1,2) = bh_posn_y(1)
     pointcoords(1,3) = 0.D0
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(5))

     data_headers(6) = 'phiyminus'
     call CCTK_VarIndex(vindex,"bssn::phi")
     horizdirn_x = 0.D0
     horizdirn_y = -100000.D0
     horizdirn_z = 0.D0
     call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
     pointcoords(1,1) = bh_posn_x(1)
     pointcoords(1,2) = bh_posn_y(1)-out_radius
     pointcoords(1,3) = 0.D0
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(6))

     data_headers(7) = 'phizplus'
     call CCTK_VarIndex(vindex,"bssn::phi")
     horizdirn_x = 0.D0
     horizdirn_y = 0.D0
     horizdirn_z = 100000.D0
     call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
     pointcoords(1,1) = bh_posn_x(1)
     pointcoords(1,2) = bh_posn_y(1)
     pointcoords(1,3) = out_radius
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(7))

     data_headers(8) = 'phixminus0.5'
     call CCTK_VarIndex(vindex,"bssn::phi")
     horizdirn_x = -100000.D0
     horizdirn_y = 0.D0
     horizdirn_z = 0.D0
     call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
     pointcoords(1,1) = bh_posn_x(1)-out_radius*0.5D0
     pointcoords(1,2) = bh_posn_y(1)
     pointcoords(1,3) = 0.D0
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(8))

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)     

     !--------------------------------------------------------------------------------!
     !Box1X_VolInt,Box1Y_VolInt,Box1Z_VolInt,Box1denom_VolInt
     filename = 'bondi.xon'
     num_cols = 3
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Box1x'
     export_data(2) = position_x(2) 
     data_headers(3) = 'Box1y'
     export_data(3) = position_y(2)

     !     data_headers(4) = 'Box2x'
     !     export_data(4) = position_x(3)
     !     data_headers(5) = 'Box2y'
     !     export_data(5) = position_y(3)

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ HAMILTONIAN CONSTRAINT ************!

     filename = "bondi-ham.con"
     num_cols = 8
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Haml ResN'
     export_data(2) = Ham_const_VolIntN*SymmFactor

     data_headers(3) = 'Haml ResD'
     export_data(3) = Ham_const_VolIntN*SymmFactor

     data_headers(4) = "Ham Res2P1.1N"
     export_data(4) = Ham_const_excised_VolIntN*SymmFactor

     data_headers(5) = "Ham Res2P1.1D"
     export_data(5) = Ham_const_excised_VolIntD*SymmFactor

     data_headers(6) = "HR1.1-NoD"
     export_data(6) = Ham_const_excised_VolIntN/Ham_const_excised_VolIntD

     data_headers(7) = "HRes2P1.1NNR"
     export_data(7) = Ham_const_excised_innerregion_VolIntN*SymmFactor

     data_headers(8) = "HRes2P1.1DNR"
     export_data(8) = Ham_const_excised_innerregion_VolIntD*SymmFactor



     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     !--------------------------------------------------------------------------------!
     !************ MOMENTUM CONSTRAINT ************!

     filename = 'bondi-mom.con'
     num_cols = 10
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = "Momx ResN"
     export_data(2) = momx_const_VolIntN*SymmFactor
     data_headers(3) = "Momy ResN"
     export_data(3) = momy_const_VolIntN*SymmFactor
     data_headers(4) = "Momz ResN"
     export_data(4) = momz_const_VolIntN*SymmFactor
     data_headers(5) = "Mom ResD"
     export_data(5) = mom_const_VolIntD*SymmFactor

     data_headers(6) = "Momx Res1.1N"
     export_data(6) = momx_const_excised_VolIntN*SymmFactor
     data_headers(7) = "Momy Res1.1N"
     export_data(7) = momy_const_excised_VolIntN*SymmFactor
     data_headers(8) = "Momz Res1.1N"
     export_data(8) = momz_const_excised_VolIntN*SymmFactor
     data_headers(9) = "Mom Res1.1D"
     export_data(9) = mom_const_excised_VolIntD*SymmFactor

     data_headers(10) = "MRx 1.1NoD"
     export_data(10) = momx_const_excised_VolIntN/mom_const_excised_VolIntD

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ GAMMA CONSTRAINT ************!
     filename = 'bondi-gam.con'
     num_cols = 4

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = "Gamx Res"
     export_data(2) = Gamx_const_VolInt*SymmFactor
     data_headers(3) = "Gamy Res"
     export_data(3) = Gamy_const_VolInt*SymmFactor
     data_headers(4) = "Gamz Res"
     export_data(4) = Gamz_const_VolInt*SymmFactor

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     ! BH local diagnostics
     ! ihf = Isolated horizon formalism
     filename = 'bondi-ihf-diags.mon'
     num_cols = 5
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     N_theta_bh_diagnostics = 50
     N_phi_bh_diagnostics = 50

     Mirr = 0.d0
     M_BH = 0.d0
     J_BH = 0.d0
     xi_err = 0.d0

     horizon_number = 1
     call BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, &
          horizon_number,N_theta_bh_diagnostics, &
          N_phi_bh_diagnostics,Symmetry,found_horizon)

     data_headers(2) = 'Mirr_bh1'
     export_data(2) = Mirr

     data_headers(3) = 'M_bh1'
     export_data(3) = M_BH

     data_headers(4) = 'J_bh1'
     export_data(4) = J_BH

     data_headers(5) = 'xi_err1'
     export_data(5) = xi_err

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ EM DIAGNOSTICS  ************!
     filename = 'bondi-emf.con'
     num_cols = 9
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     ! EM energy volume integral is computed in scheduler of diagnostics_mhd
     data_headers(2) = 'Eem'
     export_data(2) = em_energy_VolInt*SymmFactor

     ! EM energy volume integral is computed in scheduler of diagnostics_mhd
     data_headers(3) = 'Eem_outsideBH'
     export_data(3) = em_energy_outsideBH_VolInt*SymmFactor

     ! (Alt.) EM energy volume integral is computed in scheduler of diagnostics_mhd <-- CRAPPY, don't use
     !data_headers(3) = 'Eem2'
     !export_data(3) = em_energy2_VolInt*SymmFactor

     ! b_phi volume integral is computed in scheduler of diagnostics_mhd
     data_headers(4) = 'int B_phi'
     export_data(4) = b_phi_VolInt*SymmFactor

     data_headers(5) = 'Bx_max'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex (vindex, "mhd_evolve::Bx")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(5),1,vindex)

     data_headers(6) = 'By_max'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex (vindex, "mhd_evolve::By")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(6),1,vindex)

     data_headers(7) = "Monopole"
     export_data(7) = monopole_VolInt*SymmFactor

     data_headers(8) = "Monop.t=0"
     export_data(8) = initial_monopole_value*SymmFactor

     data_headers(9) = "Mpole_outBH"
     export_data(9) = monopole_outsideBH_VolInt*SymmFactor

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
  end if
end subroutine diagnostic_mag_bondi_accretion
