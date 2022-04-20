#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!------------------------------------
! Do stationary puncture diagnostics
!------------------------------------
subroutine bbh_kerr_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,ext_noghost
  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,Zmin
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20
  integer                                  :: vindex
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  !Needed for GW extraction:
  real*8                                   :: Mass_GW,hplus,hcross,moncrief_gw_extr_radius
  real*8                                   :: psi4r,psi4i
  integer                                  :: jj,nmodes_moncrief
  real*8,dimension(2)                      :: gw_amplitudes
  !End GW extraction variables  

  real*8                                   :: old_surf_radius,SymmFactor
  CCTK_POINTER :: reduction_value_pointer
  integer :: header_flag,handle
  integer :: index,num_cols
  integer :: ierr,myproc_rank,ii
  ! Variables used for BH local diagnotsics:
  integer :: found_horizon, horizon_number
  integer :: N_theta_bh_diagnostics, N_phi_bh_diagnostics
  real*8  :: Mirr, M_BH, J_BH, xi_err

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(trA_detg_enforce.ne.2) then
     !First perform sanity check.  Already done if trA_detg_enforce==2!
     call sanitycheck_restore_Aij(cctkGH,cctk_lsh,gxx,gxy,gxz,gyy,gyz,gzz, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
  end if

  if(MOD(cctk_iteration,out_every)==0) then
     call set_SymmFactor(Symmetry,SymmFactor)

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1
        ! Initialize variables for BH diagnostics
        Mirr1 = 0.d0
	Mirr2 = 0.d0
	Mirr3 = 0.d0
	M_BH1 = 0.d0
	M_BH2 = 0.d0
	M_BH3 = 0.d0
	J_BH1 = 0.d0
	J_BH2 = 0.d0
	J_BH3 = 0.d0
	xi_err1 = 0.d0
	xi_err2 = 0.d0
	xi_err3 = 0.d0
     else
        header_flag = 0
     end if

     !--------------------------------------------------------------------------------!

     filename = 'bbh_kerr.mon'
     num_cols = 11
     export_data = 0.D0
     data_headers = 'unknowndude'

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'lapsemin'
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(index,"lapse::lapm1")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(2),1,index)
     export_data(2) = export_data(2) + 1.D0

     data_headers(3) = 'phimax'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(index,"bssn::phi")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(3),1,index)

     data_headers(4) = 'Jsurf_BH1_inner'
     old_surf_radius = surf_radius
     surf_radius=BH_Vol_Excise_Radius
     call J_surf_integral_offcenter(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(4))
     surf_radius = old_surf_radius

     data_headers(5) = 'Jsurf_BH2_inner'
     old_surf_radius = surf_radius
     surf_radius=BH_Vol_Excise_Radius
     call J_surf_integral_offcenter(cctkGH,bh_posn_x(2),bh_posn_y(2),bh_posn_z(2),export_data(5))

     data_headers(6) = 'Jsurf_BH3_inner'
     old_surf_radius = surf_radius
     surf_radius=BH_Vol_Excise_Radius
     call J_surf_integral_offcenter(cctkGH,bh_posn_x(3),bh_posn_y(3),bh_posn_z(3),export_data(6))

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(7) = "J_exVol"
     export_data(7) = J_ADM_VolInt*SymmFactor

     data_headers(8) = "M_ADMsurf_BH1_inner"
     old_surf_radius = surf_radius
     surf_radius=BH_Vol_Excise_Radius
     call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(8))
     surf_radius = old_surf_radius

     data_headers(9) = "M_ADMsurf_BH2_inner"
     old_surf_radius = surf_radius
     surf_radius=BH_Vol_Excise_Radius
     call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(2),bh_posn_y(2),bh_posn_z(2),export_data(9))
     surf_radius = old_surf_radius

     data_headers(10) = "M_ADMsurf_BH3_inner"
     old_surf_radius = surf_radius
     surf_radius=BH_Vol_Excise_Radius
     call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(3),bh_posn_y(3),bh_posn_z(3),export_data(10))
     surf_radius = old_surf_radius

     ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(11) = "M_ADM_exVol"
     export_data(11) = M_ADM_VolInt*SymmFactor

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!

     filename = 'bbh_kerr.kon'
     num_cols = 10
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Xbh1'
     export_data(2) = bh_posn_x(1)

     data_headers(3) = 'Ybh1'
     export_data(3) = bh_posn_y(1)

     data_headers(4) = 'Zbh1'
     export_data(4) = bh_posn_z(1)

     data_headers(5) = 'Xbh2'
     export_data(5) = bh_posn_x(2)

     data_headers(6) = 'Ybh2'
     export_data(6) = bh_posn_y(2)

     data_headers(7) = 'Zbh2'
     export_data(7) = bh_posn_z(2)

     data_headers(8) = 'Xbh2'
     export_data(8) = bh_posn_x(3)

     data_headers(9) = 'Ybh2'
     export_data(9) = bh_posn_y(3)

     data_headers(10) = 'Zbh2'
     export_data(10) = bh_posn_z(3)

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ HAMILTONIAN CONSTRAINT ************!

     filename = 'bbh_kerr-ham.con'
     num_cols = 5
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Haml ResN'
     export_data(2) = Ham_const_VolIntN*SymmFactor

     data_headers(3) = 'Haml ResD'
     export_data(3) = Ham_const_VolIntD*SymmFactor

     data_headers(4) = "Ham Res2P1.1N"
     export_data(4) = Ham_const_excised_VolIntN*SymmFactor

     data_headers(5) = "Ham Res2P1.1D"
     export_data(5) = Ham_const_excised_VolIntD*SymmFactor

     if (enable_M_constraint==1) then 
        data_headers(6) = "M_constraint" 
        export_data(6) = M_constraint*SymmFactor
        data_headers(7) = "M_constraint_inner" 
        export_data(7) = M_constraint_inner*SymmFactor
	num_cols = 7
     end if

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ MOMENTUM CONSTRAINT ************!

     filename = 'bbh_kerr-mom.con'
     num_cols = 9
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

     if (enable_Jz_constraint==1) then 
        data_headers(num_cols+1) = "Jz_constraint"
        export_data(num_cols+1) = Jz_constraint*SymmFactor
        data_headers(num_cols+2) = "Jz_constraint_inner"
        export_data(num_cols+2) = Jz_constraint_inner*SymmFactor
        num_cols = num_cols+2
     end if

     if (enable_P_constraint==1) then
        data_headers(num_cols+1) = "P_constraint"
        export_data(num_cols+1) = P_constraint*SymmFactor
        data_headers(num_cols+2) = "P_constraint_inner"
        export_data(num_cols+2) = P_constraint_inner*SymmFactor
        num_cols = num_cols+2
     end if

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !*********** GAMMA CONSTRAINT ************!
     filename = 'bbh_kerr-gam.con'
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

     do ii=1,nsurf

        surf_radius=surfvec(ii)
        rbr=rbrvec(ii)
        drbdr=drbrvec(ii)
        ddrbddr=ddrbrvec(ii)

        if(nsurf.eq.1) then
           filename = 'bbh_kerr.jon'
        else
           write(filename,31)ii
31         FORMAT('bbh_kerr.jon.',I1)
        endif

        num_cols = 9

        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME

        data_headers(2) = "fisheye radius"
        export_data(2) = surf_radius

        data_headers(3) = "phys radius"
        export_data(3) = surf_radius*rbr

        data_headers(4) = "Mass_sur"
        call M_surf_integral(cctkGH,export_data(4))

        data_headers(5) = "Ang_mom_sur"
        call J_surf_integral_offcenter(cctkGH,0.D0,0.D0,0.D0,export_data(5))

        data_headers(6) = "Komar Mass"
        call Komar_surf_integralv2(cctkGH,export_data(6))

        data_headers(7) = "P_x2"
        data_headers(8) = "P_y2"
        data_headers(9) = "P_z2"
        call P_surf_integralv2(cctkGH,export_data(7),export_data(8),export_data(9))

        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end do


     !--------------------------------------------------------------------------------!
     ! BH local diagnostics
     !
     filename = 'bbh_kerr_BH_diagnostics.mon'
     num_cols = 16
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     N_theta_bh_diagnostics = 50
     N_phi_bh_diagnostics = 50

     horizon_number = 1
     call BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, &
                       horizon_number,N_theta_bh_diagnostics, &
                       N_phi_bh_diagnostics,Symmetry,found_horizon)

     ! Update Mirr1, M_BH1, J_BH1, and xi_err1 if found_horizon==1
     ! Otherwise, do not update them
     Mirr1 = Mirr1*(1-found_horizon) + Mirr*found_horizon
     M_BH1 = M_BH1*(1-found_horizon) + M_BH*found_horizon
     J_BH1 = J_BH1*(1-found_horizon) + J_BH*found_horizon
     xi_err1 = xi_err1*(1-found_horizon) + xi_err*found_horizon
     

     data_headers(2) = 'Mirr_bh1'
     export_data(2) = Mirr1

     data_headers(3) = 'M_bh1'
     export_data(3) = M_BH1

     data_headers(4) = 'J_bh1'
     export_data(4) = J_BH1

     data_headers(5) = 'xi_err1'
     export_data(5) = xi_err1
 
     data_headers(6) = 'found AH?'
     export_data(6) = found_horizon*1.d0

     horizon_number = 2
     call BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, &
                       horizon_number,N_theta_bh_diagnostics, &
                       N_phi_bh_diagnostics,Symmetry,found_horizon)

     ! Update Mirr2, M_BH2, J_BH2, and xi_err2 if found_horizon==1
     ! Otherwise, do not update them
     Mirr2 = Mirr2*(1-found_horizon) + Mirr*found_horizon
     M_BH2 = M_BH2*(1-found_horizon) + M_BH*found_horizon
     J_BH2 = J_BH2*(1-found_horizon) + J_BH*found_horizon
     xi_err2 = xi_err2*(1-found_horizon) + xi_err*found_horizon

     data_headers(7) = 'Mirr_bh2'
     export_data(7) = Mirr2

     data_headers(8) = 'M_bh2'
     export_data(8) = M_BH2

     data_headers(9) = 'J_bh2'
     export_data(9) = J_BH2

     data_headers(10) = 'xi_err2'
     export_data(10) = xi_err2

     data_headers(11) = 'found AH?'
     export_data(11) = found_horizon*1.d0

     ! Common horizon?
     horizon_number = 3
     call BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, &
                       horizon_number,N_theta_bh_diagnostics, &
                     N_phi_bh_diagnostics,Symmetry,found_horizon)

     ! Update Mirr3, M_BH3, J_BH3, and xi_err3 if found_horizon==1
     ! Otherwise, do not update them
     Mirr3 = Mirr3*(1-found_horizon) + Mirr*found_horizon
     M_BH3 = M_BH3*(1-found_horizon) + M_BH*found_horizon
     J_BH3 = J_BH3*(1-found_horizon) + J_BH*found_horizon
     xi_err3 = xi_err3*(1-found_horizon) + xi_err*found_horizon

     data_headers(12) = 'Mirr_bh3'
     export_data(12) = Mirr3

     data_headers(13) = 'M_bh3'
     export_data(13) = M_BH3

     data_headers(14) = 'J_bh3'
     export_data(14) = J_BH3

     data_headers(15) = 'xi_err3'
     export_data(15) = xi_err3

     data_headers(16) = 'found AH?'
     export_data(16) = found_horizon*1.d0

     if(CCTK_MyProc(CCTKGH)==0)  & 
       call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

  end if

end subroutine bbh_kerr_diagnostics
