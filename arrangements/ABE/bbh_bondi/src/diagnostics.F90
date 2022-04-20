#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!------------------------------------
! Do bbh_bondi diagnostics
!------------------------------------
subroutine bbh_bondi_diagnostics(CCTK_ARGUMENTS)

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
  real*8                                   :: F_E_fluid_BH,F_E_em_BH,F_M0_BH
  real*8                                   :: F_J_fluid_BH,F_J_em_BH
  real*8                                   :: F_E_fluid_R1,F_E_em_R1,F_M0_R1,F_J_fluid_R1,F_J_em_R1
  real*8                                   :: F_E_fluid_R2,F_E_em_R2,F_M0_R2,F_J_fluid_R2,F_J_em_R2
  real*8                                   :: F_E_fluid_R3,F_E_em_R3,F_M0_R3,F_J_fluid_R3,F_J_em_R3
  real*8                                   :: F_E_fluid_R4,F_E_em_R4,F_M0_R4,F_J_fluid_R4,F_J_em_R4
  real*8                                   :: F_E_fluid_R5,F_E_em_R5,F_M0_R5,F_J_fluid_R5,F_J_em_R5
  real*8                                   :: surf_area,ave_radius
  integer                                  :: found_horizon,horizon_number,foundflag

  ! Variables used for BH local diagnotsics:
  integer :: N_theta_bh_diagnostics, N_phi_bh_diagnostics
  real*8  :: Mirr, M_BH, J_BH, xi_err
  
  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(1),1,index)
  write(*,*) "HIIIIII zmax",export_data(1)
  
  if(MOD(cctk_iteration,out_every)==0) then
     call set_SymmFactor(Symmetry,SymmFactor)

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1
     else
        header_flag = 0
     end if

!!$     !--------------------------------------------------------------------------------!
!!$     !************ HAMILTONIAN CONSTRAINT ************!
!!$
!!$     filename = 'binary_bh-ham.con'
!!$     num_cols = 5
!!$     export_data = 0.D0
!!$
!!$     data_headers(1) = '# Time'
!!$     export_data(1) = CCTK_TIME
!!$
!!$     data_headers(2) = 'Haml ResN'
!!$     export_data(2) = Ham_const_VolIntN*SymmFactor
!!$
!!$     data_headers(3) = 'Haml ResD'
!!$     export_data(3) = Ham_const_VolIntD*SymmFactor
!!$
!!$     data_headers(4) = "Ham Res2P1.1N"
!!$     export_data(4) = Ham_const_excised_VolIntN*SymmFactor
!!$
!!$     data_headers(5) = "Ham Res2P1.1D"
!!$     export_data(5) = Ham_const_excised_VolIntD*SymmFactor
!!$
!!$     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

!!$     !--------------------------------------------------------------------------------!
!!$     !************ MOMENTUM CONSTRAINT ************!
!!$
!!$     filename = 'binary_bh-mom.con'
!!$     num_cols = 9
!!$     export_data = 0.D0
!!$
!!$     data_headers(1) = '# Time'
!!$     export_data(1) = CCTK_TIME
!!$
!!$     data_headers(2) = "Momx ResN"
!!$     export_data(2) = momx_const_VolIntN*SymmFactor
!!$     data_headers(3) = "Momy ResN"
!!$     export_data(3) = momy_const_VolIntN*SymmFactor
!!$     data_headers(4) = "Momz ResN"
!!$     export_data(4) = momz_const_VolIntN*SymmFactor
!!$     data_headers(5) = "Mom ResD"
!!$     export_data(5) = mom_const_VolIntD*SymmFactor
!!$
!!$     data_headers(6) = "Momx Res1.1N"
!!$     export_data(6) = momx_const_excised_VolIntN*SymmFactor
!!$     data_headers(7) = "Momy Res1.1N"
!!$     export_data(7) = momy_const_excised_VolIntN*SymmFactor
!!$     data_headers(8) = "Momz Res1.1N"
!!$     export_data(8) = momz_const_excised_VolIntN*SymmFactor
!!$     data_headers(9) = "Mom Res1.1D"
!!$     export_data(9) = mom_const_excised_VolIntD*SymmFactor
!!$
!!$     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

!!$     !--------------------------------------------------------------------------------!
!!$     !************ GAMMA CONSTRAINT ************!
!!$     filename = 'binary_bh-gam.con'
!!$     num_cols = 4
!!$
!!$     data_headers(1) = '# Time'
!!$     export_data(1) = CCTK_TIME
!!$
!!$     data_headers(2) = "Gamx Res"
!!$     export_data(2) = Gamx_const_VolInt*SymmFactor
!!$     data_headers(3) = "Gamy Res"
!!$     export_data(3) = Gamy_const_VolInt*SymmFactor
!!$     data_headers(4) = "Gamz Res"
!!$     export_data(4) = Gamz_const_VolInt*SymmFactor
!!$
!!$     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
!!$     !--------------------------------------------------------------------------------!
     filename = 'bondi.ron'
     num_cols = 16
     export_data = 0.D0
 
     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     ! M0 flux across 1st AH
     F_M0_BH = 0.d0
     F_E_fluid_BH = 0.d0
     F_E_em_BH = 0.d0
     F_J_fluid_BH = 0.d0
     F_J_em_BH = 0.d0
     
     found_horizon = 0
     horizon_number = 1
     if (cctk_iteration .ge. 0) then 
        if (bondi_metric_type .eq. 1) then
           call flux_through_analytic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH,&
                binary_orb_freq,xbh1_initial, &
                rah1_initial,CCTK_TIME,N_theta,N_phi,Symmetry)
        else if (bondi_metric_type .eq. 2) then
           call flux_through_analytic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH,&
                0.d0,-1.0, &
                rah1_initial,CCTK_TIME,N_theta,N_phi,Symmetry)
        else if (bondi_metric_type .eq. 3) then
           call flux_through_nonstatic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH, horizon_number, & 
                cctk_time,t_last_bbh_bondi,xh1_last_bbh_bondi,yh1_last_bbh_bondi,zh1_last_bbh_bondi, &
                ah1_radii_last_bbh_bondi,N_theta,N_phi,Symmetry,found_horizon)
           foundflag = HorizonLocalCoordinateOrigin(1,xh1_last_bbh_bondi,yh1_last_bbh_bondi,zh1_last_bbh_bondi)   
        endif
     end if
     data_headers(2) = "F_M0_BH1" 
     export_data(2) = F_M0_BH
     data_headers(3) = "F_E_fluid_BH1"
     export_data(3) = F_E_fluid_BH
     data_headers(4) = "F_E_em_BH1"
     export_data(4) = F_E_em_BH 
     data_headers(5) = "F_J_fluid_BH1"
     export_data(5) = F_J_fluid_BH
     data_headers(6) = "F_J_em_BH1"
     export_data(6) = F_J_em_BH 
     
     ! M0 flux across 2nd AH
     F_M0_BH = 0.d0
     F_E_fluid_BH = 0.d0
     F_E_em_BH = 0.d0
     F_J_fluid_BH = 0.d0
     F_J_em_BH = 0.d0
     
     found_horizon = 0
     horizon_number = 2
     if (cctk_iteration .ge. 0) then
        if (bondi_metric_type .eq. 1) then
           call flux_through_analytic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH,&
                binary_orb_freq,xbh2_initial,&
                rah2_initial,CCTK_TIME,N_theta,N_phi,Symmetry)
        else if (bondi_metric_type .eq. 3) then
           call flux_through_nonstatic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH, horizon_number, &
                cctk_time,t_last_bbh_bondi,xh2_last_bbh_bondi,yh2_last_bbh_bondi,zh2_last_bbh_bondi, & 
                ah2_radii_last_bbh_bondi,N_theta,N_phi,Symmetry,found_horizon) 
           foundflag = HorizonLocalCoordinateOrigin(2,xh2_last_bbh_bondi,yh2_last_bbh_bondi,zh2_last_bbh_bondi)
        endif
     end if
     
     data_headers(7) = "F_M0_BH2"
     export_data(7) = F_M0_BH
     data_headers(8) = "F_E_fluid_BH2"
     export_data(8) = F_E_fluid_BH 
     data_headers(9) = "F_E_em_BH2"
     export_data(9) = F_E_em_BH 
     data_headers(10) = "F_J_fluid_BH2"
     export_data(10) = F_J_fluid_BH 
     data_headers(11) = "F_J_em_BH2"
     export_data(11) = F_J_em_BH 

     
     ! M0 flux across 3rd AH
     F_M0_BH = 0.d0
     F_E_fluid_BH = 0.d0
     F_E_em_BH = 0.d0
     found_horizon = 0
     horizon_number = 3
     if (cctk_iteration .ge. 0) then 
        if (bondi_metric_type .eq. 3) then
           call flux_through_nonstatic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH, horizon_number, &
                cctk_time,t_last_bbh_bondi,xh3_last_bbh_bondi,yh3_last_bbh_bondi,zh3_last_bbh_bondi, &
                ah3_radii_last_bbh_bondi,N_theta,N_phi,Symmetry,found_horizon)
           foundflag = HorizonLocalCoordinateOrigin(3,xh3_last_bbh_bondi,yh3_last_bbh_bondi,zh3_last_bbh_bondi) 
        endif
     end if
     data_headers(12) = "F_M0_BH2"
     export_data(12) = F_M0_BH
     data_headers(13) = "F_E_fluid_BH2"
     export_data(13) = F_E_fluid_BH
     data_headers(14) = "F_E_em_BH2"
     export_data(14) = F_E_em_BH
     data_headers(15) = "F_J_fluid_BH2"
     export_data(15) = F_J_fluid_BH
     data_headers(16) = "F_J_em_BH2"
     export_data(16) = F_J_em_BH

     
     write(*,*) "end of .ron"
     
     t_last_bbh_bondi=cctk_time

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
      filename = 'bondi.son'
     num_cols = 16
!     num_cols = 10
     export_data = 0.D0
 
     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     F_M0_R1 = 0.d0
     F_E_fluid_R1 = 0.d0 
     F_E_em_R1 = 0.d0 
     F_J_fluid_R1 = 0.d0 
     F_J_em_R1 = 0.d0 
     
     call flux_through_analytic_BH(cctkGH,F_E_fluid_R1,F_E_em_R1,F_M0_R1,F_J_fluid_R1,F_J_em_R1,&
          0.d0,0.d0, &
          6.d0,CCTK_TIME,N_theta,N_phi,Symmetry)
     
     data_headers(2) = "F_M0_R1"
     export_data(2) = F_M0_R1
     data_headers(3) = "F_E_fluid_R1"
     export_data(3) = F_E_fluid_R1 
     data_headers(4) = "F_J_fluid_R1"
     export_data(4) = F_J_fluid_R1 
     
     F_M0_R2 = 0.d0
     F_E_fluid_R2 = 0.d0 
     F_E_em_R2 = 0.d0 
     F_J_fluid_R2 = 0.d0 
     F_J_em_R2 = 0.d0 
     
     call flux_through_analytic_BH(cctkGH,F_E_fluid_R2,F_E_em_R2,F_M0_R2,F_J_fluid_R2,F_J_em_R2,&
          0.d0,0.d0, &
          10.d0,CCTK_TIME,N_theta,N_phi,Symmetry)
     
     data_headers(5) = "F_M0_R2"
     export_data(5) = F_M0_R2
     data_headers(6) = "F_E_fluid_R2"
     export_data(6) = F_E_fluid_R2 
     data_headers(7) = "F_J_fluid_R2"
     export_data(7) = F_J_fluid_R2 
     
     
     F_M0_R3 = 0.d0
     F_E_fluid_R3 = 0.d0 
     F_E_em_R3 = 0.d0 
     F_J_fluid_R3 = 0.d0 
     F_J_em_R3 = 0.d0 
     
     call flux_through_analytic_BH(cctkGH,F_E_fluid_R3,F_E_em_R3,F_M0_R3,F_J_fluid_R3,F_J_em_R3,&
          0.d0,0.d0, &
          20.d0,CCTK_TIME,N_theta,N_phi,Symmetry)
     
     data_headers(8) = "F_M0_R3"
     export_data(8) = F_M0_R3
     data_headers(9) = "F_E_fluid_R3"
     export_data(9) = F_E_fluid_R3 
     data_headers(10) = "F_J_fluid_R3"
     export_data(10) = F_J_fluid_R3 
     
     F_M0_R4 = 0.d0
     F_E_fluid_R4 = 0.d0 
     F_E_em_R4 = 0.d0 
     F_J_fluid_R4 = 0.d0 
     F_J_em_R4 = 0.d0 
     
     call flux_through_analytic_BH(cctkGH,F_E_fluid_R4,F_E_em_R4,F_M0_R4,F_J_fluid_R4,F_J_em_R4,&
          0.d0,0.d0, &
          50.d0,CCTK_TIME,N_theta,N_phi,Symmetry)
     
     data_headers(11) = "F_M0_R4"
     export_data(11) = F_M0_R4
     data_headers(12) = "F_E_fluid_R4"
     export_data(12) = F_E_fluid_R4 
     data_headers(13) = "F_J_fluid_R4"
     export_data(13) = F_J_fluid_R4 
     
     F_M0_R5 = 0.d0
     F_E_fluid_R5 = 0.d0 
     F_E_em_R5 = 0.d0 
     F_J_fluid_R5 = 0.d0 
     F_J_em_R5 = 0.d0 
     
     call flux_through_analytic_BH(cctkGH,F_E_fluid_R5,F_E_em_R5,F_M0_R5,F_J_fluid_R5,F_J_em_R5,&
          0.d0,0.d0, &
          100.d0,CCTK_TIME,N_theta,N_phi,Symmetry)
     
     data_headers(14) = "F_M0_R5"
     export_data(14) = F_M0_R5
     data_headers(15) = "F_E_fluid_R5"
     export_data(15) = F_E_fluid_R5 
     data_headers(16) = "F_J_fluid_R5"
     export_data(16) = F_J_fluid_R5 
     
     

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!



     
     do ii=1,nsurf

        surf_radius=surfvec(ii)
        rbr=rbrvec(ii)
        drbdr=drbrvec(ii)
        ddrbddr=ddrbrvec(ii)

        if(nsurf.eq.1) then
           filename = 'bondi.jon'
        else
           write(filename,21)ii
21         FORMAT("bondi.jon.",I1)
        endif

        num_cols = 13

        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME

        data_headers(2) = "fisheye radius"
        export_data(2) = surf_radius

        data_headers(3) = "phys radius"
        export_data(3) = surf_radius*rbr

        data_headers(4) = "Mass_sur"
        call M_surf_integral(cctkGH,export_data(4))

        data_headers(5) = "Ang_mom_surf"
        call J_surf_integral_offcenter(cctkGH,0.D0,0.D0,0.D0,export_data(5))

        data_headers(6) = "Komar Mass"
        call Komar_surf_integralv2(cctkGH,export_data(6))

        data_headers(7) = "P_x2"
        data_headers(8) = "P_y2"
        data_headers(9) = "P_z2"
        call P_surf_integralv2(cctkGH,export_data(7),export_data(8),export_data(9))
	data_headers(10) = "M_surf_mov"
	call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(10))

        data_headers(11) = "P_x"
        data_headers(12) = "P_y"
        data_headers(13) = "P_z"
        call P_surf_integral_offcenter(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(11),export_data(12),export_data(13))


        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end do


 !--------------------------------------------------------------------------------!
     ! BH local diagnostics
     !
     filename = 'binary_bh_BH_diagnostics.mon'
     num_cols = 13
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

     horizon_number = 2
     call BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, &
                       horizon_number,N_theta_bh_diagnostics, & 
        	       N_phi_bh_diagnostics,Symmetry,found_horizon)

     data_headers(6) = 'Mirr_bh2'
     export_data(6) = Mirr

     data_headers(7) = 'M_bh2'
     export_data(7) = M_BH

     data_headers(8) = 'J_bh2'
     export_data(8) = J_BH

     data_headers(9) = 'xi_err2'
     export_data(9) = xi_err

     ! Common horizon?
     horizon_number = 3
     call BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, &
                       horizon_number,N_theta_bh_diagnostics, & 
                     N_phi_bh_diagnostics,Symmetry,found_horizon)

     data_headers(10) = 'Mirr_bh3'
     export_data(10) = Mirr

     data_headers(11) = 'M_bh3'
     export_data(11) = M_BH

     data_headers(12) = 'J_bh3'
     export_data(12) = J_BH

     data_headers(13) = 'xi_err3'
     export_data(13) = xi_err
     
     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------
     filename = 'bondi.lum'
     num_cols = 4
     export_data = 0.d0
     data_headers = 'unknowndud'

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME
     
     data_headers(2) = 'brem_qei'
     export_data(2) = brem_qei_VolInt*SymmFactor
     
     data_headers(3) = 'brem_qee'
     export_data(3) = brem_qee_VolInt*SymmFactor
     
     data_headers(4) = 'synchrotron'
     export_data(4) = synch_VolInt*SymmFactor                                                                                                  
       

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     !--------------------------------------------------------------------------------
     
  end if

end subroutine bbh_bondi_diagnostics
