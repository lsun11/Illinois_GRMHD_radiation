#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!------------------------------------
! Do disk_powerlaw diagnostics
!------------------------------------
subroutine disk_powerlaw_diagnostics(CCTK_ARGUMENTS)

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


     !--------------------------------------------------------------------------------!
     filename = 'disk_powerlaw.ron'
     write(*,*) "beginning of .ron"
     num_cols = 16
!     num_cols = 10
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
        if (metric_type .eq. 1) then
           call flux_through_analytic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH,&
                binary_orb_freq,xbh1_initial, &
                rah1_initial,CCTK_TIME,N_theta,N_phi,Symmetry)
        else if (metric_type .eq. 2) then
           call flux_through_analytic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH,&
                0.d0,-1.0, &
                rah1_initial,CCTK_TIME,N_theta,N_phi,Symmetry)
        else if (metric_type .eq. 3) then
           call flux_through_nonstatic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH, horizon_number, & 
                cctk_time,t_last_disk_powerlaw,xh1_last_disk_powerlaw,yh1_last_disk_powerlaw,zh1_last_disk_powerlaw, &
                ah1_radii_last_disk_powerlaw,N_theta,N_phi,Symmetry,found_horizon)
           foundflag = HorizonLocalCoordinateOrigin(1,xh1_last_disk_powerlaw,yh1_last_disk_powerlaw,zh1_last_disk_powerlaw)   
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
        if (metric_type .eq. 1) then
           call flux_through_analytic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH,&
                binary_orb_freq,xbh2_initial,&
                rah2_initial,CCTK_TIME,N_theta,N_phi,Symmetry)
        else if (metric_type .eq. 3) then
           call flux_through_nonstatic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH, horizon_number, &
                cctk_time,t_last_disk_powerlaw,xh2_last_disk_powerlaw,yh2_last_disk_powerlaw,zh2_last_disk_powerlaw, & 
                ah2_radii_last_disk_powerlaw,N_theta,N_phi,Symmetry,found_horizon) 
           foundflag = HorizonLocalCoordinateOrigin(2,xh2_last_disk_powerlaw,yh2_last_disk_powerlaw,zh2_last_disk_powerlaw)
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
     F_J_fluid_BH = 0.d0
     F_J_em_BH = 0.d0
     
     found_horizon = 0
     horizon_number = 3
     if (cctk_iteration .ge. 0) then 
        if (metric_type .eq. 3) then
           call flux_through_nonstatic_BH(cctkGH,F_E_fluid_BH,F_E_em_BH,F_M0_BH,F_J_fluid_BH,F_J_em_BH, horizon_number, &
                cctk_time,t_last_disk_powerlaw,xh3_last_disk_powerlaw,yh3_last_disk_powerlaw,zh3_last_disk_powerlaw, & 
                ah3_radii_last_disk_powerlaw,N_theta,N_phi,Symmetry,found_horizon) 
           foundflag = HorizonLocalCoordinateOrigin(3,xh3_last_disk_powerlaw,yh3_last_disk_powerlaw,zh3_last_disk_powerlaw) 
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
     
     t_last_disk_powerlaw=cctk_time

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     
     !--------------------------------------------------------------------------------!


     filename = 'disk_powerlaw.son'
!     num_cols = 16
     num_cols = 10
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
          10.d0,CCTK_TIME,N_theta,N_phi,Symmetry)
     
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
          25.d0,CCTK_TIME,N_theta,N_phi,Symmetry)
     
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
          40.d0,CCTK_TIME,N_theta,N_phi,Symmetry)
     
     data_headers(8) = "F_M0_R3"
     export_data(8) = F_M0_R3
     data_headers(9) = "F_E_fluid_R3"
     export_data(9) = F_E_fluid_R3 
     data_headers(10) = "F_J_fluid_R3"
     export_data(10) = F_J_fluid_R3 
     
!     F_M0_R4 = 0.d0
!     F_E_fluid_R4 = 0.d0 
!     F_E_em_R4 = 0.d0 
!     F_J_fluid_R4 = 0.d0 
!     F_J_em_R4 = 0.d0 
     
!     call flux_through_analytic_BH(cctkGH,F_E_fluid_R4,F_E_em_R4,F_M0_R4,F_J_fluid_R4,F_J_em_R4,&
!          0.d0,0.d0, &
!          20.d0,CCTK_TIME,N_theta,N_phi,Symmetry)
     
!     data_headers(11) = "F_M0_R4"
!     export_data(11) = F_M0_R4
!     data_headers(12) = "F_E_fluid_R4"
!     export_data(12) = F_E_fluid_R4 
!     data_headers(13) = "F_J_fluid_R4"
!     export_data(13) = F_J_fluid_R4 
     
!     F_M0_R5 = 0.d0
!     F_E_fluid_R5 = 0.d0 
!     F_E_em_R5 = 0.d0 
!     F_J_fluid_R5 = 0.d0 
!     F_J_em_R5 = 0.d0 
     
 !    call flux_through_analytic_BH(cctkGH,F_E_fluid_R5,F_E_em_R5,F_M0_R5,F_J_fluid_R5,F_J_em_R5,&
 !         0.d0,0.d0, &
 !         25.d0,CCTK_TIME,N_theta,N_phi,Symmetry)
     
 !    data_headers(14) = "F_M0_R5"
 !    export_data(14) = F_M0_R5
 !    data_headers(15) = "F_E_fluid_R5"
 !    export_data(15) = F_E_fluid_R5 
 !    data_headers(16) = "F_J_fluid_R5"
 !    export_data(16) = F_J_fluid_R5 
     
     

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!



     filename='disk_powerlaw.don'
     
     num_cols = 24
     export_data = 0.D0
     
     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME
     
     ! M0 volume integral is computed in scheduler of diagnostics_mhd
     data_headers(2) = "M0"
     export_data(2) = M0_VolInt*SymmFactor 
     data_headers(3) = "M0_horiz"
     export_data(3) = M0_horiz_VolInt*SymmFactor
     data_headers(4) = "M0_r1"
     export_data(4) = M0_r1_VolInt*SymmFactor
     data_headers(5) = "M0_r2"
     export_data(5) = M0_r2_VolInt*SymmFactor
     data_headers(6) = "M0_r3"
     export_data(6) = M0_r3_VolInt*SymmFactor
     data_headers(7) = "fluid_energy_horiz"
     export_data(7) = fluid_energy_horiz_VolInt*SymmFactor
     data_headers(8) = "fluid_energy_r1"
     export_data(8) = fluid_energy_r1_VolInt*SymmFactor
     data_headers(9) = "fluid_energy_r2"
     export_data(9) = fluid_energy_r2_VolInt*SymmFactor
     data_headers(10) = "fluid_energy_r3"
     export_data(10) = fluid_energy_r3_VolInt*SymmFactor

!     data_headers(5) = "em_energy_between"
!     export_data(5) = em_energy_between_VolInt*SymmFactor
     data_headers(11) = "fluid_J_horiz"
     export_data(11) = fluid_J_horiz_VolInt*SymmFactor
     data_headers(12) = "fluid_J_r1"
     export_data(12) = fluid_J_r1_VolInt*SymmFactor
     data_headers(13) = "fluid_J_r2"
     export_data(13) = fluid_J_r2_VolInt*SymmFactor
     data_headers(14) = "fluid_J_r3"
     export_data(14) = fluid_J_r3_VolInt*SymmFactor

     data_headers(15) = "minternal_horiz"
     export_data(15) = minternal_horiz_VolInt*SymmFactor
     data_headers(16) = "minternal_r1" 
     export_data(16) = minternal_r1_VolInt*SymmFactor 
     data_headers(17) = "minternal_r2" 
     export_data(17) = minternal_r2_VolInt*SymmFactor 
     data_headers(18) = "minternal_r3"
     export_data(18) = minternal_r3_VolInt*SymmFactor 
     
     data_headers(19) = "minternal_cold_horiz"
     export_data(19) = minternal_cold_horiz_VolInt*SymmFactor 
     data_headers(20) = "minternal_cold_r1"
     export_data(20) = minternal_cold_r1_VolInt*SymmFactor 
     data_headers(21) = "minternal_cold_r2"
     export_data(21) = minternal_cold_r2_VolInt*SymmFactor 
     data_headers(22) = "minternal_cold_r3"
     export_data(22) = minternal_cold_r3_VolInt*SymmFactor    
     
     data_headers(23) = "minternal"
     export_data(23) = minternal_VolInt*SymmFactor 
     data_headers(24) = "minternal_cold" 
     export_data(24) = minternal_cold_VolInt*SymmFactor

!     data_headers(7) = "em_J_between"
!     export_data(7) = em_J_between_VolInt*SymmFactor
     
     
     !     ! M0 inside AH volume integral is computed in scheduler of diagnostics_mhd
     !     data_headers(3) = "M0 inside AH"
     !     export_data(3) = M0_AH_VolInt*SymmFactor
     
     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     
     !--------------------------------------------------------------------------------
     if((MOD(cctk_iteration,16*out_every)==0) .and. (1==0)) then!!! .and.(CCTK_MyProc(CCTKGH)==0)) then
        !call surface_dens_1d(cctkGH,CCTK_TIME,binary_orb_freq,surfxnum,surfphinum,surfznum,surfzmax,surfxmin,surfxmax,Symmetry,CCTK_MyProc(CCTKGH))
        !call surface_dens_2d(cctkGH,CCTK_TIME,surfxnum,surfznum,surfxmax,surfzmax,Symmetry,CCTK_MyProc(CCTKGH))
        call surface_dens_1d(cctkGH,CCTK_TIME,binary_orb_freq,surfxnum,surfphinum,surfznum,surfzmax,surfxmin,surfxmax,Symmetry,CCTK_MyProc(CCTKGH))
        call surface_dens_2d(cctkGH,CCTK_TIME,surfxnum,surfznum,surfxmax,surfzmax,Symmetry,CCTK_MyProc(CCTKGH))  
     endif
     !--------------------------------------------------------------------------------

     write(*,*) "before lum"
        
     !--------------------------------------------------------------------------------
     filename = 'disk_powerlaw.lum'
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
     write(*,*) "after lum"
  end if

end subroutine disk_powerlaw_diagnostics
