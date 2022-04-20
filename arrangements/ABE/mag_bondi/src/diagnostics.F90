!------------------------------------------------
! Diagnostics driver routine for disk evolutions
!------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine diagnostic_mag_bondi_accretion(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)    :: ext
  real*8, dimension(20)    :: export_data
  character, dimension(20) :: data_headers*20
  character                :: filename*50
  real*8                   :: Mass0,Jtot,SymmFactor,temp
  real*8                   :: dIntegral,multfactor,dIntegral_Mon
  real*8                   :: F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em
  real*8                   :: out_M0_flux1, out_fluid_energy_flux1
  real*8                   :: out_EM_energy_flux1, out_fluid_J_flux1
  real*8                   :: out_EM_J_flux1, out_M0_flux2,out_fluid_energy_flux2
  real*8                   :: out_EM_energy_flux2, out_fluid_J_flux2, out_EM_J_flux2
  real*8                   :: out_M0_flux3, out_fluid_energy_flux3, out_EM_energy_flux3
  real*8                   :: out_fluid_J_flux3, out_EM_J_flux3
  CCTK_REAL :: reduction_value
  CCTK_INT :: red_tmp
  integer :: header_flag,handle,dummy,vindex
  integer :: index,num_cols,i
  integer :: ierr,myproc_rank
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  integer, parameter :: AXISYM_FULL = 5
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  !--------------------------------------------------------------------------------!
  ! Update movingbox center
  !---------------------------------------------------------------------------------
  if (enable_movingbox==1) then
     position_x(1) = 2.d0*sin(0.1d0*cctk_time)
     position_y(1) = 1.6d0*(1.d0-cos(0.1d0*cctk_time))
     position_z(1) = 0.d0
  end if

  if (mod(cctk_iteration,out_every) .eq. 0) then 

     !header_flag == 1 -> write file headers
     if(cctk_iteration==0) then
        header_flag = 1
     else
        header_flag = 0
     end if

     call set_SymmFactor(Symmetry,SymmFactor)

     !--------------------------------------------------------------------------------!
     ! Calculations for the mon file.
     !--------------------------------------------------------------------------------!

     filename = 'mbondi_accretion.don'
     num_cols = 3
     export_data = 0.D0
     data_headers = 'unknowndude'

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME
     
     data_headers(2) = 'rhocent'
     call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
     call interp_center_carp(cctkGH,vindex,export_data(2))

     data_headers(3) = 'rho_b_max'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(3),1,vindex)

     if(CCTK_MyProc(CCTKGH) .eq. 0.D0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if


     !--------------------------------------------------------------------------------!
     ! Calculations for the mon file.
     !--------------------------------------------------------------------------------!

     call mag_bondi_accretion_fluxes_ana(cctkGH, Symmetry, N_theta, N_phi, & 
          xn_d, yn_d, zn_d, nx_d,ny_d,nz_d, &
          F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em)

     filename = 'mbondi_accretion.mon'
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     i = 2
     data_headers(i) = 'F_M0'
     export_data(i) = F_M0

     i = i+1
     data_headers(i) = 'F_E_fluid'
     export_data(i) = F_E_fluid

     i = i+1
     data_headers(i) = 'F_E_em'
     export_data(i) = F_E_em

     i = i+1
     data_headers(i) = 'F_J_fluid'
     export_data(i) = F_J_fluid

     i = i+1
     data_headers(i) = 'F_J_em'
     export_data(i) = F_J_em

     i = i+1  
     data_headers(i) = 'Monpole'
     export_data(i) = monopole_VolInt*SymmFactor

     num_cols = i

     if(CCTK_MyProc(CCTKGH) .eq. 0.D0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !--------------------------------------------------------------------------------!
     ! Calculations for the fon file
     !--------------------------------------------------------------------------------!

     filename = 'mbondi_accretion.fon'

     ! Calculate the fluxes through the various shells.
     temp = surf_radius
     surf_radius = r_out_flux1
     call mhd_surface_fluxes(cctkGH,F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em)
     out_M0_flux1           = F_M0
     out_fluid_energy_flux1 = F_E_fluid
     out_EM_energy_flux1    = F_E_em
     out_fluid_J_flux1      = F_J_fluid
     out_EM_J_flux1         = F_J_em

     surf_radius = r_out_flux2
     call mhd_surface_fluxes(cctkGH,F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em)
     out_M0_flux2           = F_M0
     out_fluid_energy_flux2 = F_E_fluid
     out_EM_energy_flux2    = F_E_em
     out_fluid_J_flux2      = F_J_fluid
     out_EM_J_flux2         = F_J_em

     surf_radius = r_out_flux3
     call mhd_surface_fluxes(cctkGH,F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em)
     out_M0_flux3           = F_M0
     out_fluid_energy_flux3 = F_E_fluid
     out_EM_energy_flux3    = F_E_em
     out_fluid_J_flux3      = F_J_fluid
     out_EM_J_flux3         = F_J_em

     num_cols = 16
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'fM0_1'
     export_data(2) = out_M0_flux1
     data_headers(3) = 'fE_fluid_1'
     export_data(3) =  out_fluid_energy_flux1
     data_headers(4) = 'fE_em_1'
     export_data(4) = out_EM_energy_flux1
     data_headers(5) = 'fJ_fluid_1' 
     export_data(5) = out_fluid_J_flux1
     data_headers(6) = 'fJ_em_1'
     export_data(6) = out_EM_J_flux1

     data_headers(7) = 'fM0_2'
     export_data(7) = out_M0_flux2
     data_headers(8) = 'fE_fluid_2'
     export_data(8) =  out_fluid_energy_flux2
     data_headers(9) = 'fE_em_2'
     export_data(9) = out_EM_energy_flux2
     data_headers(10) = 'fJ_fluid_2' 
     export_data(10) = out_fluid_J_flux2
     data_headers(11) = 'fJ_em_2'
     export_data(11) = out_EM_J_flux2

     data_headers(12) = 'fM0_3'
     export_data(12) = out_M0_flux3
     data_headers(13) = 'fE_fluid_3'
     export_data(13) =  out_fluid_energy_flux3
     data_headers(14) = 'fE_em_3'
     export_data(14) = out_EM_energy_flux3
     data_headers(15) = 'fJ_fluid_3' 
     export_data(15) = out_fluid_J_flux3
     data_headers(16) = 'fJ_em_3'
     export_data(16) = out_EM_J_flux3

     surf_radius = temp

     if(CCTK_MyProc(CCTKGH) .eq. 0.D0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !--------------------------------------------------------------------------------!
     filename = 'mbondi_accretion.xon'
     num_cols = 3
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Box1x'
     export_data(2) = position_x(1) 
     data_headers(3) = 'Box1y'
     export_data(3) = position_y(1)

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

  end if

end subroutine diagnostic_mag_bondi_accretion
