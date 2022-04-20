#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!------------------------------------
! Do disk_moncrief diagnostics
!------------------------------------
subroutine disk_moncrief_diagnostics(CCTK_ARGUMENTS)

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
  real*8                                   :: F_M0_BH,surf_area,ave_radius
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
     filename = 'bondi.ron'
  
!     num_cols = 12
     num_cols = 4
     export_data = 0.D0
 
     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     ! M0 flux across 1st AH
     F_M0_BH = 0.d0
     found_horizon = 0
     horizon_number = 1
   
     if (cctk_iteration .gt. 0) then 
        call M0_flux_through_cookpfeiffer_rot_BH(cctkGH,F_M0_BH,horizon_number,binary_orb_freq,binary_separation, &
             ah_rad,CCTK_TIME,N_theta,N_phi,Symmetry)
     end if
     data_headers(2) = "F_M0_BH1" 
     export_data(2) = F_M0_BH
     
     ! M0 flux across 2nd AH
     F_M0_BH = 0.d0
     found_horizon = 0
     horizon_number = 2
     if (cctk_iteration .gt. 0) then
        call M0_flux_through_cookpfeiffer_rot_BH(cctkGH,F_M0_BH,horizon_number,binary_orb_freq,binary_separation, &
             ah_rad,CCTK_TIME,N_theta,N_phi,Symmetry)
     end if
     data_headers(3) = "F_M0_BH2"
     export_data(3) = F_M0_BH
     
     if (1.eq.0) then
        ! M0 flux across 3rd AH
        F_M0_BH = 0.d0
        found_horizon = 0
        horizon_number = 3
        if (cctk_iteration .gt. 0) then 
           call M0_flux_through_cookpfeiffer_rot_BH(cctkGH,F_M0_BH,horizon_number,binary_orb_freq,binary_separation, &
                ah_rad,CCTK_TIME,N_theta,N_phi,Symmetry)
        end if
        data_headers(4) = "F_M0_BH3" 
        export_data(4) = F_M0_BH
     else
        data_headers(4) = "F_M0_BH3" 
        export_data(4) = 0.d0
     end if
     
     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     
     !--------------------------------------------------------------------------------!
   
!!$     !--------------------------------------------------------------------------------
!!$     filename = 'bondi.lum'
!!$     num_cols = 4
!!$     export_data = 0.d0
!!$     data_headers = 'unknowndud'
!!$
!!$     data_headers(1) = '# Time'
!!$     export_data(1) = CCTK_TIME
!!$     
!!$     data_headers(2) = 'brem_qei'
!!$     export_data(2) = brem_qei_VolInt*SymmFactor
!!$     
!!$     data_headers(3) = 'brem_qee'
!!$     export_data(3) = brem_qee_VolInt*SymmFactor
!!$     
!!$     data_headers(4) = 'synchrotron'
!!$     export_data(4) = synch_VolInt*SymmFactor                                                                                                  
!!$       
!!$
!!$     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
!!$     !--------------------------------------------------------------------------------
!!$     
  end if

end subroutine disk_moncrief_diagnostics
