#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for shocktests thorn
!-----------------------------------------------------------------------------
subroutine shocktests_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,Zmin
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20

  integer :: num_rows
  parameter(num_rows = 50)
  real*8, dimension(num_rows)              :: R_out,rho_out,P_out,v_out
  integer                                  :: vindex
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  real*8                                   :: old_surf_radius,SymmFactor,gammashock,xmin,xmax
  real*8                                   :: rhoshock_out
  CCTK_POINTER :: reduction_value_pointer
  integer :: header_flag,handle
  integer :: index,num_cols
  integer :: ierr,myproc_rank,ii

  ext = cctk_lsh

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

     filename = 'stests.mon'
     num_cols = 4
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'rho_err VolInt'
     export_data(2) = rho_error_VolInt

     data_headers(3) = 'P_err VolInt'
     export_data(3) = P_error_VolInt

     data_headers(4) = 'vx_err VolInt'
     export_data(4) = vx_error_VolInt

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     filename = 'stests.don'

     num_cols = 5
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Mo 1D VolInt'
     export_data(2) = shock_restmass_VolInt

     data_headers(3) = 'Mo Analytic'
     export_data(3) = shock_restmass_analytic

     data_headers(4) = '|1-M0/M0_ana(t=0)|'
     export_data(4) = abs(1.D0-shock_restmass_VolInt/shock_restmass_analytic_initial)

     data_headers(5) = '|1-M0/M0(t=0)|'
     if(cctk_iteration==0) M0_initial=shock_restmass_VolInt
     export_data(5) = abs(1.D0-shock_restmass_VolInt/M0_initial)

     
     write(*,*) "HIEEEE",shock_restmass_analytic_initial,shock_restmass_VolInt

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     filename = 'stests.ron'
     num_cols = 4
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'rho_err reflect'
     export_data(2) = rho_error_VolInt_reflection

     data_headers(3) = 'P_err reflect'
     export_data(3) = P_error_VolInt_reflection

     data_headers(4) = 'vx_err reflect'
     export_data(4) = vx_error_VolInt_reflection

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     !--------------------------------------------------------------------------------!
  end if
end subroutine shocktests_diagnostics
