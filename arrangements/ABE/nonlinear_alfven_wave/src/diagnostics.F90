#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for shocktests thorn
!-----------------------------------------------------------------------------
subroutine nonlinear_alfven_wave_diagnostics(CCTK_ARGUMENTS)

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

     !!filename = 's.mon'

     !!num_cols = 3
     !!export_data = 0.D0

     !!data_headers(1) = '# Time'
     !!export_data(1) = CCTK_TIME

     !!data_headers(2) = 'Mo 1D VolInt'
     !!export_data(2) = shock_restmass_VolInt

     !!if(cctk_iteration==0) M0_initial=shock_restmass_VolInt

     !!data_headers(3) = '|1-M0/M0(t=0)|'
     !!export_data(3) = abs(1.D0-shock_restmass_VolInt/M0_initial)

     !!
     !!if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
  end if
end subroutine nonlinear_alfven_wave_diagnostics
