!-----------------------------------------------------------------------------
!
!$Id: predict_mhd.F90  $
!
!-----------------------------------------------------------------------------
!
! Predict & Correct RHS routines for mhd variables
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
!
! Predictor
!
!-----------------------------------------------------------------------------
subroutine bondi_diag_integrals(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3)                    :: ext,ext_noghost
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20
  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,xmax2
  real*8                                   :: RestMass,ADMMass,Jtot,temp,resid,resid_norm
  real*8                                   :: Mass_GW,hpx_alt
  real*8,dimension(2)                      :: gw_amplitudes
  real*8                                   :: resid_x,resid_y,resid_z,resid_mom_norm
  real*8                                   :: dIntegral
  real*8                                   :: gxres,gyres,gzres
  real*8                                   :: multfactor
  real*8                                   :: rho_tiny
  integer                                  :: interpolate_order
integer :: ii
  integer                                  :: adjimin, adjjmin, adjkmin
  integer                                  :: adjimax, adjjmax, adjkmax
  CCTK_REAL :: reduction_value
  CCTK_INT :: red_tmp
  integer :: header_flag,handle,dummy
  integer :: index,num_cols
  integer :: ierr,myproc_rank
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

!  write(6,*)'bondi diagnostics!'

  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)
  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax2,1,index)

  call compute_adj_minmax(ext,X,Y,Z,Symmetry, &
       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax, &
       multfactor,cctk_nghostzones,xmax2)
 
  dtrhob=(rho_b-dtrhob)/dT

!  call setup_global_coord_arrays(CCTK_PASS_FTOF)

  !header_flag == 1 -> write file headers
  if(CCTK_TIME == 0.D0) then
     header_flag = 1
  else
     header_flag = 0
  end if

  do ii=1,nsurfrho

  !--------------------------------------------------------------------------------!
     
     if(nsurfrho.eq.1) then
        filename = 'bondi.mon'
     else
        if(ii.le.9) then
           write(filename,31)ii
        else
           write(filename,32)ii
        endif
31      FORMAT('bondi.mon.',I1)
32      FORMAT('bondi.mon.',I2)
     endif

     rhosurf=rhosurfvec(ii)
     u0sch=u0vec(ii)

     num_cols = 6
     
     data_headers(1) = "Time"
     export_data(1) = CCTK_TIME
     
     call rhoflux(CCTK_PASS_FTOF)

     data_headers(2) = "Mass flux"
     export_data(2) = out_surf_int_sum
     data_headers(3) = "isosurface area"
     export_data(3) = F_M0
     data_headers(4) = "z axis"
     export_data(4) = out_surf_px
     data_headers(5) = "x axis"
     export_data(5) = out_surf_py
     data_headers(6) = "xz"
     export_data(6) = out_surf_pz
!     write(6,*)'mass flux:',out_surf_int_sum
!     write(6,*)'isosurface area:',F_M0
     
     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
     
     !--------------------------------------------------------------------------------!
     
  enddo

  if(ahf_active == .true.) then
     filename = 'bondi.hon'
     
     num_cols = 2
     
     data_headers(1) = "Time"
     export_data(1) = CCTK_TIME
     
     call ahfluxes(CCTK_PASS_FTOF)
     data_headers(2) = "Horizon flux"
     export_data(2) = F_M0
!     write(6,*)'Horizon flux:',F_M0
     
     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
  endif
  
  dtrhob = rho_b

end subroutine bondi_diag_integrals
