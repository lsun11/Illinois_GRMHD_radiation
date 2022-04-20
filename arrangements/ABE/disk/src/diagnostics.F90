!------------------------------------------------
! Diagnostics driver routine for disk evolutions
!------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine diagnostic_integrals_disk(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)    :: ext
  real*8, dimension(50)    :: export_data
  character, dimension(50) :: data_headers*20
  character                :: filename*50
  real*8                   :: dT,dX,dY,dZ,xmax
  real*8                   :: Mass0,Jtot,temp
  real*8                   :: dIntegral,multfactor,dIntegral_Mon
  real*8                   :: dIntegral_M, dIntegral_J
  real*8                   :: Eem
  real*8                   :: max_b2p, max_b2rho, beta_m1
  real*8                   :: Bx_max, By_max, Bz_max
  real*8                   :: Bx_min, By_min, Bz_min
  real*8                   :: Monopole,Br_max,Bphi_max,Btheta_max
  real*8                   :: Br_min, Bphi_min, Btheta_min
  real*8                   :: F_M0_l, F_E_fluid_l, F_E_em_l
  real*8                   :: F_J_fluid_l, F_J_em_l
  real*8                   :: out_M0_flux1, out_fluid_energy_flux1
  real*8                   :: out_EM_energy_flux1, out_fluid_J_flux1
  real*8                   :: out_EM_J_flux1, out_M0_flux2,out_fluid_energy_flux2
  real*8                   :: out_EM_energy_flux2, out_fluid_J_flux2, out_EM_J_flux2
  real*8                   :: out_M0_flux3, out_fluid_energy_flux3, out_EM_energy_flux3
  real*8                   :: out_fluid_J_flux3, out_EM_J_flux3
  integer                  :: adjimin, adjjmin, adjkmin
  integer                  :: adjimax, adjjmax, adjkmax
  real*8                   :: w_eh, m0
  CCTK_REAL :: reduction_value
  CCTK_INT :: red_tmp
  integer :: header_flag,handle,dummy
  integer :: index,num_cols
  integer :: ierr,myproc_rank
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)

  call compute_adj_minmax(ext,X,Y,Z,Symmetry, &
       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax, &
       multfactor,cctk_nghostzones,xmax)

  !header_flag == 1 -> write file headers
  if(CCTK_TIME == 0.D0) then
     header_flag = 1
  else
     header_flag = 0
  end if

  !--------------------------------------------------------------------------------!
  ! Calculations for the mon file.
  !--------------------------------------------------------------------------------!

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::rho_b")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,rho_b_max,1,index)

  call compute_b2p_b2rho(ext,sbt,sbx,sby,sbz,rho_b,P,lapm1,shiftx,shifty, &
       shiftz, phi, gxx,gxy,gxz,gyy,gyz,gzz, rho_br,rho_bl,Pr,Symmetry)

  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')

  if (Symmetry == AXISYM .and. ext(2)==3 ) then
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::Pr')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::rho_br')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::rho_bl')
  end if

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::rho_br")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,max_b2rho,1,index)

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::rho_bl")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,max_b2p,1,index)

  call  compute_absdivb_disk(ext,dX,dY,dZ,X,Z,Bx,By,Bz,phi,Pr,Sym_Bz,Symmetry,excision_radius)

  ! I want rest mass and J outside the horizon.
  ! This flagrantly assumes axisymmetry.
  ! store rho_star and x^2*\tilde{S}_y in temporary gfs
  rho_br = rho_star
  rho_bl = X*mhd_st_y

  ! mask the integrands
  rho_br = rho_br * mskf
  rho_bl = rho_bl * mskf

  call mass_0_integral(ext, dIntegral_M, &
       X, Y, Z, &
       rho_br, Symmetry, &
       adjimin, adjimax, &
       adjjmin, adjjmax, &
       adjkmin, adjkmax,CCTK_MyProc(cctkGH)) 
  dIntegral_M = multfactor*dIntegral_M
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral_M,Mass0,CCTK_VARIABLE_REAL)

  call mass_0_integral(ext, dIntegral_J, &
       X, Y, Z, &
       rho_bl, Symmetry, &
       adjimin, adjimax, &
       adjjmin, adjjmax, &
       adjkmin, adjkmax,CCTK_MyProc(cctkGH)) 
  dIntegral_J = multfactor*dIntegral_J
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral_J,Jtot,CCTK_VARIABLE_REAL)

  call mass_0_integral(ext, dIntegral_Mon, &
       X, Y, Z, &
       Pr, Symmetry, &
       adjimin, adjimax, &
       adjjmin, adjjmax, &
       adjkmin, adjkmax,CCTK_MyProc(cctkGH)) 
  dIntegral_Mon = multfactor*dIntegral_Mon
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral_Mon,Monopole,CCTK_VARIABLE_REAL)

  call em_energy_disk(ext,dIntegral, &
       X,Y,Z,phi,gxx,gxy,gxz, &
       gyy,gyz,gzz, &
       Ex,Ey,Ez,Bx,By,Bz, &
       Symmetry, mskf, &
       adjimin,adjimax, &
       adjjmin,adjjmax, &
       adjkmin,adjkmax)
  dIntegral = multfactor*dIntegral
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,Eem,CCTK_VARIABLE_REAL)

  ! convert Bi to physical coordinates and store in temporary GFs
!!$     if (fisheye_enable == 1) then
!!$        call bi_phys(ext,X,Y,Z,PhysicalRadius,RadiusDerivative,Bx,By,Bz,vxr,vyr,vzr,Symmetry)
!!$     else
  vxr = Bx
  vyr = By
  vzr = Bz
  !     end if

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::vxr")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Bx_max,1,index)

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::vyr")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,By_max,1,index)

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::vzr")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Bz_max,1,index)

  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"mhd_evolve::vxr")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Bx_min,1,index)

  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"mhd_evolve::vyr")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,By_min,1,index)

  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"mhd_evolve::vzr")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Bz_min,1,index)

  if (abs(Bx_min) .gt. abs(Bx_max)) Bx_max = Bx_min
  if (abs(By_min) .gt. abs(By_max)) By_max = By_min
  if (abs(Bz_min) .gt. abs(Bz_max)) Bz_max = Bz_min

  call fluxes_ana(CCTK_PASS_FTOF)

  filename = 'disk.mon'
  export_data = 0.D0

  num_cols = 16

  data_headers(1) = 'Time'
  export_data(1) = CCTK_TIME

  data_headers(2) = 'Mass0'
  export_data(2) = Mass0

  data_headers(3) = 'Jtot'
  export_data(3) = Jtot

  data_headers(4) = 'rho_b_max'
  export_data(4) = rho_b_max

  data_headers(5) = 'Eem'
  export_data(5) = Eem

  data_headers(6) = 'max_b2rho'
  export_data(6) = max_b2rho

  data_headers(7) = 'max_b2p'
  export_data(7) = max_b2p

  data_headers(8) = 'Bx_max'
  export_data(8) = Bx_max

  data_headers(9) = 'By_max'
  export_data(9) = By_max

  data_headers(10) = 'Bz_max'
  export_data(10) = Bz_max

  data_headers(11) = 'F_M0'
  export_data(11) = F_M0

  data_headers(12) = 'F_E_fluid'
  export_data(12) = F_E_fluid

  data_headers(13) = 'F_E_em'
  export_data(13) = F_E_em

  data_headers(14) = 'F_J_fluid'
  export_data(14) = F_J_fluid

  data_headers(15) = 'F_J_em'
  export_data(15) = F_J_em

  data_headers(16) = 'Monopole'
  export_data(16) = Monopole

  if(CCTK_MyProc(CCTKGH) .eq. 0.D0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !--------------------------------------------------------------------------------!
  ! Calculations for the fon file
  !--------------------------------------------------------------------------------!

  filename = 'disk.fon'

  ! Calculate the fluxes through the various shells.
  surf_radius = r_out_flux1
  call fluxes(CCTK_PASS_FTOF)
  out_M0_flux1           = F_M0
  out_fluid_energy_flux1 = F_E_fluid
  out_EM_energy_flux1    = F_E_em
  out_fluid_J_flux1      = F_J_fluid
  out_EM_J_flux1         = F_J_em

  surf_radius = r_out_flux2
  call fluxes(CCTK_PASS_FTOF)
  out_M0_flux2           = F_M0
  out_fluid_energy_flux2 = F_E_fluid
  out_EM_energy_flux2    = F_E_em
  out_fluid_J_flux2      = F_J_fluid
  out_EM_J_flux2         = F_J_em

  surf_radius = r_out_flux3
  call fluxes(CCTK_PASS_FTOF)
  out_M0_flux3           = F_M0
  out_fluid_energy_flux3 = F_E_fluid
  out_EM_energy_flux3    = F_E_em
  out_fluid_J_flux3      = F_J_fluid
  out_EM_J_flux3         = F_J_em

  num_cols = 16
  export_data = 0.D0

  data_headers(1) = 'Time'
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

  if(CCTK_MyProc(CCTKGH) .eq. 0.D0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

end subroutine diagnostic_integrals_disk
