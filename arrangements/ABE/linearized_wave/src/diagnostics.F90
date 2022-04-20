#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------
! Linearized wave diagnostic output
!-----------------------------------
subroutine lw_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,ext_noghost
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20
  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,xmax,ymax,zmax,Zmin
  real*8                                   :: RestMass,ADMMass,Jtot,temp,resid,resid_norm
  real*8                                   :: Mass_GW,hplus,hcross
  real*8                                   :: hplus_anal,hcross_anal,E_GW_anal
  real*8                                   :: psi4_22_ana,psi4_21_ana,psi4_20_ana
  real*8,dimension(2)                      :: gw_amplitudes
  real*8,dimension(1,3)                    :: pointcoords
  real*8,dimension(1)                      :: interp_output
  real*8                                   :: resid_x,resid_y,resid_z,resid_mom_norm,gxx_a
  real*8                                   :: dIntegral
  real*8                                   :: gxres,gyres,gzres
  real*8                                   :: SymmFactor
  real*8                                   :: rho_tiny,psi4r,psi4i
  real*8                                   :: radius_gworig,theta_gworig,phi_gworig
  integer                                  :: interpolate_order
  integer                                  :: ii,yoriginindex

  CCTK_REAL :: reduction_value
  CCTK_INT :: red_tmp
  integer :: header_flag,handle,dummy
  integer :: vindex,num_cols
  integer :: ierr,myproc_rank
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  ext = cctk_lsh
  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(MOD(cctk_iteration,out_every)==0) then
     call set_SymmFactor(Symmetry,SymmFactor)

     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"grid::X")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,vindex)
     call CCTK_VarIndex(vindex,"grid::Y")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymax,1,vindex)
     call CCTK_VarIndex(vindex,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax,1,vindex)

     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(vindex,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Zmin,1,vindex)

     !  call setup_global_coord_arrays(CCTK_PASS_FTOF)

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1
     else
        header_flag = 0
     end if

     !--------------------------------------------------------------------------------!
     filename = 'lw.wave'
     num_cols = 11

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = "rad_flux"
     Mass_GW = 0.0D0

     ! Initialize several variables at t=0 
     if (cctk_iteration==0) then
        momentsr_old = 0.d0
        momentsi_old = 0.d0
        odd_momentsr_old = 0.d0
        odd_momentsi_old = 0.d0
        int_momentsr = 0.d0
        int_momentsi = 0.d0
        momentsr = 0.D0
        momentsi = 0.D0
     end if

     gw_amplitudes = 0.D0
     hplus = 0.d0
     hcross = 0.d0
     hplus_anal = 0.d0
     hcross_anal = 0.d0
     E_GW_anal = 0.d0

     !radius_GW=5.05049502524257576
     !theta_gw= 1.53121499085584
     !phi_gw =  9.90066658798857D-003

     !  radius_GW = sqrt(X(10,1,1)**2 + dY*dY*0.25D0 + dZ*dZ*0.25D0)
     !  theta_gw = atan(sqrt(X(10,1,1)**2 + dY*dY*0.25D0)/dZ*0.5D0)
     !  phi_gw = atan(dY*0.5D0/X(10,1,1))  
     !  write(*,*) "ll",radius_GW,theta_GW,phi_GW

     call gw_wave_flux(cctkGH,radius_GW,Mass_GW, &
          dT,momentsr,momentsi, &
          momentsr_old,momentsi_old,odd_momentsr,odd_momentsi, &
          odd_momentsr_old,odd_momentsi_old,int_momentsr,int_momentsi, &
          hplus,hcross,theta_GW,phi_GW, &
          Zmin,dZ)

     call gw_anal(CCTK_TIME+time_shift,radius_GW_phys,theta_GW,phi_GW, &
          hplus_anal,hcross_anal,E_GW_anal,psi4r,psi4i, & 
          psi4_22_ana, psi4_21_ana,psi4_20_ana, & 
          amplitude,width,mode)

     data_headers(2) = "hplus"
     data_headers(3) = "hcross"
     export_data(2) = hplus
     export_data(3) = hcross
     !hplus_anal and hcross_anal are multiplied by radius_gw_phys in the above subroutines
     data_headers(4) = "hp_ana"
     data_headers(5) = "hc_ana" 
     export_data(4) = hplus_anal/radius_GW_phys
     export_data(5) = hcross_anal/radius_GW_phys

     data_headers(6) = "E_GW_ana"
     export_data(6) = E_GW_anal

     data_headers(7) = "Re(psi4)" 
     export_data(7) = psi4r

     data_headers(8) = "Im(psi4)"
     export_data(8) = psi4i

!     call leq2_meq2_mode_average_lw(CCTK_PASS_FTOF)
!     data_headers(9) = "Re(psi4)_22" 
!     export_data(9) = Psi4resumlw
     
     data_headers(9) = "psi4_22_ana"
     export_data(9) = psi4_22_ana

     data_headers(10) = "psi4_21_ana"
     export_data(10) = psi4_21_ana

     data_headers(11) = "psi4_20_ana"
     export_data(11) = psi4_20_ana

!     data_headers(10) = "Im(psi4)_22"
!     export_data(10) = Psi4imsumlw

!     data_headers(11) = "momentsr(1)"
!     export_data(11) = momentsr(1)
!     data_headers(12) = "momentsr(2)"
!     export_data(12) = momentsr(2)

!     data_headers(13) = "momentsi(1)"
!     export_data(13) = momentsi(1)
!     data_headers(14) = "momentsi(2)"
!     export_data(14) = momentsi(2)

!     data_headers(15) = "odd_momentsr(1)"
!     export_data(15) = odd_momentsr(1)
!     data_headers(16) = "odd_momentsi(1)"
!     export_data(16) = odd_momentsi(1)

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !--------------------------------------------------------------------------------!

     !Following is here so we can easily compare against analytic solution.
     filename = 'lw.mon'
     num_cols = 4

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     pointcoords(1,1) = Xmax/4.D0
     pointcoords(1,2) = Ymax/4.D0
     pointcoords(1,3) = 0.D0

     data_headers(2) = "gxx_analytic"
     !     call lin_wave_analytic_gxx(ext, X, Y, Z, cctk_time+time_shift, amplitude, width, PsiTau, &
     !          PhysicalRadius,RadiusDerivative, mode)
     call CCTK_VarIndex(vindex,"bssn::PsiTau")
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,interp_output)

     export_data(2) = interp_output(1)
     gxx_a = interp_output(1)

     data_headers(3) = "gxx_numeric"
     call CCTK_VarIndex(vindex,"bssn::gxx")
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,interp_output)
     export_data(3) = interp_output(1)

     data_headers(4) = "gxx_a-gxx_n"
     export_data(4) = gxx_a - interp_output(1)

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !     PsiTau = abs(PsiTau - gxx)

     !  call lin_wave_analytic(ext, X, Y, Z, cctk_time+time_shift, amplitude, width, &
     !       gxx_analy,gxy_analy,gxz_analy,gyy_analy,gyz_analy,gzz_analy,Kxx_analy,Kxy_analy,Kxz_analy,Kyy_analy,Kyz_analy,Kzz_analy, &
     !       PhysicalRadius,RadiusDerivative, mode)

  end if

end subroutine lw_diagnostics
