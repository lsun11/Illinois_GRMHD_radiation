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
subroutine OS_rad_diag_integrals(CCTK_ARGUMENTS)

  implicit none
  interface
     subroutine OS_rad_analytic_new(Fcoeff,z,R_a,tau,RoM,M,Eo4PiB,tau_a,tau_s,eta,Q,rho_b,Ec,Fc,E,F)
       use nr
       use nrtype
       implicit none

       real(DP), intent(in) :: R_a,tau,RoM,M,Eo4PiB,tau_a,tau_s
       real(DP), intent(out) :: E,F,rho_b,z,Fcoeff
       real(DP) :: eta,tau_t,sinchi0,chi0,Q,Ec,Fc,x1,x2,xacc,term
       integer :: n
       interface
          FUNCTION func(eta,tau,RoM,M)
            USE nrtype
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: eta,tau,RoM,M
            REAL(DP) :: func
          END FUNCTION func
       end interface
     end subroutine OS_rad_analytic_new
  end interface

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,ext_noghost,global_ext
  !  real*8, dimension(cctk_gsh(3)-1,3)                     :: pointcoords
  !  real*8, dimension(cctk_gsh(3)-1)                     :: outputdata
  real*8, dimension(20)                     :: export_data
  character, dimension(20)                     :: data_headers*20
  character                                :: filename*25,file_suffix*10
  real*8                                   :: dT,dX,dY,dZ,xmax
  real*8                                   :: RestMass,ADMMass,Jtot,temp,resid,resid_norm
  real*8                                   :: resid_x,resid_y,resid_z,resid_mom_norm
  real*8                                   :: dIntegral
  real*8                                   :: gxres,gyres,gzres
  real*8                                   :: multfactor
  real*8                                   :: rho_tiny
  integer                                  :: interpolate_order

  !Needed for GW extraction:
  real*8                                   :: Mass_GW,E_GW,hplus,hcross,J_GW
  real*8                                   :: psi4r,psi4i
  real*8,dimension(2)                      :: gw_amplitudes
  !End GW extraction variables  

  integer                                  :: adjimin, adjjmin, adjkmin
  integer                                  :: adjimax, adjjmax, adjkmax
  CCTK_REAL :: reduction_value
  CCTK_INT :: red_tmp
  integer :: header_flag,handle,dummy
  integer :: index,num_cols
  integer :: ierr,myproc_rank
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8                                     :: PI
  real*8                                     :: rho_b_initial
  integer :: ii,jj
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  PI = acos(-1.d0)

  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)
  global_ext = cctk_gsh

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

  !  call setup_global_coord_arrays(CCTK_PASS_FTOF)

  !header_flag == 1 -> write file headers
  if(CCTK_TIME == 0.D0) then
     header_flag = 1
  else
     header_flag = 0
  end if

  !--------------------------------------------------------------------------------!

  filename = 'OS_rad.mon'
  num_cols = 7
  export_data = 0.D0

  data_headers(1) = 'Time'
  export_data(1) = CCTK_TIME

  data_headers(2) = 'Phicent'
  call interp_center(cctkGH,cctk_nghostzones,phi,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(2))

  data_headers(3) = 'Kcent'
  call interp_center(cctkGH,cctk_nghostzones,trK,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(3)) 

  data_headers(4) = 'gxxcent-1'
  call interp_center(cctkGH,cctk_nghostzones,gxx,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(4)) 
  export_data(4) = export_data(4) - 1

  data_headers(5) = 'rhocent'
  call interp_center(cctkGH,cctk_nghostzones,rho_b,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(5)) 
  write(*,*) "rho center = ",export_data(5)

  data_headers(6) = 'lapsecent'
  call interp_center(cctkGH,cctk_nghostzones,lapm1,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(6)) 
  export_data(6) = export_data(6) + 1

  data_headers(7) = 'Omega_c'
  Pr = (vy*X - vx*Y) / (X*X + Y*Y)
  call interp_center(cctkGH,cctk_nghostzones,Pr,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(7))

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !--------------------------------------------------------------------------------!

  filename = 'OS_rad.don'

  num_cols = 10
  export_data = 0.D0

  data_headers(1) = 'Time'
  export_data(1) = CCTK_TIME

  data_headers(2) = 'M0'
  call mass_0_integral(ext, dIntegral, &
       X, Y, Z, &
       rho_star, Symmetry, &
       adjimin, adjimax, &
       adjjmin, adjjmax, &
       adjkmin, adjkmax,CCTK_MyProc(cctkGH)) 
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(2),CCTK_VARIABLE_REAL)

  data_headers(3) = 'Eint'
  call minternal_hybrid(ext,dIntegral, X,Y,Z, &
       rho_b_atm,rho_star,rho_b,P, &
       neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
       Symmetry,adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax); 
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(3),CCTK_VARIABLE_REAL)

  data_headers(4) = 'T'
  call tint(ext,dIntegral, &
       X,Y,Z,phi, &
       vx,vy,st_x,st_y, &
       Symmetry, &
       adjimin,adjimax, &
       adjjmin,adjjmax, &
       adjkmin,adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(4),CCTK_VARIABLE_REAL)

  data_headers(5) = 'M_ADM'
  call mass_integral(ext,dIntegral,X,Y,Z, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Axx,Axy,Axz,Ayy,Ayz,Azz,trK, &
       trRtilde,rho,phi, &
       Gammaxxx,Gammaxxy,Gammaxxz, &
       Gammaxyy,Gammaxyz,Gammaxzz, &
       Gammayxx,Gammayxy,Gammayxz, &
       Gammayyy,Gammayyz,Gammayzz, &
       Gammazxx,Gammazxy,Gammazxz, &
       Gammazyy,Gammazyz,Gammazzz, &
       Symmetry, &
       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(5),CCTK_VARIABLE_REAL)

  data_headers(6) = 'M_Komar'
  call komar_mass_integral(ext,dIntegral,X,Y,Z,phi,lapm1,rho,S, &
       shiftx,shifty,shiftz, Sx,Sy,Sz, &
       Symmetry, adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(6),CCTK_VARIABLE_REAL)

  data_headers(7) = 'Eem'
  call em_energy(ext,dIntegral, &
       X,Y,Z,phi,gxx,gxy,gxz, &
       gyy,gyz,gzz, &
       Ex,Ey,Ez,Bx,By,Bz, &
       Symmetry, &
       adjimin,adjimax, &
       adjjmin,adjjmax, &
       adjkmin,adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(7),CCTK_VARIABLE_REAL)

  data_headers(8) = 'Eem2'
  call em_energy2(ext,dIntegral, &
       X,Y,Z,lapm1, &
       vx,vy,vz, &
       shiftx,shifty,shiftz, &
       phi,gxx,gxy,gxz, &
       gyy,gyz,gzz, &
       sbt,sbx,sby,sbz, &
       Symmetry, &
       adjimin,adjimax, &
       adjjmin,adjjmax, &
       adjkmin,adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(8),CCTK_VARIABLE_REAL)

  data_headers(9) = 'int B_phi'
  call b_phi(ext,dIntegral, rho_star, &
       X,Y,Z,phi, &
       Bx,By,Bz, &
       gxy,gyy,gyz, &
       Symmetry, &
       adjimin,adjimax, &
       adjjmin,adjjmax, &
       adjkmin,adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(9),CCTK_VARIABLE_REAL)

  data_headers(10) = 'By_max'
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::By")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(10),1,index)

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !I don't know what's wrong here.  I'll try and fix it later

  !--------------------------------------------------------------------------------!
!!$  filename = 'OS_rad.jon'
!!$  num_cols = 11
!!$
!!$  data_headers(1) = "Time"
!!$  export_data(1) = CCTK_TIME
!!$
!!$  data_headers(2) = "AngMom"
!!$  !we use Pr as a temporary variable below:
!!$  call jdens_cpp(cctkGH,cctk_lsh,cctk_nghostzones, Symmetry, &
!!$       Pr,  &
!!$       dx,dy,dz, &
!!$       X,Y,Z,  &
!!$       phi, Sx, Sy, trK, & 
!!$       Axx, Axy, Axz, Ayy, Ayz, Azz,  &
!!$       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
!!$
!!$!  call jdens(ext,Pr, &
!!$!       X,Y,Z,phi,Sx,Sy,trK, &
!!$!       Axx,Axy,Axz,Ayy,Ayz,Azz, &
!!$!       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
!!$!       Symmetry)
!!$
!!$  call mass_0_integral(ext,dIntegral,X,Y,Z,Pr,Symmetry, &
!!$       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
!!$  dIntegral = dIntegral*multfactor
!!$  call CCTK_ReductionHandle(handle,"sum")
!!$  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(2),CCTK_VARIABLE_REAL)
!!$
!!$
!!$  data_headers(3) = "J_out"
!!$  call J_surf_integral(CCTK_PASS_FTOF)
!!$  export_data(3) = out_surf_int_sum
!!$
!!$  data_headers(4) = "Mass_sur"
!!$  call surf_Mass_fisheye(CCTK_PASS_FTOF)
!!$  export_data(4) = out_surf_int_sum
!!$
!!$  data_headers(5) = "Komar Mass"
!!$  call Komar_Mass(CCTK_PASS_FTOF)
!!$  export_data(5) = out_surf_int_sum
!!$
!!$  if(excision_enable==0) then
!!$     data_headers(6) = "out_M0_flux"
!!$     data_headers(7) = "out_fluid_Eng_flux"
!!$     data_headers(8) = "out_EM_Energy_flux"
!!$     data_headers(9) = "out_fluid_J_flux"
!!$     data_headers(10) = "out_EM_J_flux"
!!$     call fluxes(CCTK_PASS_FTOF)
!!$     export_data(6) = F_M0
!!$     export_data(7) = F_E_fluid
!!$     export_data(8) = F_E_em
!!$     export_data(9) = F_J_fluid
!!$     export_data(10) = F_J_em  
!!$  else
!!$     data_headers(6) = "out_M0_flux"
!!$     data_headers(7) = "out_fluid_Eng_flux"
!!$     data_headers(8) = "out_EM_Energy_flux"
!!$     data_headers(9) = "out_fluid_J_flux"
!!$     data_headers(10) = "out_EM_J_flux"
!!$     call ahfluxes(CCTK_PASS_FTOF)
!!$     export_data(6) = F_M0
!!$     export_data(7) = F_E_fluid
!!$     export_data(8) = F_E_em
!!$     export_data(9) = F_J_fluid
!!$     export_data(10) = F_J_em  
!!$  end if
!!$
!!$
!!$
!!$  data_headers(11) = "rho_isoflux"
!!$  out_surf_int_sum = 0.D0
!!$!  call rhoflux(CCTK_PASS_FTOF)
!!$  export_data(11) = out_surf_int_sum
!!$
!!$  if(CCTK_MyProc(CCTKGH) .eq. 0) then
!!$     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
!!$  end if
  !--------------------------------------------------------------------------------!

  filename = 'OS_rad.con'
  num_cols = 8

  data_headers(1) = "Time"
  export_data(1) = CCTK_TIME


  data_headers(2) = "Haml Res"
  !************ HAMILTONIAN CONSTRAINT ************!
  PsiTau = PsiRes
  ! Next the L2 norm.  Note I don't need to Sync since I don't care
  ! about the ghost zones.
  dIntegral = 0.D0
  call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
  resid = sqrt(resid)

  dIntegral = 0.D0
  call norm2_dv(ext,PsiNorm,X,Y,Z,dIntegral,adjimin, &
       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry);
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
  resid_norm = sqrt(resid_norm)
  export_data(2) = sqrt(multfactor)*resid/resid_norm
  !************ END HAMILTONIAN CONSTRAINT ************!

  data_headers(3) = "Momx Res"
  data_headers(4) = "Momy Res"
  data_headers(5) = "Momz Res"
  !************ MOMENTUM CONSTRAINT ************!
  call norm2_dv(ext,MRsx,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_x,CCTK_VARIABLE_REAL)
  resid_x = sqrt(resid_x)

  call norm2_dv(ext,MRsy,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_y,CCTK_VARIABLE_REAL)
  resid_y = sqrt(resid_y)

  dIntegral = dIntegral*multfactor
  call norm2_dv(ext,MRsz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_z,CCTK_VARIABLE_REAL)
  resid_z = sqrt(resid_z)

  call norm2_dv(ext,MNorm,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_mom_norm,CCTK_VARIABLE_REAL)
  resid_mom_norm = sqrt(resid_mom_norm)

  export_data(3)=sqrt(multfactor)*resid_x/resid_mom_norm
  export_data(4)=sqrt(multfactor)*resid_y/resid_mom_norm
  export_data(5)=sqrt(multfactor)*resid_z/resid_mom_norm

  !************ END MOMENTUM CONSTRAINT ************!

  data_headers(6) = "Gamx Res"
  data_headers(7) = "Gamy Res"
  data_headers(8) = "Gamz Res"
  !************ GAMMA CONSTRAINT ************!
  call gamcheck(cctkGH,cctk_lsh, cctk_nghostzones, Symmetry, &
       dx,  dy,  dz,  &
       gconx, gcony, gconz, &
       Gammax, Gammay, Gammaz,  &
       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
  call norm2_dv(ext,gconx,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gxres,CCTK_VARIABLE_REAL)
  dIntegral = dIntegral*multfactor
  export_data(6) = sqrt(multfactor)*sqrt(gxres)
  call norm2_dv(ext,gcony,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gyres,CCTK_VARIABLE_REAL)
  export_data(7) = sqrt(multfactor)*sqrt(gyres)
  call norm2_dv(ext,gconz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gzres,CCTK_VARIABLE_REAL)
  export_data(8) = sqrt(multfactor)*sqrt(gzres)
  !************ END GAMMA CONSTRAINT ************!

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if
  !--------------------------------------------------------------------------------!
  filename = 'OS_rad.ron'
  num_cols = 6

  data_headers(1) = "Time"
  export_data(1) = CCTK_TIME

  !# time         Gamx_max         rho_b_max        Ang_Mom          lapse_min        Kpoly_avg        Rc               Rp               max(b^2/rho0)    max(b^2/P)       
  data_headers(2) = "Gammaxmax"
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"bssn::Gammax")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(2),1,index)

  data_headers(3) = "rho_bmax"
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::rho_b")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(3),1,index)

  data_headers(4) = "AngMom"
  !we use Pr as a temporary variable below:
  call jdens_cpp(cctkGH,cctk_lsh,cctk_nghostzones, Symmetry, &
       Pr,  &
       dx,dy,dz, &
       X,Y,Z,  &
       phi, Sx, Sy, trK, & 
       Axx, Axy, Axz, Ayy, Ayz, Azz,  &
       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)


  !  call jdens(ext,Pr, &
  !       X,Y,Z,phi,Sx,Sy,trK, &
  !       Axx,Axy,Axz,Ayy,Ayz,Azz, &
  !       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
  !       Symmetry)

  call mass_0_integral(ext,dIntegral,X,Y,Z,Pr,Symmetry, &
       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(4),CCTK_VARIABLE_REAL)

  data_headers(5) = "lapse_min"
  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"lapse::lapm1")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(5),1,index)
  export_data(5) = export_data(5) + 1.D0

  data_headers(6) = "Rp"
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::rho_b")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,rho_tiny,1,index)
  rho_tiny = rho_tiny * 1.D-3
  call find_equatorial_radius(ext,cctkgh,X,Y,Z,dX,dY,dZ,rho_tiny,rho_star,phi,gxx,gxy,gyy,Symmetry,export_data(6),CCTK_MyProc(CCTKGH))

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !--------------------------------------------------------------------------------!
  filename = 'OS_rad.xon'
  num_cols = 5

  data_headers(1) = "Time"
  export_data(1) = CCTK_TIME

  data_headers(2) = "Xcen"
  data_headers(3) = "Ycen"
  data_headers(4) = "Zcen"
  !First (re)compute the Rest Mass (recompute if necessary; this fileio module should NOT depend on others).
  call mass_0_integral(ext, dIntegral, &
       X, Y, Z, &
       rho_star, Symmetry, &
       adjimin, adjimax, &
       adjjmin, adjjmax, &
       adjkmin, adjkmax,CCTK_MyProc(cctkGH)) 
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,RestMass,CCTK_VARIABLE_REAL)

  if (Symmetry==NO_SYMM) then
     call wavg(ext,dIntegral,X,Y,Z,X,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(2),CCTK_VARIABLE_REAL)

     call wavg(ext,dIntegral,X,Y,Z,Y,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(3),CCTK_VARIABLE_REAL)

     call wavg(ext,dIntegral,X,Y,Z,Z,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(4),CCTK_VARIABLE_REAL)
  else if (Symmetry==EQUATORIAL) then
     call wavg(ext,dIntegral,X,Y,Z,X,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(2),CCTK_VARIABLE_REAL)

     call wavg(ext,dIntegral,X,Y,Z,Y,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(3),CCTK_VARIABLE_REAL)

     export_data(4) = 0.D0
  else if (Symmetry==AXISYM .and. Z(2,2,2) < 0.0) then
     export_data(2) = 0.D0
     export_data(3) = 0.D0

     call wavg(ext,dIntegral,X,Y,Z,Z,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(4),CCTK_VARIABLE_REAL)

  else 
     export_data(2) = 0.D0
     export_data(3) = 0.D0
     export_data(4) = 0.D0
  end if

  data_headers(5) = "Rcen"
  export_data(5) = sqrt(export_data(2)*export_data(2) + export_data(3)*export_data(3) + export_data(4)*export_data(4))

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !--------------------------------------------------------------------------------!
  if (mod(CCTK_iteration,out_every) .eq. 0) then

     rho_b_initial=1.d0/(4.d0/3.d0*PI*RoM**3)

     do ii=1,3
        do jj=1,narr
           pos(jj,ii) = coord(jj,ii+1)
        end do
     end do

     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,ext,global_ext,&
          dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,&
          pos,gxx,gxx_pcle) 
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,ext,global_ext,&
          dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,&
          pos,gxz,gxz_pcle) 
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,ext,global_ext,&
          dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,&
          pos,gzz,gzz_pcle) 
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,ext,global_ext,&
          dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,&
          pos,phi,phi_pcle) 
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,ext,global_ext,&
          dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,&
          pos,rho_b,rho_b_pcle)
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,ext,global_ext,&
          dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,&
          pos,P,P_pcle)
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,ext,global_ext,&
          dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,&
          pos,E_rad,E_rad_pcle)
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,ext,global_ext,&
          dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,&
          pos,F,F_pcle)


     areal_radius=exp(2.d0*phi_pcle)*sqrt(coord(:,4)*coord(:,4)*gxx_pcle-2.d0*coord(:,2)*coord(:,4)*gxz_pcle+coord(:,2)*coord(:,2)*gzz_pcle)

     do ii=1,narr
        call OS_rad_analytic_new(F_coeff(ii),z_anal(ii),areal_radius(ii),&
             coord(ii,1),RoM,1.d0,1.d0,opt_depth_a,opt_depth_s,eta(ii),Q(ii),&
             rho_b_anal(ii),Ec_anal(ii),Fc_anal(ii),E_anal(ii),F_anal(ii))

        Ec_pcle(ii) = E_rad_pcle(ii)*Q(ii)**4
        Fc_pcle(ii) = F_pcle(ii)*Q(ii)**4

        !E_rad_gradient(1)=F_coeff(2)*(E_rad_pcle(2)-E_rad_pcle(1))/(z_anal(2)-z_anal(1))
        !if (ii .gt. 1) E_rad_gradient(ii)=F_coeff(ii)*(E_rad_pcle(ii)-E_rad_pcle(ii-1))/(z_anal(ii)-z_anal(ii-1))
     end do

     filename = 'OS_rad_mid.bon'
 
     num_cols = 16

     data_headers(1) = "Coord time"                                             
     export_data(1)=CCTK_TIME
     data_headers(2) = "tau"                        
     export_data(2)=coord(narr/2,1)
     data_headers(3) = "tc/tcoll"        
     export_data(3)=eta(narr/2)/PI
     data_headers(4) = "R_a"
     export_data(4)=areal_radius(narr/2)
     data_headers(5) = "rho_b"   
     export_data(5)=rho_b_pcle(narr/2)
     data_headers(6) = "P"   
     export_data(6)=P_pcle(narr/2)
     data_headers(7) = "E_rad"   
     export_data(7)=E_rad_pcle(narr/2)
     data_headers(8) = "F"   
     export_data(8)=F_pcle(narr/2)
     data_headers(9) = "Ec"   
     export_data(9)=Ec_pcle(narr/2)
     data_headers(10) = "Fc"   
     export_data(10)=Fc_pcle(narr/2)
     data_headers(11) = "rho_b_anal"   
     export_data(11)=rho_b_anal(narr/2)*rho_b_initial
     data_headers(12) = "P_anal"   
     export_data(12)=rho_b_anal(narr/2)*E_anal(narr/2)**0.25d0*PoRho*rho_b_initial
     data_headers(13) = "E_rad_anal"   
     export_data(13)=E_anal(narr/2)*rho_b_initial*PoRho/Po4PiB
     data_headers(14) = "F_anal"   
     export_data(14)=F_anal(narr/2)*rho_b_initial*PoRho/Po4PiB
     data_headers(15) = "Ec_anal"   
     export_data(15)=Ec_anal(narr/2)*rho_b_initial*PoRho/Po4PiB
     data_headers(16) = "Fc_anal"   
     export_data(16)=Fc_anal(narr/2)*rho_b_initial*PoRho/Po4PiB

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
     !--------------------------------------------------------------------------------!

     filename = 'OS_rad_edge.bon'

     num_cols = 16

     data_headers(1) = "Coord time"                                             
     export_data(1)=CCTK_TIME
     data_headers(2) = "tau"                        
     export_data(2)=coord(narr,1)
     data_headers(3) = "tc/tcoll"        
     export_data(3)=eta(narr)/PI
     data_headers(4) = "R_a"
     export_data(4)=areal_radius(narr)
     data_headers(5) = "rho_b"   
     export_data(5)=rho_b_pcle(narr)
     data_headers(6) = "P"   
     export_data(6)=P_pcle(narr)
     data_headers(7) = "E_rad"   
     export_data(7)=E_rad_pcle(narr)
     data_headers(8) = "F"   
     export_data(8)=F_pcle(narr)
     data_headers(9) = "Ec"   
     export_data(9)=Ec_pcle(narr)
     data_headers(10) = "Fc"   
     export_data(10)=Fc_pcle(narr)
     data_headers(11) = "rho_b_anal"   
     export_data(11)=rho_b_anal(narr)*rho_b_initial
     data_headers(12) = "P_anal"   
     export_data(12)=rho_b_anal(narr)*E_anal(narr)**0.25d0*PoRho*rho_b_initial
     data_headers(13) = "E_rad_anal"   
     export_data(13)=E_anal(narr)*rho_b_initial*PoRho/Po4PiB
     data_headers(14) = "F_anal"   
     export_data(14)=F_anal(narr)*rho_b_initial*PoRho/Po4PiB
     data_headers(15) = "Ec_anal"   
     export_data(15)=Ec_anal(narr)*rho_b_initial*PoRho/Po4PiB
     data_headers(16) = "Fc_anal"   
     export_data(16)=Fc_anal(narr)*rho_b_initial*PoRho/Po4PiB

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
  endif
  !--------------------------------------------------------------------------------!

  !particle tracer stuff for analytic solutions

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     filename = 'analytic_profiles.xg'
     if (CCTK_TIME .eq. 0) then
        open(UNIT=15,FILE=filename,STATUS="REPLACE")
        write(15,*) ""
        write(15,*) '"analytic profiles'
     else
        open(UNIT=15,FILE=filename,STATUS="OLD",POSITION="APPEND")
     endif
     write(*,*) "iteration: ",CCTK_iteration
     if (mod(CCTK_iteration,out_every) .eq. 0) then
        write(15,*) ""
        write(15,*) ""
        write(15,*) '"Time =',CCTK_TIME

        write(15,'(t2,a6,t17,a1,t31,a5,t45,a1,t59,a1,t73,a1,t87,a2,t101,a2)') '"coord','z','rho_b','P','E','F','Ec','Fc'
        do ii=1,narr
           write(15,'(8e18.10)') coord(ii,2),z_anal(ii),rho_b_anal(ii)*rho_b_initial,&
                rho_b_anal(ii)*E_anal(ii)**0.25d0*PoRho*rho_b_initial,&
                E_anal(ii)*rho_b_initial*PoRho/Po4PiB,&
                F_anal(ii)*rho_b_initial*PoRho/Po4PiB,&
                Ec_anal(ii)*rho_b_initial*PoRho/Po4PiB,Fc_anal(ii)*rho_b_initial*PoRho/Po4PiB

        end do
     endif
     close(15)
  endif

  !--------------------------------------------------------------------------------!

  !numerical solutions at the positions of particles 

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     filename = 'numerical_profiles.xg'
     if (CCTK_TIME .eq. 0) then
        open(UNIT=17,FILE=filename,STATUS="REPLACE")
        write(17,*) ""
        write(17,*) '"numerical profiles'
     else
        open(UNIT=17,FILE=filename,STATUS="OLD",POSITION="APPEND")
     endif
     if (mod(CCTK_iteration,out_every) .eq. 0) then
        write(17,*) ""
        write(17,*) ""
        write(17,*) '"Time =',CCTK_TIME

        write(17,'(t2,a6,t17,a1,t31,a5,t45,a1,t59,a1,t73,a1,t87,a2,t101,a2)') '"coord','z','rho_b','P','E','F','Ec','Fc'
        do ii=1,narr
           write(17,'(8e18.10)') coord(ii,2),z_anal(ii),rho_b_pcle(ii),P_pcle(ii),E_rad_pcle(ii),F_pcle(ii),Ec_pcle(ii),Fc_pcle(ii)
        end do
     endif
     close(17)
  endif

  !----------------------------------------------------------------------------------!

  !convergence testing stuff
  write(*,*) "out_every: ",out_every
  write(*,*) "global_ext(1): ",global_ext(1)
  write(*,*) "coarse_res: ", coarse_res

  if (mod(CCTK_iteration,out_every) .eq. 0) then
     filename = 'convergence.xg'

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        if (CCTK_TIME .eq. 0) then
           open(UNIT=16,FILE=filename,STATUS="REPLACE")
           write(16,*) '"rho for convergence test'
           write(16,*) ""
           write(16,*) ""
           write(16,*) ""         
        else
           open(UNIT=16,FILE=filename,STATUS="OLD",POSITION="APPEND")
        endif
     endif
     write(*,*) "hello"
     write(16,*) ""
     write(16,*) ""
     write(16,*) '"Time =',CCTK_TIME
     write(*,*) "global_ext(1): ",global_ext
     do ii=1,coarse_res+2
        pos_int(ii,1)=(ii-1.5d0)*dX*(global_ext(1)-2)/coarse_res 
        pos_int(ii,2)=0.d0 
      ! write(*,*) "ii: ",ii
      ! write(*,*) "global_ext(3): ",global_ext(3)
      ! write(*,*) "coarse_res: ",coarse_res
      ! write(*,*) "dZ: ",dZ
       pos_int(ii,3)=0.5d0*dZ*(global_ext(3)-2)/coarse_res 
     end do

     rho_b_conv=0.d0
     P_conv=0.d0
     E_rad_conv=0.d0
     F_radx_conv=0.d0


     call interpolate_pointset(cctkGH,cctk_nghostzones,coarse_res,2,&
          ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
          Xglobal,Yglobal,Zglobal,pos_int,rho_b,rho_b_conv) 
     call interpolate_pointset(cctkGH,cctk_nghostzones,coarse_res,2,&
          ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
          Xglobal,Yglobal,Zglobal,pos_int,P,P_conv) 
     call interpolate_pointset(cctkGH,cctk_nghostzones,coarse_res,2,&
          ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
          Xglobal,Yglobal,Zglobal,pos_int,E_rad,E_rad_conv) 
     call interpolate_pointset(cctkGH,cctk_nghostzones,coarse_res,2,&
          ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
          Xglobal,Yglobal,Zglobal,pos_int,F_radx,F_radx_conv) 


     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        do ii=1,coarse_res+2
           write(16,'(5e18.10)') pos_int(ii,1),rho_b_conv(ii),P_conv(ii),E_rad_conv(ii),F_radx_conv(ii)
        end do
        close(16)
     endif
  endif

write(*,*) "END OF OS_rad_diag_integrals F_rad0(20,2,2): ",F_rad0(20,2,2)

end subroutine OS_rad_diag_integrals
