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
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3)                    :: ext,ext_noghost
!  real*8, dimension(cctk_gsh(3)-1,3)                     :: pointcoords
!  real*8, dimension(cctk_gsh(3)-1)                     :: outputdata
  real*8, dimension(20)                     :: export_data
  character, dimension(20)                     :: data_headers*20
  character                                :: filename*25,file_suffix*10
  real*8                                   :: dT,dX,dY,dZ,xmax
  real*8                                   :: RestMass,ADMMass,Jtot,temp,resid,resid_norm
  real*8                                   :: Mass_GW,E_GW,hpx_alt
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

  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)
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
  num_cols = 9
  export_data = 0.D0
  
  data_headers(1) = 'Time'
  export_data(1) = CCTK_TIME
  
  data_headers(3) = 'Phicent'
  call interp_center(cctkGH,phi,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(3))

!  call interp_center_old(ext,cctkGH,X,Y,Z,dX,dY,dZ,phi,Symmetry,export_data(2),CCTK_MyProc(CCTKGH))
!  call CCTK_VarIndex(index,"bssn::phi")
!  call interp_center(cctkGH,mag_interp_operator,index,export_data(2)) 

  data_headers(4) = 'Kcent'
  call interp_center(cctkGH,trK,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(4)) 
!  call interp_center(ext,cctkGH,X,Y,Z,dX,dY,dZ,trK,Symmetry,export_data(3),CCTK_MyProc(CCTKGH))
!  call CCTK_VarIndex(index,"bssn::trK")
!  call interp_center(cctkGH,mag_interp_operator,index,export_data(3))

  data_headers(5) = 'gxxcent-1'
  call interp_center(cctkGH,gxx,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(5)) 
!  call interp_center_old(ext,cctkGH,X,Y,Z,dX,dY,dZ,gxx,Symmetry,export_data(4),CCTK_MyProc(CCTKGH))
!  call CCTK_VarIndex(index,"bssn::gxx")
!  call interp_center(cctkGH,mag_interp_operator,index,export_data(4)) 
  export_data(5) = export_data(5) - 1

  data_headers(6) = 'rhocent'
  call interp_center(cctkGH,rho_b,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(6)) 
  !  call interp_center(ext,cctkGH,X,Y,Z,dX,dY,dZ,rho_b,Symmetry,export_data(5),CCTK_MyProc(CCTKGH))
  !  call CCTK_VarIndex(index,"mhd_evolve::rho_b")
  !  call interp_center(cctkGH,mag_interp_operator,index,export_data(5)) 

  data_headers(7) = 'lapsecent'
  call interp_center(cctkGH,lapm1,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(7)) 
  !  call interp_center(ext,cctkGH,X,Y,Z,dX,dY,dZ,lapm1,Symmetry,export_data(6),CCTK_MyProc(CCTKGH))
  !call CCTK_VarIndex(index,"lapse::lapm1")
  !  call interp_center(cctkGH,mag_interp_operator,index,export_data(6)) 
  export_data(7) = export_data(7) + 1
  
  data_headers(2) = "Proper time"
  write(*,*) "dT: ",dT
  write(*,*) "lapse: ",export_data(7)
  if (CCTK_TIME > 0) then
     proper_time = proper_time + dT*export_data(7)
  end if
  export_data(2)=proper_time
  
  data_headers(8) = 'Omega_c'
  Pr = (vy*X - vx*Y) / (X*X + Y*Y)
  call interp_center(cctkGH,Pr,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(8))
  !  call interp_center(ext,cctkGH,X,Y,Z,dX,dY,dZ,Pr,Symmetry,export_data(7),CCTK_MyProc(CCTKGH))
  !call CCTK_VarIndex(index,"mhd_evolve::Pr")
  !  call interp_center(cctkGH,mag_interp_operator,index,export_data(7)) 

!  write(*,*) "massblah0 = ",index
  data_headers(9) = "E_radcent"
  call interp_center(cctkGH,E_rad,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(9))
  
  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if
  
  !--------------------------------------------------------------------------------!
  
  filename = 'OS_rad.don'

  num_cols = 11
  export_data = 0.D0

  data_headers(1) = 'Time'
  export_data(1) = CCTK_TIME

  data_headers(2) = "Proper time"
  export_data(2)=proper_time   
  
  data_headers(3) = 'M0'
  call mass_0_integral(ext, dIntegral, &
       X, Y, Z, &
       rho_star, Symmetry, &
       adjimin, adjimax, &
       adjjmin, adjjmax, &
       adjkmin, adjkmax,CCTK_MyProc(cctkGH)) 
  !call CCTK_VarIndex(index,"magnetar::dIntegral")
!  call CCTK_ReductionHandle(handle,"sum")
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(3),CCTK_VARIABLE_REAL)
  !  call CCTK_Reduce (ierr,cctkGH, 0,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)

  data_headers(4) = 'Eint'
  call minternal_hybrid(ext,dIntegral, X,Y,Z, &
       rho_b_atm,rho_star,rho_b,P, &
       neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
       Symmetry,adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax); 
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(4),CCTK_VARIABLE_REAL)

  data_headers(5) = 'T'
  call tint(ext,dIntegral, &
       X,Y,Z,phi, &
       vx,vy,st_x,st_y, &
       Symmetry, &
       adjimin,adjimax, &
       adjjmin,adjjmax, &
       adjkmin,adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(5),CCTK_VARIABLE_REAL)

  data_headers(6) = 'M_ADM'
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
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(6),CCTK_VARIABLE_REAL)

  data_headers(7) = 'M_Komar'
  call komar_mass_integral(ext,dIntegral,X,Y,Z,phi,lapm1,rho,S, &
       shiftx,shifty,shiftz, Sx,Sy,Sz, &
       Symmetry, adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(7),CCTK_VARIABLE_REAL)

  data_headers(8) = 'Eem'
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
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(8),CCTK_VARIABLE_REAL)

  data_headers(9) = 'Eem2'
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
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(9),CCTK_VARIABLE_REAL)

  data_headers(10) = 'int B_phi'
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
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(10),CCTK_VARIABLE_REAL)

  data_headers(11) = 'By_max'
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::By")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(11),1,index)

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !--------------------------------------------------------------------------------!
  filename = 'OS_rad.jon'
  num_cols = 11

  data_headers(1) = "Time"
  export_data(1) = CCTK_TIME
  
  data_headers(2) = "Proper time"
  export_data(2)=proper_time  
  
  data_headers(3) = "AngMom"
  !we use Pr as a temporary variable below:
  call jdens(ext,Pr, &
       X,Y,Z,phi,Sx,Sy,trK, &
       Axx,Axy,Axz,Ayy,Ayz,Azz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Symmetry) !,imin,jmin,kmin,sNx,sNz, &
  !       PI_SYMM_K_gz,PI_SYMM_gupxx_gz,PI_SYMM_gupxy_gz,PI_SYMM_gupxz_gz, &
  !       PI_SYMM_gupyy_gz,PI_SYMM_gupyz_gz,PI_SYMM_gupzz_gz)
  call mass_0_integral(ext,dIntegral,X,Y,Z,Pr,Symmetry, &
       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(3),CCTK_VARIABLE_REAL)

 
!if(cctk_iteration.gt.10000) then
   data_headers(4) = "J_out"
  call J_surf_integral(CCTK_PASS_FTOF)
  export_data(4) = out_surf_int_sum
!end if
  
!if(cctk_iteration.gt.10000) then

  data_headers(5) = "Mass_sur"
  call surf_Mass_fisheye(CCTK_PASS_FTOF)
  export_data(5) = out_surf_int_sum

  data_headers(6) = "Komar Mass"
!  write(*,*) "Now computing",data_headers(5)
  call Komar_Mass(CCTK_PASS_FTOF)
!  write(*,*) "Finished computing",data_headers(5)
!end if

  export_data(6) = out_surf_int_sum
  data_headers(7) = "out_M0_flux"
  data_headers(8) = "out_fluid_Eng_flux"
  data_headers(9) = "out_EM_Energy_flux"
  data_headers(10) = "out_fluid_J_flux"
  data_headers(11) = "out_EM_J_flux"
!  write(*,*) "Now computing",data_headers(6)
  call fluxes(CCTK_PASS_FTOF)
!  write(*,*) "finished computing",data_headers(6)
  export_data(7) = F_M0
  export_data(8) = F_E_fluid
  export_data(9) = F_E_em
  export_data(10) = F_J_fluid
  export_data(11) = F_J_em  

!  data_headers(12) = "rho_isoflux"
 ! out_surf_int_sum = 0.D0
 ! call rhoflux(CCTK_PASS_FTOF)
 ! export_data(12) = out_surf_int_sum

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if
  !--------------------------------------------------------------------------------!

  filename = 'OS_rad.con'
  num_cols = 9
  
  data_headers(1) = "Time"
  
  export_data(1) = CCTK_TIME
  
  data_headers(2) = "Proper time"                                                                                                            
  export_data(2)=proper_time
  
  data_headers(3) = "Haml Res"
  !************ HAMILTONIAN CONSTRAINT ************!
  call pre_con(ext,phi,Psi,KK, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Axx,Axy,Axz,Ayy,Ayz,Azz,trK)
  call ham_residual(ext,X,Y,Z, &
       Psi,PsiRes,PsiTau,KK,trRtilde,rho, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Gammax,Gammay,Gammaz,Symmetry)
  !imin,jmin,kmin,sNx,sNz,PI_SYMM_Psi_gz)
  PsiTau = PsiRes
  !
  ! Next the L2 norm.  Note I don't need to Sync since I don't care
  ! about the ghost zones.
  !
  dIntegral = 0.D0
  call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
!  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
  resid = sqrt(resid)

  call ham_norm(ext,X,Y,Z, &
       Psi,PsiNorm,PsiTau,KK,trRtilde,rho, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Gammax,Gammay,Gammaz,Symmetry) !, &
  !imin,jmin,kmin,sNx,sNz,PI_SYMM_Psi_gz);
  dIntegral = 0.D0
  call norm2_dv(ext,PsiNorm,X,Y,Z,dIntegral,adjimin, &
       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry);
!  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
  resid_norm = sqrt(resid_norm)
  export_data(3) = sqrt(multfactor)*resid/resid_norm
  !write(*,*) "Resid = ",resid,resid_norm,resid/resid_norm
  !************ END HAMILTONIAN CONSTRAINT ************!

  data_headers(4) = "Momx Res"
  data_headers(5) = "Momy Res"
  data_headers(6) = "Momz Res"
  !************ MOMENTUM CONSTRAINT ************!
  call raise_a(ext,Axx, Axy, Axz, Ayy, Ayz, Azz, &
       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
       Rxx, Rxy, Rxz, Ryy, Ryz, Rzz)
  call mom_residual(ext,X,Y,Z,MRsx,MRsy,MRsz,MNorm, phi, &
       Axx, Axy, Axz, Ayy, Ayz, Azz, &
       Gammaxxx,Gammaxxy, &
       Gammaxxz,Gammaxyy,Gammaxyz,Gammaxzz, &
       Gammayxx,Gammayxy,Gammayxz,Gammayyy, &
       Gammayyz,Gammayzz,Gammazxx,Gammazxy, &
       Gammazxz,Gammazyy,Gammazyz,Gammazzz, &
       trK, Sx, Sy, Sz, &
       gupxx, gupxy, gupxz, gupyy, gupyz, & 
       gupzz, gxxx, gxyx, gxzx, gyyx, & 
       gyzx, gzzx, gxxy, gxyy, gxzy, & 
       gyyy, gyzy, gzzy, gxxz, gxyz, & 
       gxzz, gyyz, gyzz, gzzz, Rxx, & 
       Rxy, Rxz, Ryy, Ryz, Rzz, &
       Symmetry)
  call norm2_dv(ext,MRsx,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
!  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_x,CCTK_VARIABLE_REAL)
  resid_x = sqrt(resid_x)
  call norm2_dv(ext,MRsy,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
!  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_y,CCTK_VARIABLE_REAL)
  resid_y = sqrt(resid_y)
  dIntegral = dIntegral*multfactor
  call norm2_dv(ext,MRsz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
!  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_z,CCTK_VARIABLE_REAL)
  resid_z = sqrt(resid_z)
  call norm2_dv(ext,MNorm,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
!  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_mom_norm,CCTK_VARIABLE_REAL)
  resid_mom_norm = sqrt(resid_mom_norm)
  export_data(4)=sqrt(multfactor)*resid_x/resid_mom_norm
  export_data(5)=sqrt(multfactor)*resid_y/resid_mom_norm
  export_data(6)=sqrt(multfactor)*resid_z/resid_mom_norm

  !************ END MOMENTUM CONSTRAINT ************!

  data_headers(7) = "Gamx Res"
  data_headers(8) = "Gamy Res"
  data_headers(9) = "Gamz Res"
  !************ GAMMA CONSTRAINT ************!
  call gamcheck(ext,X,Y,Z,gconx,gcony,gconz,Gammax, &
       Gammay,Gammaz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Symmetry)
  call norm2_dv(ext,gconx,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gxres,CCTK_VARIABLE_REAL)
  dIntegral = dIntegral*multfactor
  export_data(7) = sqrt(multfactor)*sqrt(gxres)
  call norm2_dv(ext,gcony,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gyres,CCTK_VARIABLE_REAL)
  export_data(8) = sqrt(multfactor)*sqrt(gyres)
  call norm2_dv(ext,gconz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gzres,CCTK_VARIABLE_REAL)
  export_data(9) = sqrt(multfactor)*sqrt(gzres)
  !************ END GAMMA CONSTRAINT ************!

!  file_suffix = '001.dat'
!  filename = 'OS_rad' // file_suffix

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if
  !--------------------------------------------------------------------------------!

  filename = 'OS_rad.ron'
  num_cols = 7

  data_headers(1) = "Time"
  export_data(1) = CCTK_TIME

  data_headers(2) = "Proper time"                                                                                                          
  export_data(2)=proper_time
  
!# time         Gamx_max         rho_b_max        Ang_Mom          lapse_min        Kpoly_avg        Rc               Rp               max(b^2/rho0)    max(b^2/P)       
  data_headers(3) = "Gammaxmax"
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"bssn::Gammax")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(3),1,index)

  data_headers(4) = "rho_bmax"
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::rho_b")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(4),1,index)

  data_headers(5) = "AngMom"
  !we use Pr as a temporary variable below:
  call jdens(ext,Pr, &
       X,Y,Z,phi,Sx,Sy,trK, &
       Axx,Axy,Axz,Ayy,Ayz,Azz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Symmetry) !,imin,jmin,kmin,sNx,sNz, &
  !       PI_SYMM_K_gz,PI_SYMM_gupxx_gz,PI_SYMM_gupxy_gz,PI_SYMM_gupxz_gz, &
  !       PI_SYMM_gupyy_gz,PI_SYMM_gupyz_gz,PI_SYMM_gupzz_gz)
  call mass_0_integral(ext,dIntegral,X,Y,Z,Pr,Symmetry, &
       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(5),CCTK_VARIABLE_REAL)

  data_headers(6) = "lapse_min"
  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"lapse::lapm1")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(6),1,index)
  export_data(6) = export_data(6) + 1.D0

  data_headers(7) = "Rp"
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::rho_b")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,rho_tiny,1,index)
  rho_tiny = rho_tiny * 1.D-3
  call find_equatorial_radius(ext,cctkgh,X,Y,Z,dX,dY,dZ,rho_tiny,rho_star,phi,gxx,gxy,gyy,Symmetry,export_data(7),CCTK_MyProc(CCTKGH))

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !--------------------------------------------------------------------------------!

  filename = 'OS_rad.xon'
  num_cols = 6

  data_headers(1) = "Time"
  export_data(1) = CCTK_TIME

  data_headers(2) = "Proper time"                                                                                                          
  export_data(2)=proper_time
  
  data_headers(3) = "Xcen"
  data_headers(4) = "Ycen"
  data_headers(5) = "Zcen"
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
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(3),CCTK_VARIABLE_REAL)
     
     call wavg(ext,dIntegral,X,Y,Z,Y,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(4),CCTK_VARIABLE_REAL)

     call wavg(ext,dIntegral,X,Y,Z,Z,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(5),CCTK_VARIABLE_REAL)
  else if (Symmetry==EQUATORIAL) then
     call wavg(ext,dIntegral,X,Y,Z,X,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(3),CCTK_VARIABLE_REAL)
     
     call wavg(ext,dIntegral,X,Y,Z,Y,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(4),CCTK_VARIABLE_REAL)
     
     export_data(5) = 0.D0
  else if (Symmetry==AXISYM .and. Z(2,2,2) < 0.0) then
     export_data(3) = 0.D0
     export_data(4) = 0.D0

     call wavg(ext,dIntegral,X,Y,Z,Z,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(5),CCTK_VARIABLE_REAL)

  else 
     export_data(3) = 0.D0
     export_data(4) = 0.D0
     export_data(5) = 0.D0
  end if

  data_headers(6) = "Rcen"
  export_data(6) = sqrt(export_data(3)*export_data(3) + export_data(4)*export_data(4) + export_data(5)*export_data(5))

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  
  !--------------------------------------------------------------------------------!

end subroutine OS_rad_diag_integrals
