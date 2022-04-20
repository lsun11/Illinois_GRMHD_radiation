!THIS IS BROKEN.  MUST BE FIXED







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
  real*8, dimension(20)                    :: export_data
  character, dimension(20)                 :: data_headers*20
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
  integer :: ii,jj
  integer                                  :: adjimin, adjjmin, adjkmin
  integer                                  :: adjimax, adjjmax, adjkmax
  CCTK_REAL :: reduction_value
  CCTK_INT :: red_tmp
  integer :: header_flag,handle,dummy
  integer :: index,num_cols
  integer :: ierr,myproc_rank
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8, dimension(4,3)                     :: pointcoords
  real*8, dimension(narr,3)                  :: pos
  real*8, dimension(coarse_res,3)            :: pos_int
  real*8, dimension(4)                       :: rho_b_int,P_int,E_rad_int
  real*8, dimension(4)                       :: F_radx_int
  real*8, dimension(4)                       :: areal_radius_target
  real*8                                     :: areal_radius_old,areal_radius_new
  real*8                                     :: frac,r_interp
  real *8, dimension(4)                      :: blah
  real *8, dimension(narr)                   :: areal_radius,gxx_pcle,gxz_pcle,gzz_pcle,phi_pcle,rho_b_pcle,u0_pcle,E_rad_pcle,Ec_pcle,eta,Q,F_pcle,Fc_pcle
  real *8, dimension(narr)                   :: z_anal,rho_b_anal,Ec_anal,Fc_anal,E_anal,F_anal,E_rad_gradient,F_coeff
  real*8                                     :: PI
  real*8                                     :: rho_b_initial
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
  !this is the new method I am using because I don't like interpolating at the center. It can be misleading.
  if (Symmetry==AXISYM) then
     filename = 'OS_rad.mon'
     num_cols = 17
     export_data = 0.D0
     
     data_headers(1) = 'Time'
     export_data(1) = CCTK_TIME
     
     data_headers(3) = 'Phicent'
     export_data(3) = phi(2,2,2)
     
     data_headers(4) = 'Kcent'
     export_data(4) = trK(2,2,2)

     data_headers(5) = 'gxxcent-1'
     export_data(5) = gxx(2,2,2) - 1

     data_headers(6) = 'rhocent'
     !call interp_center(cctkGH,rho_b,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(6)) 
     export_data(6) = rho_b(2,2,2)

     data_headers(7) = 'lapsecent'
     export_data(7) = lapm1(2,2,2) + 1
     
     data_headers(2) = "Proper time"
     if (CCTK_TIME > 0) then
        proper_time = proper_time + dT*export_data(7)
     end if
     export_data(2)=proper_time
     
     data_headers(8) = 'Omega_c'
     Pr = (vy*X - vx*Y) / (X*X + Y*Y)
     export_data(8) = Pr(2,2,2)

     data_headers(9) = "E_radcent"
     export_data(9) = E_rad(2,2,2)
     
     data_headers(10) = 'Pcent'
     export_data(10) = P(2,2,2)

     data_headers(11) = 'mhd_st_y_cent'
     export_data(11) = mhd_st_y(2,2,2)
     
     data_headers(12) = 'vy_cent'
     export_data(12) = vy(2,2,2)
     
     data_headers(13) = 'shifty_cent'
     export_data(13) = shifty(2,2,2)
     
     data_headers(14) = 'Gammay_cent'
     export_data(14) = Gammay(2,2,2)
     
     data_headers(15) = 'gxy_cent'
     export_data(15) = gxy(2,2,2)
     
     data_headers(16) = 'gyy_cent'
     export_data(16) = gyy(2,2,2)

     data_headers(16) = 'S_rad_y_cent'
     export_data(16) = S_rad_y(2,2,2)
 
     data_headers(17) = 'F_rady_cent'
     export_data(17) = F_rady(2,2,2)

     
     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
  else
     
     !This is the old way that zach had set up, with minor additions
     filename = 'OS_rad.mon'
     num_cols = 10
     export_data = 0.D0
     
     data_headers(1) = 'Time'
     export_data(1) = CCTK_TIME
     
     data_headers(3) = 'Phicent'
     call interp_center(cctkGH,phi,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(3))
     
     
     data_headers(4) = 'Kcent'
     call interp_center(cctkGH,trK,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(4)) 
     
     data_headers(5) = 'gxxcent-1'
     call interp_center(cctkGH,gxx,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(5)) 
     export_data(5) = export_data(5) - 1
     
     data_headers(6) = 'rhocent'
     call interp_center(cctkGH,rho_b,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(6)) 
     
     data_headers(7) = 'lapsecent'
     call interp_center(cctkGH,lapm1,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(7)) 
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
     data_headers(9) = "E_radcent"
     
     call interp_center(cctkGH,E_rad,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(9))
     
     data_headers(10) = 'Pcent'
     call interp_center(cctkGH,P,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(10)) 
     
     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
  endif
  !--------------------------------------------------------------------------------!
  
  filename = 'OS_rad.don'

  num_cols = 10
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

!I'll fix this later



!  filename = 'OS_rad.con'
!  num_cols = 9
!  
!  data_headers(1) = "Time"
!  
!  export_data(1) = CCTK_TIME
!  
!  data_headers(2) = "Proper time"                                                                                                            
!  export_data(2)=proper_time
!  
!  data_headers(3) = "Haml Res"
!  !************ HAMILTONIAN CONSTRAINT ************!
!  PsiTau = PsiRes
!  !
!  ! Next the L2 norm.  Note I don't need to Sync since I don't care
!  ! about the ghost zones.
!  !
!  dIntegral = 0.D0
!  call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
!       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
!!  dIntegral = dIntegral*multfactor
!  call CCTK_ReductionHandle(handle,"sum")
!  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
!  resid = sqrt(resid)
!
!  dIntegral = 0.D0
!  call norm2_dv(ext,PsiNorm,X,Y,Z,dIntegral,adjimin, &
!       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry);
!!  dIntegral = dIntegral*multfactor
!  call CCTK_ReductionHandle(handle,"sum")
!  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
!  resid_norm = sqrt(resid_norm)
!  export_data(3) = sqrt(multfactor)*resid/resid_norm
!  !write(*,*) "Resid = ",resid,resid_norm,resid/resid_norm
!  !************ END HAMILTONIAN CONSTRAINT ************!
!
!  data_headers(4) = "Momx Res"
!  data_headers(5) = "Momy Res"
!  data_headers(6) = "Momz Res"
!  !************ MOMENTUM CONSTRAINT ************!
!  call norm2_dv(ext,MRsx,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
!       adjjmax,adjkmin,adjkmax,Symmetry)
!!  dIntegral = dIntegral*multfactor
!  call CCTK_ReductionHandle(handle,"sum")
!  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_x,CCTK_VARIABLE_REAL)
!  resid_x = sqrt(resid_x)
!  call norm2_dv(ext,MRsy,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
!       adjjmax,adjkmin,adjkmax,Symmetry)
!!  dIntegral = dIntegral*multfactor
!  call CCTK_ReductionHandle(handle,"sum")
!  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_y,CCTK_VARIABLE_REAL)
!  resid_y = sqrt(resid_y)
!  dIntegral = dIntegral*multfactor
!  call norm2_dv(ext,MRsz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
!       adjjmax,adjkmin,adjkmax,Symmetry)
!!  dIntegral = dIntegral*multfactor
!  call CCTK_ReductionHandle(handle,"sum")
!  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_z,CCTK_VARIABLE_REAL)
!  resid_z = sqrt(resid_z)
!  call norm2_dv(ext,MNorm,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
!       adjjmax,adjkmin,adjkmax,Symmetry)
!!  dIntegral = dIntegral*multfactor
!  call CCTK_ReductionHandle(handle,"sum")
!  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_mom_norm,CCTK_VARIABLE_REAL)
!  resid_mom_norm = sqrt(resid_mom_norm)
!  export_data(4)=sqrt(multfactor)*resid_x/resid_mom_norm
!  export_data(5)=sqrt(multfactor)*resid_y/resid_mom_norm
!  export_data(6)=sqrt(multfactor)*resid_z/resid_mom_norm
!
!  !************ END MOMENTUM CONSTRAINT ************!
!
!  data_headers(7) = "Gamx Res"
!  data_headers(8) = "Gamy Res"
!  data_headers(9) = "Gamz Res"
!!  call gamcheck(ext,X,Y,Z,gconx,gcony,gconz,Gammax, &
!       Gammay,Gammaz, &
!       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
!       Symmetry)
!  call norm2_dv(ext,gconx,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
!       adjjmax,adjkmin,adjkmax,Symmetry)
!  call CCTK_ReductionHandle(handle,"sum")
!  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gxres,CCTK_VARIABLE_REAL)
!  dIntegral = dIntegral*multfactor
!  export_data(7) = sqrt(multfactor)*sqrt(gxres)
!  call norm2_dv(ext,gcony,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
!       adjjmax,adjkmin,adjkmax,Symmetry)
!  dIntegral = dIntegral*multfactor
!  call CCTK_ReductionHandle(handle,"sum")
!  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gyres,CCTK_VARIABLE_REAL)
!  export_data(8) = sqrt(multfactor)*sqrt(gyres)
!  call norm2_dv(ext,gconz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
!       adjjmax,adjkmin,adjkmax,Symmetry)
!  dIntegral = dIntegral*multfactor
!  call CCTK_ReductionHandle(handle,"sum")
!  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gzres,CCTK_VARIABLE_REAL)
!  export_data(9) = sqrt(multfactor)*sqrt(gzres)
!  !************ END GAMMA CONSTRAINT ************!
!
!!  file_suffix = '001.dat'
!!  filename = 'OS_rad' // file_suffix
!
!  if(CCTK_MyProc(CCTKGH) .eq. 0) then
!     write(*,*) "filename: ",filename
!     write(*,*) "num_cols: ",num_cols
!     write(*,*) "data_headers: ",data_headers
!     write(*,*) "export_data: ",export_data
!     write(*,*) "header_flag: ",header_flag
!
!     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
!  end if
!


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
!  filename = 'OS_rad.bon'
!  num_cols = 16
!  
!  !note: there are issues here with ghostzones, symmetries, etc.  Use with caution!
!
!  
!
!  areal_radius_target(1)=2.d0
!  areal_radius_target(2)=4.d0
!  areal_radius_target(3)=6.d0
!  areal_radius_target(4)=8.d0
!
!
!  areal_radius_old=exp(2.d0*phi(2,2,2))*sqrt(Z(2,2,2)*Z(2,2,2)*gxx(2,2,2)-2.d0*X(2,2,2)*Z(2,2,2)*gxz(2,2,2)+X(2,2,2)*X(2,2,2)*gzz(ii,2,2))
!  areal_radius_new=areal_radius_old
!
!  rho_b_int = 0.d0
!  P_int = 0.d0
!  E_rad_int = 0.d0
!  F_radx_int = 0.d0
!
!  write(*,*) "areal_radius_old: ",areal_radius_old
!  write(*,*) "areal_radius_new: ",areal_radius_new
!
!  
!  call CCTK_VarIndex(index,'mhd_evolve::rho_b')
!  call BndCartoon2DVI(dummy, cctkGH, 0, index)
!  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_vars')
! 
!  call CCTK_VarIndex(index,'mhd_evolve::P')
!  call BndCartoon2DVI(dummy, cctkGH, 0, index)
!  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_vars')
! 
!  blah=0.d0
!  pointcoords=0.d0
!
!  do jj=1,4
!     do ii=2,ext(1)-1
!        areal_radius_new = exp(2.d0*phi(ii,2,2))*sqrt(Z(ii,2,2)*Z(ii,2,2)*gxx(ii,2,2)-2.d0*X(ii,2,2)*Z(ii,2,2)*gxz(ii,2,2)+X(ii,2,2)*X(ii,2,2)*gzz(ii,2,2))
!    
!     if (areal_radius_new.gt. areal_radius_target(jj) .and. areal_radius_old .lt. areal_radius_target(jj) .and. Zlocal1d(1) .eq. Zglobal(1)) then
!     
!        frac = (areal_radius_target(jj)- areal_radius_old)/(areal_radius_new-areal_radius_old)
!        pointcoords(jj,1) = X(ii-1,2,2)+frac*dX
!        pointcoords(jj,2) = Y(ii-1,2,2)
!        pointcoords(jj,3) = Z(ii-1,2,2)
!        exit
!     else
!        areal_radius_old=areal_radius_new
!     endif
!  end do
!end do
! 
!call CCTK_ReductionHandle(handle,"sum")
!call  CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,pointcoords,pointcoords,4,CCTK_VARIABLE_REAL)
!
!call interpolate_pointset(cctkGH,4,2,ext,global_ext,dX,dY,dZ, &
!     Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,rho_b,rho_b_int) 
!call interpolate_pointset(cctkGH,4,2,ext,global_ext,dX,dY,dZ, &
!     Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,P,P_int) 
!call interpolate_pointset(cctkGH,4,2,ext,global_ext,dX,dY,dZ, &
!     Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,E_rad,E_rad_int) 
!call interpolate_pointset(cctkGH,4,2,ext,global_ext,dX,dY,dZ, &
!     Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,F_radx,F_radx_int) !
!
!call CCTK_ReductionHandle(handle,"sum")
!
!!call  CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,rho_b_int,rho_b_int,4,CCTK_VARIABLE_REAL)
!!call  CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,P_int,P_int,4,CCTK_VARIABLE_REAL)
!!call  CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,E_rad_int,E_rad_int,4,CCTK_VARIABLE_REAL)
!!call  CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,F_radx_int,F_radx_int,4,CCTK_VARIABLE_REAL)
!
!data_headers(1) = "rho_b at r=2"                                             
!export_data(1)=rho_b_int(1)
!data_headers(2) = "rho_b at r=4"                                             
!export_data(2)=rho_b_int(2)
!data_headers(3) = "rho_b at r=6"                                             
!export_data(3)=rho_b_int(3)
!data_headers(4) = "rho_b at r=8"                                             
!export_data(4)=rho_b_int(4)
!data_headers(5) = "P at r=2"                                             
!export_data(5)=P_int(1)
!data_headers(6) = "P at r=4"                                             
!export_data(6)=P_int(2)
!data_headers(7) = "P at r=6"                                             
!export_data(7)=P_int(3)
!data_headers(8) = "P at r=8"                                             
!export_data(8)=P_int(4)
!data_headers(9) = "E_rad at r=2"                                             
!export_data(9)=E_rad_int(1)
!data_headers(10) = "E_rad at r=4"                                             
!export_data(10)=E_rad_int(2)
!data_headers(11) = "E_rad at r=6"                                             
!export_data(11)=E_rad_int(3)
!data_headers(12) = "E_rad at r=8"                                             
!export_data(12)=E_rad_int(4)
!data_headers(13) = "F_radx at r=2"                                             
!export_data(13)=F_radx_int(1)
!data_headers(14) = "F_radx at r=4"                                             
!export_data(14)=F_radx_int(2)
!data_headers(15) = "F_radx at r=6"                                             
!export_data(15)=F_radx_int(3)
!data_headers(16) = "F_radx at r=8"                                             
!export_data(16)=F_radx_int(4)
!
! if(CCTK_MyProc(CCTKGH) .eq. 0) then
!     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
!  end if
!!endif
!!--------------------------------------------------------------------------------!
  filename = 'OS_rad.bbon'

  num_cols = 17

  rho_b_initial=1.d0/(4.d0/3.d0*PI*RoM**3)

  do ii=1,3
     do jj=1,narr
        pos(jj,ii) = coord(jj,ii+1)
     end do
  end do
 
  call interpolate_pointset(cctkGH,narr,2,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos,gxx,gxx_pcle) 
  call interpolate_pointset(cctkGH,narr,2,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos,gxz,gxz_pcle) 
  call interpolate_pointset(cctkGH,narr,2,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos,gzz,gzz_pcle) 
  call interpolate_pointset(cctkGH,narr,2,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos,phi,phi_pcle) 
  call interpolate_pointset(cctkGH,narr,2,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos,rho_b,rho_b_pcle)
  call interpolate_pointset(cctkGH,narr,2,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos,E_rad,E_rad_pcle)
  call interpolate_pointset(cctkGH,narr,2,ext,global_ext,dX,dY,dZ, &
       Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos,F,F_pcle)
 ! call interpolate_pointset(cctkGH,narr,2,ext,global_ext,dX,dY,dZ, &
 !      Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos,u0,u0_pcle)


  areal_radius=exp(2.d0*phi_pcle)*sqrt(coord(:,4)*coord(:,4)*gxx_pcle-2.d0*coord(:,2)*coord(:,4)*gxz_pcle+coord(:,2)*coord(:,2)*gzz_pcle)

  do ii=1,narr
     call OS_rad_analytic_new(F_coeff(ii),z_anal(ii),areal_radius(ii),&
          coord(ii,1),RoM,1.d0,1.d0,opt_depth_a,opt_depth_s,eta(ii),Q(ii),&
          rho_b_anal(ii),Ec_anal(ii),Fc_anal(ii),E_anal(ii),F_anal(ii))
     
     Ec_pcle(ii) = E_rad_pcle(ii)*Q(ii)**4
     Fc_pcle(ii) = F_pcle(ii)*Q(ii)**4

     E_rad_gradient(1)=F_coeff(2)*(E_rad_pcle(2)-E_rad_pcle(1))/(z_anal(2)-z_anal(1))
     if (ii .gt. 1) E_rad_gradient(ii)=F_coeff(ii)*(E_rad_pcle(ii)-E_rad_pcle(ii-1))/(z_anal(ii)-z_anal(ii-1))
     
  end do

  data_headers(1) = "Coord time"                                             
  export_data(1)=CCTK_TIME

  do ii=0,narr-1,narr/2
     data_headers(8*ii/(narr/2)+2) = "tau at pcle"                        
     export_data(8*ii/(narr/2)+2)=coord(narr-ii,1)
     data_headers(8*ii/(narr/2)+3) = "tc/tcoll at pcle"        
     export_data(8*ii/(narr/2)+3)=eta(narr-ii)/PI
     data_headers(8*ii/(narr/2)+4) = "R_a at pcle"
     export_data(8*ii/(narr/2)+4)=areal_radius(narr-ii)
     data_headers(8*ii/(narr/2)+5) = "rho_b at pcle"   
     export_data(8*ii/(narr/2)+5)=rho_b_pcle(narr-ii)
     data_headers(8*ii/(narr/2)+6) = "Fc at pcle"   
     export_data(8*ii/(narr/2)+6)=Fc_pcle(narr-ii)
     data_headers(8*ii/(narr/2)+7) = "Fc_anal at pcle"   
     export_data(8*ii/(narr/2)+7)=Fc_anal(narr-ii)*rho_b_initial*PoRho/Po4PiB
     data_headers(8*ii/(narr/2)+8) = "F at pcle"   
     export_data(8*ii/(narr/2)+8)=F_pcle(narr-ii)
     data_headers(8*ii/(narr/2)+9) = "F_anal at pcle"   
     export_data(8*ii/(narr/2)+9)=F_anal(narr-ii)*rho_b_initial*PoRho/Po4PiB
 !    data_headers(3*ii+2) = "u0 at particle"   
 !    export_data(3*ii+2)= u0_pcle(ii)
  end do

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if
  

!--------------------------------------------------------------------------------!

!particle tracer stuff

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
   
   write(15,'(t2,a6,t17,a1,t31,a5,t45,a1,t59,a1,t73,a1,t87,a6,t101,a2,t115,a7,t129,a2,t143,a7)') '"coord','z','rho_b','P','E','F','grad_E','Ec','Ec_anal','Fc','Fc_anal'
   do ii=1,narr
      write(15,'(12e18.10)') coord(ii,2),z_anal(ii),rho_b_anal(ii)*rho_b_initial,&
           rho_b_anal(ii)*E_anal(ii)**0.25d0*PoRho*rho_b_initial,&
           E_anal(ii)*rho_b_initial*PoRho/Po4PiB,&
           F_anal(ii)*rho_b_initial*PoRho/Po4PiB,&
           E_rad_gradient(ii),&
           Ec_pcle(ii),Ec_anal(ii)*rho_b_initial*PoRho/Po4PiB,Fc_pcle(ii),Fc_anal(ii)*rho_b_initial*PoRho/Po4PiB,F_pcle(ii)
   end do
endif
close(15)
endif

!----------------------------------------------------------------------------------!

!convergence testing stuff
write(*,*) "out_every: ",out_every
write(*,*) "global_ext(1): ",global_ext(1)
write(*,*) "coarse_res: ", coarse_res

if (mod(CCTK_iteration,out_every*(global_ext(1)-1)/coarse_res) .eq. 0) then
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
   write(16,*) ""
   write(*,*) "global_ext(1): ",global_ext
   do ii=1,coarse_res
      pos_int(ii,1)=(ii-0.5d0)*dX*(global_ext(1)-1)/coarse_res 
      pos_int(ii,2)=0.d0 
      pos_int(ii,3)=0.5d0*dZ*(global_ext(3)-1)/coarse_res 
   end do

   rho_b_conv=0.d0
   P_conv=0.d0
   E_rad_conv=0.d0
   F_radx_conv=0.d0

   !call CCTK_ReductionHandle(handle,"sum")
   !call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,pos_int,pos_int,coarse_res,CCTK_VARIABLE_REAL)
   
   call interpolate_pointset(cctkGH,coarse_res,2,ext,global_ext,dX,dY,dZ, &
        Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos_int,rho_b,rho_b_conv) 
   call interpolate_pointset(cctkGH,coarse_res,2,ext,global_ext,dX,dY,dZ, &
        Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos_int,P,P_conv) 
   call interpolate_pointset(cctkGH,coarse_res,2,ext,global_ext,dX,dY,dZ, &
        Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos_int,E_rad,E_rad_conv) 
   call interpolate_pointset(cctkGH,coarse_res,2,ext,global_ext,dX,dY,dZ, &
        Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pos_int,F_radx,F_radx_conv) 
   
   !call  CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,rho_b_conv,rho_b_conv,coarse_res,CCTK_VARIABLE_REAL)
   !call  CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,P_conv,P_conv,coarse_res,CCTK_VARIABLE_REAL)
   !call  CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,E_rad_conv,E_rad_conv,coarse_res,CCTK_VARIABLE_REAL)
   !call  CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,F_radx_conv,F_radx_conv,coarse_res,CCTK_VARIABLE_REAL)

   if(CCTK_MyProc(CCTKGH) .eq. 0) then
      do ii=1,coarse_res
         write(16,'(5e18.10)') pos_int(ii,1),rho_b_conv(ii),P_conv(ii),E_rad_conv(ii),F_radx_conv(ii)
      end do
  !    write(16,*) ""
      close(16)
   endif
endif



end subroutine OS_rad_diag_integrals
