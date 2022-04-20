#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!---------------------
! Do bns diagnostics
!---------------------
subroutine bns_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,ext_noghost
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20
  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,xmax
  real*8                                   :: RestMass,ADMMass,Jtot,temp,resid,resid_norm
  real*8                                   :: Mass_GW,hpx_alt
  real*8                                   :: hplus,hcross
  real*8                                   :: resid_x,resid_y,resid_z,resid_mom_norm
  real*8                                   :: dIntegral,pmV,pV
  real*8                                   :: gxres,gyres,gzres
  real*8                                   :: multfactor
  real*8                                   :: rho_tiny
  real*8                                   :: rbr_GW,drbdr_GW,ddrbddr_GW
  integer                                  :: interpolate_order
  integer :: ii,l,m,nmodes
  integer                                  :: adjimin, adjjmin, adjkmin
  integer                                  :: adjimax, adjjmax, adjkmax
  CCTK_REAL :: reduction_value
  CCTK_INT :: red_tmp
  real*8  :: surf_radius2
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

  !  call setup_global_coord_arrays(CCTK_PASS_FTOF)

  !header_flag == 1 -> write file headers
  if(CCTK_TIME == 0.D0) then
     header_flag = 1
  else
     header_flag = 0
  end if

  !--------------------------------------------------------------------------------!

  filename = 'bns.mon'
  num_cols = 5
  export_data = 0.D0

  data_headers(1) = '# Time'
  export_data(1) = CCTK_TIME

  data_headers(2) = 'Phicent'
  call interp_center(cctkGH,cctk_nghostzones,phi,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(2))

  data_headers(3) = 'Kcent'
  call interp_center(cctkGH,cctk_nghostzones,trK,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(3)) 

  data_headers(4) = 'Mass_sur1'
  out_surf_int_sum = 0.d0
  if (fisheye_enable==1) then 
     rbr_GW = radius_GW_phys/radius_GW
     drbdr_GW = dR_GW 
     ddrbddr_GW = ddR_GW
  else
     rbr_GW = 1.d0
     drbdr_GW = 1.d0
     ddrbddr_GW = 0.d0
  end if
  call surf_Mass3_fisheye(cctkGH,cctk_lsh,cctk_gsh,cctk_nghostzones, & 
       X,Y,Z,dX,dY,dZ,ntot, &
       radius_GW,N_theta,N_phi,Xlocal1d,Ylocal1d,Zlocal1d, &
       Xglobal,Yglobal,Zglobal,phi,gupxx,gupxy,gupxz,gupyy, &
       gupyz,gupzz,phix,phiy,phiz,Gammax,Gammay,Gammaz,rbr_GW,drbdr_GW, &
       ddrbddr_GW,out_surf_int_sum,sym_factor,dcostheta,dphi)
  export_data(4) = out_surf_int_sum

  data_headers(5) = 'J_surf1' 
  out_surf_int_sum = 0.d0
  call  J_surf_integral2(cctkGH,cctk_lsh,cctk_gsh,cctk_nghostzones, &
       X,Y,Z,dX,dY,dZ,ntot,radius_GW,N_theta,N_phi, &
       Xlocal1d,Ylocal1d,Zlocal1d, &
       Xglobal,Yglobal,Zglobal,phi,gupxx,gupxy,gupxz,gupyy, &
       gupyz,gupzz,Axx,Axy,Axz,Ayy,Ayz,Azz, &
       out_surf_int_sum,sym_factor,dcostheta,dphi)
  export_data(5) = out_surf_int_sum

  if (radius_GW2 .gt. 0.d0) then 
     num_cols = num_cols +2
     if (fisheye_enable==1) then
        rbr_GW = radius_GW2_phys/radius_GW2
        drbdr_GW = dR2_GW
        ddrbddr_GW = ddR2_GW
     else
        rbr_GW = 1.d0
        drbdr_GW = 1.d0
        ddrbddr_GW = 0.d0
     end if
     data_headers(num_cols-1) = 'Mass_sur2'
     out_surf_int_sum = 0.d0
     call surf_Mass3_fisheye(cctkGH,cctk_lsh,cctk_gsh,cctk_nghostzones, & 
          X,Y,Z,dX,dY,dZ,ntot, &
          radius_GW2,N_theta,N_phi,Xlocal1d,Ylocal1d,Zlocal1d, &
          Xglobal,Yglobal,Zglobal,phi,gupxx,gupxy,gupxz,gupyy, &
          gupyz,gupzz,phix,phiy,phiz,Gammax,Gammay,Gammaz,rbr_GW,drbdr_GW, &
          ddrbddr_GW,out_surf_int_sum,sym_factor,dcostheta,dphi)
     export_data(num_cols-1) = out_surf_int_sum

     data_headers(num_cols) = 'J_surf2'
     out_surf_int_sum = 0.d0
     call  J_surf_integral2(cctkGH,cctk_lsh,cctk_gsh,cctk_nghostzones, &
          X,Y,Z,dX,dY,dZ,ntot,radius_GW2,N_theta,N_phi, &
          Xlocal1d,Ylocal1d,Zlocal1d, &
          Xglobal,Yglobal,Zglobal,phi,gupxx,gupxy,gupxz,gupyy, &
          gupyz,gupzz,Axx,Axy,Axz,Ayy,Ayz,Azz, &
          out_surf_int_sum,sym_factor,dcostheta,dphi)
     export_data(num_cols) = out_surf_int_sum
  end if

  if (radius_GW3 .gt. 0.d0) then
     num_cols = num_cols +2
     if (fisheye_enable==1) then
        rbr_GW = radius_GW3_phys/radius_GW3
        drbdr_GW = dR3_GW
        ddrbddr_GW = ddR3_GW
     else
        rbr_GW = 1.d0
        drbdr_GW = 1.d0
        ddrbddr_GW = 0.d0
     end if
     data_headers(num_cols-1) = 'Mass_sur3'
     out_surf_int_sum = 0.d0
     call surf_Mass3_fisheye(cctkGH,cctk_lsh,cctk_gsh,cctk_nghostzones, &
          X,Y,Z,dX,dY,dZ,ntot, &
          radius_GW3,N_theta,N_phi,Xlocal1d,Ylocal1d,Zlocal1d, &
          Xglobal,Yglobal,Zglobal,phi,gupxx,gupxy,gupxz,gupyy, &
          gupyz,gupzz,phix,phiy,phiz,Gammax,Gammay,Gammaz,rbr_GW,drbdr_GW, &
          ddrbddr_GW,out_surf_int_sum,sym_factor,dcostheta,dphi)
     export_data(num_cols-1) = out_surf_int_sum

     data_headers(num_cols) = 'J_surf3'
     out_surf_int_sum = 0.d0
     call  J_surf_integral2(cctkGH,cctk_lsh,cctk_gsh,cctk_nghostzones, &
          X,Y,Z,dX,dY,dZ,ntot,radius_GW3,N_theta,N_phi, &
          Xlocal1d,Ylocal1d,Zlocal1d, &
          Xglobal,Yglobal,Zglobal,phi,gupxx,gupxy,gupxz,gupyy, &
          gupyz,gupzz,Axx,Axy,Axz,Ayy,Ayz,Azz, &
          out_surf_int_sum,sym_factor,dcostheta,dphi)
     export_data(num_cols) = out_surf_int_sum
  end if

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !--------------------------------------------------------------------------------!

  filename = 'bns.don'

  num_cols = 10
  export_data = 0.D0

  data_headers(1) = '# Time'
  export_data(1) = CCTK_TIME

  data_headers(2) = 'M0'
  call mass_0_integral(ext, dIntegral, &
       X, Y, Z, &
       rho_star, Symmetry, &
       adjimin, adjimax, &
       adjjmin, adjjmax, &
       adjkmin, adjkmax,CCTK_MyProc(cctkGH)) 
  !call CCTK_VarIndex(index,"bns::dIntegral")
  !  call CCTK_ReductionHandle(handle,"sum")
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(2),CCTK_VARIABLE_REAL)
  !  call CCTK_Reduce (ierr,cctkGH, 0,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)

  data_headers(3) = 'Eint'
  call minternal_hybrid(ext,dIntegral, X,Y,Z, &
       rho_b_atm,rho_star,rho_b,P, &
       neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
       Symmetry,adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax) 
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
  if (fisheye_enable==1) then
     call surf_Mass2_fisheye(CCTK_PASS_FTOF)
  else
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
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
  end if
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

  !--------------------------------------------------------------------------------!
  filename = 'bns.jon'
  num_cols = 2

  data_headers(1) = "# Time"
  export_data(1) = CCTK_TIME

  data_headers(2) = "AngMom"
  !we use Pr as a temporary variable below:
  call jdens_cpp(cctkGH,cctk_lsh,cctk_nghostzones, Symmetry, &
       Pr,  &
       dx,dy,dz, &
       X,Y,Z,  &
       phi, Sx, Sy, trK, &
       Axx, Axy, Axz, Ayy, Ayz, Azz,  &
       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
  !call jdens(ext,Pr, &
  !     X,Y,Z,phi,Sx,Sy,trK, &
  !     Axx,Axy,Axz,Ayy,Ayz,Azz, &
  !     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
  !     Symmetry) 
  call mass_0_integral(ext,dIntegral,X,Y,Z,Pr,Symmetry, &
       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(2),CCTK_VARIABLE_REAL)

  if(1==0) then
     !if(cctk_iteration.gt.10000) then
     data_headers(3) = "J_out"
     call J_surf_integral(CCTK_PASS_FTOF)
     export_data(3) = out_surf_int_sum
     !end if

     !if(cctk_iteration.gt.10000) then

     data_headers(4) = "Mass_sur"
     call surf_Mass2_fisheye(CCTK_PASS_FTOF)
     export_data(4) = out_surf_int_sum

     data_headers(5) = "Komar Mass"
     !  write(*,*) "Now computing",data_headers(5)
     call Komar_Mass(CCTK_PASS_FTOF)
     !  write(*,*) "Finished computing",data_headers(5)
     !end if

     export_data(5) = out_surf_int_sum
     data_headers(6) = "out_M0_flux"
     data_headers(7) = "out_fluid_Eng_flux"
     data_headers(8) = "out_EM_Energy_flux"
     data_headers(9) = "out_fluid_J_flux"
     data_headers(10) = "out_EM_J_flux"
     !  write(*,*) "Now computing",data_headers(6)
     call fluxes(CCTK_PASS_FTOF)
     !  write(*,*) "finished computing",data_headers(6)
     export_data(6) = F_M0
     export_data(7) = F_E_fluid
     export_data(8) = F_E_em
     export_data(9) = F_J_fluid
     export_data(10) = F_J_em  

     data_headers(11) = "rho_isoflux"
     out_surf_int_sum = 0.D0
     !  write(*,*) "Now computing",data_headers(11)
     call rhoflux(CCTK_PASS_FTOF)
     !  write(*,*) "finished computing",data_headers(11)
     export_data(11) = out_surf_int_sum
  end if
  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if
  !--------------------------------------------------------------------------------!

  filename = 'bns.con'
  num_cols = 8

  data_headers(1) = "# Time"
  export_data(1) = CCTK_TIME


  data_headers(2) = "Haml Res"
  !************ HAMILTONIAN CONSTRAINT ************!
  !
  ! Next the L2 norm.  Note I don't need to Sync since I don't care
  ! about the ghost zones.
  !
  dIntegral = 0.D0
  call norm2_dv(ext,PsiRes,X,Y,Z,dIntegral,adjimin, &
       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
  resid = sqrt(resid)

  dIntegral = 0.D0
  call norm2_dv(ext,PsiNorm,X,Y,Z,dIntegral,adjimin, &
       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
  resid_norm = sqrt(resid_norm)
  export_data(2) = resid/resid_norm
  !write(*,*) "Resid = ",resid,resid_norm,resid/resid_norm
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
  export_data(3)=resid_x/resid_mom_norm
  export_data(4)=resid_y/resid_mom_norm
  export_data(5)=resid_z/resid_mom_norm

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
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gxres,CCTK_VARIABLE_REAL)
  export_data(6) = sqrt(gxres)
  call norm2_dv(ext,gcony,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gyres,CCTK_VARIABLE_REAL)
  export_data(7) = sqrt(gyres)
  call norm2_dv(ext,gconz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gzres,CCTK_VARIABLE_REAL)
  export_data(8) = sqrt(gzres)

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if
  !************ END GAMMA CONSTRAINT ************!
  !--------------------------------------------------------------------------------!

  filename = 'bns.ron'
  num_cols = 7

  data_headers(1) = "# Time"
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
  !call jdens(ext,Pr, &
  !     X,Y,Z,phi,Sx,Sy,trK, &
  !     Axx,Axy,Axz,Ayy,Ayz,Azz, &
  !     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
  !     Symmetry) 
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

  data_headers(7) = "betam1" 
  ! set rho_br = exp(6 phi), rho_bl = exp(4 phi), v02l = b^2/2
  rho_br = exp(phi)
  rho_bl = rho_br**4
  rho_br = rho_br**6
  call bns_compute_b2(ext,v02l,rho_bl,rho_br,gxx,gxy,gxz,gyy,gyz,gzz, &
       Bx,By,Bz,rho_star,&
       rho_b,lapm1,vx,vy,vz,shiftx,shifty,shiftz)
  call wavg(ext,dIntegral,X,Y,Z,v02l,rho_br,Symmetry, &
       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,pmV,CCTK_VARIABLE_REAL)
  call wavg(ext,dIntegral,X,Y,Z,P,rho_br,Symmetry, &
       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,pV,CCTK_VARIABLE_REAL)
  export_data(7) = pmV/pV

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !--------------------------------------------------------------------------------!

  filename = 'bns.xon'
  num_cols = 5

  data_headers(1) = "# Time"
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
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(2),CCTK_VARIABLE_REAL)

     call wavg(ext,dIntegral,X,Y,Z,Y,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(3),CCTK_VARIABLE_REAL)

     call wavg(ext,dIntegral,X,Y,Z,Z,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(4),CCTK_VARIABLE_REAL)
  else if (Symmetry==EQUATORIAL) then
     call wavg(ext,dIntegral,X,Y,Z,X,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(2),CCTK_VARIABLE_REAL)

     call wavg(ext,dIntegral,X,Y,Z,Y,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
     dIntegral = dIntegral*multfactor/RestMass
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(3),CCTK_VARIABLE_REAL)

     export_data(4) = 0.D0
  else if (Symmetry==AXISYM .and. Z(2,2,2) < 0.0) then
     export_data(2) = 0.D0
     export_data(3) = 0.D0

     call wavg(ext,dIntegral,X,Y,Z,Z,rho_star,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
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

  !gw_amplitudes = 0.D0
  hplus = 0.d0
  hcross = 0.d0

  nmodes = ubound(momentsr,1)

  call gw_wave_flux(cctkGH,radius_GW,Mass_GW, &
       dT,momentsr,momentsi, &
       momentsr_old,momentsi_old,odd_momentsr,odd_momentsi, &
       odd_momentsr_old,odd_momentsi_old,int_momentsr,int_momentsi, &
       hplus,hcross,theta_GW,phi_GW, &
       Zmin,dZ)

  filename = 'bns.mome1'
  num_cols = 2*nmodes + 1

  data_headers(1) = "# Time"
  export_data(1) = CCTK_TIME

  do ii=1,nmodes
     data_headers(2*ii) = " " 
     data_headers(2*ii+1) = " " 
     export_data(2*ii) = momentsr(ii)
     export_data(2*ii+1) = momentsi(ii)
  end do

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if
  !--------------------------------------------------------------------------------!

  filename = 'bns.momo1'
  num_cols = 2*nmodes + 1

  data_headers(1) = "# Time"
  export_data(1) = CCTK_TIME

  do ii=1,nmodes
     data_headers(2*ii) = " "
     data_headers(2*ii+1) = " "
     export_data(2*ii) = odd_momentsr(ii)
     export_data(2*ii+1) = odd_momentsi(ii)
  end do

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if
  !--------------------------------------------------------------------------------!

  if (radius_GW2 .gt. 0.d0) then 
     call gw_wave_flux(cctkGH,radius_GW2,Mass_GW, &
          dT,momentsr,momentsi, &
          momentsr_old,momentsi_old,odd_momentsr,odd_momentsi, &
          odd_momentsr_old,odd_momentsi_old,int_momentsr,int_momentsi, &
          hplus,hcross,theta_GW,phi_GW, &
          Zmin,dZ)

     filename = 'bns.mome2'
     num_cols = 2*nmodes + 1
     data_headers(1) = "# Time"
     export_data(1) = CCTK_TIME

     do ii=1,nmodes
        data_headers(2*ii) = " "
        data_headers(2*ii+1) = " "
        export_data(2*ii) = momentsr(ii)
        export_data(2*ii+1) = momentsi(ii)
     end do

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     filename = 'bns.momo2'
     num_cols = 2*nmodes + 1

     data_headers(1) = "# Time"
     export_data(1) = CCTK_TIME

     do ii=1,nmodes
        data_headers(2*ii) = " "
        data_headers(2*ii+1) = " "
        export_data(2*ii) = odd_momentsr(ii)
        export_data(2*ii+1) = odd_momentsi(ii)
     end do

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

  end if

  !--------------------------------------------------------------------------------!

  if (radius_GW3 .gt. 0.d0) then
     call gw_wave_flux(cctkGH,radius_GW3,Mass_GW, &
          dT,momentsr,momentsi, &
          momentsr_old,momentsi_old,odd_momentsr,odd_momentsi, &
          odd_momentsr_old,odd_momentsi_old,int_momentsr,int_momentsi, &
          hplus,hcross,theta_GW,phi_GW, &
          Zmin,dZ)

     filename = 'bns.mome3'
     num_cols = 2*nmodes + 1
     data_headers(1) = "# Time"
     export_data(1) = CCTK_TIME

     do ii=1,nmodes
        data_headers(2*ii) = " "
        data_headers(2*ii+1) = " "
        export_data(2*ii) = momentsr(ii)
        export_data(2*ii+1) = momentsi(ii)
     end do

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     filename = 'bns.momo3'
     num_cols = 2*nmodes + 1

     data_headers(1) = "# Time"
     export_data(1) = CCTK_TIME

     do ii=1,nmodes
        data_headers(2*ii) = " "
        data_headers(2*ii+1) = " "
        export_data(2*ii) = odd_momentsr(ii)
        export_data(2*ii+1) = odd_momentsi(ii)
     end do

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

  end if

end subroutine bns_diagnostics

! Compute b^2/2
subroutine bns_compute_b2(ext,sb2f,psi4,psi6,gxx,gxy,gxz,gyy,gyz,gzz, & 
     Bx,By,Bz,rho_star,& 
     rho_b,alpm1,vx,vy,vz,betax,betay,betaz)
  implicit none
  integer, dimension(3)					:: ext
  real*8, dimension(ext(1),ext(2),ext(3))		:: gxx,gxy,gxz,alpm1
  real*8, dimension(ext(1),ext(2),ext(3))               :: psi6,psi4,sb2f
  real*8, dimension(ext(1),ext(2),ext(3))               :: gyy,gyz,gzz,rho_star
  real*8, dimension(ext(1),ext(2),ext(3))               :: Bx,By,Bz,vx,vy,vz
  real*8, dimension(ext(1),ext(2),ext(3))               :: betax,betay,betaz,rho_b
  ! Aux vars
  real*8, dimension(ext(1),ext(2),ext(3))               :: u0,B_x,B_y,B_z
  real*8, parameter :: pi = 3.14159265358979323846d0
  !
  u0 = rho_star / (rho_b*(1.d0+alpm1)*psi6)
  B_x = psi4*(gxx*Bx + gxy*By + gxz*Bz)
  B_y = psi4*(gxy*Bx + gyy*By + gyz*Bz)
  B_z = psi4*(gxz*Bz + gyz*By + gzz*Bz)

  sb2f = ( Bx*B_x + By*B_y + Bz*B_z + & 
       ( u0*( B_x*(vx+betax) + B_y*(vy+betay) + B_z*(vz+betaz) ) )**2 ) & 
       / (8.d0*pi* ( (1.d0+alpm1)*u0 )**2 )

end subroutine bns_compute_b2
