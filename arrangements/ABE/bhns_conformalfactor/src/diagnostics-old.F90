#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!---------------------
! Do bhns diagnostics
!---------------------
subroutine bhns_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,ext_noghost
  real*8, dimension(20)                    :: export_data
  character, dimension(20)                 :: data_headers*20
  character                                :: filename*25,file_suffix*10
  real*8                                   :: dT,dX,dY,dZ,xmax,ymax,zmax,old_surf_radius,r_AH
  real*8                                   :: RestMass,ADMMass,Jtot,temp,resid,resid_norm
  real*8                                   :: resid_x,resid_y,resid_z,resid_mom_norm
  real*8                                   :: dIntegral
  real*8                                   :: gxres,gyres,gzres
  real*8                                   :: multfactor
  real*8                                   :: detmin,detmax,Kmin,Kmax
  integer                                  :: interpolate_order

  !Needed for GW extraction:
  real*8                                   :: Mass_GW,hplus,hcross,moncrief_gw_extr_radius
  real*8                                   :: psi4r,psi4i
  integer                                  :: jj,nmodes_moncrief
  real*8,dimension(2)                      :: gw_amplitudes
  !End GW extraction variables  

  real*8                                   :: x_coord,y_coord,z_coord,variable_extremum_value

  integer                                  :: adjimin, adjjmin, adjkmin
  integer                                  :: adjimax, adjjmax, adjkmax
  CCTK_REAL :: reduction_value
  CCTK_INT :: red_tmp
  integer :: header_flag,handle,dummy
  integer :: index,num_cols
  integer :: ierr,myproc_rank,ii
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(trA_detg_enforce.ne.2) then
     !First perform sanity check.  Already done if trA_detg_enforce==2!
     call sanitycheck_restore_Aij(cctkGH,cctk_lsh,gxx,gxy,gxz,gyy,gyz,gzz, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
  end if

  if(MOD(cctk_iteration,out_every)==0) then
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

     filename = 'bhns.mon'
     num_cols = 17
     export_data = 0.D0

     data_headers(1) = 'Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Phicent'
     call interp_center(cctkGH,cctk_nghostzones,phi,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(2))

     data_headers(3) = 'Kcent'
     call interp_center(cctkGH,cctk_nghostzones,trK,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(3)) 

     data_headers(4) = 'lapsecent'
     call interp_center(cctkGH,cctk_nghostzones,lapm1,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(4)) 
     export_data(4) = export_data(4) + 1.D0

     data_headers(5) = 'psicent'
     Psi = exp(phi)
     call interp_center(cctkGH,cctk_nghostzones,Psi,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,export_data(5)) 

     data_headers(6) = 'lapsemin'
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(index,"lapse::lapm1")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(6),1,index)
     export_data(6) = export_data(6) + 1.D0

     data_headers(7) = 'phimax'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(index,"bssn::phi")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(7),1,index)

     data_headers(8) = 'rho_b_max'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(index,"mhd_evolve::rho_b")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(8),1,index)

     data_headers(9) = 'rest mass'
     call mass_0_integral(ext, dIntegral, &
          X, Y, Z, &
          rho_star, Symmetry, &
          adjimin, adjimax, &
          adjjmin, adjjmax, &
          adjkmin, adjkmax,CCTK_MyProc(cctkGH))
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(9),CCTK_VARIABLE_REAL)

     data_headers(10) = 'restmassinAH'
     PsiTau = rho_star
     where(abs(emask-1.D0) .lt. 1.D-8) 
        PsiTau = 0.D0
     end where
     call mass_0_integral(ext, dIntegral, &
          X, Y, Z, &
          PsiTau, Symmetry, &
          adjimin, adjimax, &
          adjjmin, adjjmax, &
          adjkmin, adjkmax,CCTK_MyProc(cctkGH))
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(10),CCTK_VARIABLE_REAL)

     data_headers(11) = "AngMom_vol"
     !we use Pr as a temporary variable below:
     call jdens_cpp(cctkGH,cctk_lsh,cctk_nghostzones, Symmetry, &
          Pr,  &
          dx,dy,dz, &
          X,Y,Z,  &
          phi, Sx, Sy, trK, &
          Axx, Axy, Axz, Ayy, Ayz, Azz,  &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
     call mass_0_integral(ext,dIntegral,X,Y,Z,Pr,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(11),CCTK_VARIABLE_REAL)


     data_headers(12) = 'Jsurfpunc1.1M'
     old_surf_radius = surf_radius
     surf_radius=1.1D0*BigM
     call J_surf_integral_punc(CCTK_PASS_FTOF)
     surf_radius = old_surf_radius
     export_data(12) = out_surf_int_sum

     data_headers(13) = "JVol_punc1.1M"
     !we use Pr as a temporary variable below:
     call jdens_cpp(cctkGH,cctk_lsh,cctk_nghostzones, Symmetry, &
          Pr,  &
          dx,dy,dz, &
          X,Y,Z,  &
          phi, Sx, Sy, trK, & 
          Axx, Axy, Axz, Ayy, Ayz, Azz,  &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz)
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,1.1D0*BigM,Xbh1,Ybh1,Zbh1,1.1D0*BigM)
     call mass_0_integral(ext,dIntegral,X,Y,Z,Pr,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(13),CCTK_VARIABLE_REAL)

     data_headers(14) = "M_ADM1.1Surf"
     old_surf_radius = surf_radius
     surf_radius=1.1D0*BigM
     call surf_Mass2_fisheye_punc(CCTK_PASS_FTOF)
     surf_radius = old_surf_radius
     export_data(14) = out_surf_int_sum

     data_headers(15) = "M_ADM1.1exVol"
     call mass_integral_exclude_bh(ext,dIntegral,X,Y,Z, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          Axx,Axy,Axz,Ayy,Ayz,Azz,trK, &
          trRtilde,rho,phi, &
          Gammaxxx,Gammaxxy,Gammaxxz, &
          Gammaxyy,Gammaxyz,Gammaxzz, &
          Gammayxx,Gammayxy,Gammayxz, &
          Gammayyy,Gammayyz,Gammayzz, &
          Gammazxx,Gammazxy,Gammazxz, &
          Gammazyy,Gammazyz,Gammazzz, &
          Symmetry, xbh1, ybh1, zbh1, 1.1D0*BigM, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax);
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(15),CCTK_VARIABLE_REAL)

     data_headers(16) = 'restmassinAH11'
     PsiTau = rho_star
     r_AH = 0.11D0
     where(sqrt((X-Xbh1)*(X-Xbh1) + (Y-Ybh1)*(Y-Ybh1) + (Z-Zbh1)*(Z-Zbh1)) .gt. r_AH) 
        PsiTau = 0.D0
     end where
     call mass_0_integral(ext, dIntegral, &
          X, Y, Z, &
          PsiTau, Symmetry, &
          adjimin, adjimax, &
          adjjmin, adjjmax, &
          adjkmin, adjkmax,CCTK_MyProc(cctkGH))
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(16),CCTK_VARIABLE_REAL)

     data_headers(17) = 'restmassinAH09'
     PsiTau = rho_star
     r_AH = 0.09D0
     where(sqrt((X-Xbh1)*(X-Xbh1) + (Y-Ybh1)*(Y-Ybh1) + (Z-Zbh1)*(Z-Zbh1)) .gt. r_AH) 
        PsiTau = 0.D0
     end where
     call mass_0_integral(ext, dIntegral, &
          X, Y, Z, &
          PsiTau, Symmetry, &
          adjimin, adjimax, &
          adjjmin, adjjmax, &
          adjkmin, adjkmax,CCTK_MyProc(cctkGH))
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(17),CCTK_VARIABLE_REAL)


     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
     !--------------------------------------------------------------------------------!

     filename = 'bhns.bon'

     num_cols = 2
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'UNBOUNDrestmass'
     ! Store -1-u_0 in Pl
     call bhns_compute_neg_one_minus_u_0(ext,Pl, &
          u0,vx,vy,vz, &
          lapm1,shiftx,shifty,shiftz, &
          phi,gxx,gxy,gxz,gyy,gyz,gzz)
     ! Set Pl=0 both inside AH, AND where -1-u_0 < 0
     where(abs(1.D0-emask) .gt. 1.D-8 .or. Pl.lt.0.D0) 
        Pl = 0.D0
     elsewhere
        Pl = 1.D0
     end where
     rho_bl = rho_star * Pl
     call mass_0_integral(ext, dIntegral, &
          X, Y, Z, &
          rho_bl, Symmetry, &
          adjimin, adjimax, &
          adjjmin, adjjmax, &
          adjkmin, adjkmax,CCTK_MyProc(cctkGH))
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,export_data(2),CCTK_VARIABLE_REAL)

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
     !--------------------------------------------------------------------------------!

     filename = 'bhns.kon'
     num_cols = 4
     export_data = 0.D0

     data_headers(1) = 'Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Xbh1'
     export_data(2) = Xbh1

     data_headers(3) = 'Ybh1'
     export_data(3) = Ybh1

     data_headers(4) = 'Zbh1'
     export_data(4) = Zbh1

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !--------------------------------------------------------------------------------!

     do ii=1,nsurf

        surf_radius=surfvec(ii)
        rbr=rbrvec(ii)
        drbdr=drbrvec(ii)
        ddrbddr=ddrbrvec(ii)

        if(nsurf.eq.1) then
           filename = 'bhns.jon'
        else
           write(filename,31)ii
31         FORMAT('bhns.jon.',I1)
        endif

        num_cols = 9

        data_headers(1) = "Time"
        export_data(1) = CCTK_TIME

        data_headers(2) = "fisheye radius"
        export_data(2) = surf_radius

        data_headers(3) = "phys radius"
        export_data(3) = surf_radius*rbr

        data_headers(4) = "Mass_sur"
        call surf_Mass2_fisheye(CCTK_PASS_FTOF)
        export_data(4) = out_surf_int_sum

        data_headers(5) = "Ang_mom_sur"
        call J_surf_integral(CCTK_PASS_FTOF)
        export_data(5) = out_surf_int_sum

        data_headers(6) = "Komar Mass"
        call Komar_Mass(CCTK_PASS_FTOF)
        export_data(6) = out_surf_int_sum

        data_headers(7) = "P_x"
        data_headers(8) = "P_y"
        data_headers(9) = "P_z"
        call P_surf_integral(CCTK_PASS_FTOF)
        export_data(7) = out_surf_px
        export_data(8) = out_surf_py
        export_data(9) = out_surf_pz

        if(CCTK_MyProc(CCTKGH) .eq. 0) then
           call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
        end if
     end do

     !--------------------------------------------------------------------------------!

     filename = 'bhns-ham.con'
     num_cols = 11

     data_headers(1) = "Time"
     export_data(1) = CCTK_TIME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Ham. Constraint: bbh-excluding version.
!!! WARNING: Do not turn on constraint damping (ch2 parameter)
!!!          if this is setting PsiRes!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     data_headers(2) = "Haml ResN"
     !************ HAMILTONIAN CONSTRAINT ************!
     PsiTau = PsiRes
     ! First compute L2 norm of PsiRes.  
     !   Note I don't need to Sync since I don't care
     !        about the ghost zones.
     dIntegral = 0.D0
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
          adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
     export_data(2) = sqrt(resid)

     data_headers(3) = "Haml ResD"
     !Next compute the norm of PsiNorm:
     dIntegral = 0.D0
     call norm2_dv(ext,PsiNorm,X,Y,Z,dIntegral,adjimin, &
          adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry);
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
     export_data(3) = sqrt(resid_norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     data_headers(4) = "Ham Res2P1.1N"
     PsiTau = PsiRes
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,1.1D0*BigM,Xbh1,Ybh1,Zbh1,1.1D0*BigM)
     ! Next the L2 norm.  Note I don't need to Sync since I don't care
     ! about the ghost zones.
     dIntegral = 0.D0
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
          adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
     export_data(4) = sqrt(resid)

     data_headers(5) = "Ham Res2P1.1D"
     PsiTau = PsiNorm
     dIntegral = 0.D0
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,1.1D0*BigM,Xbh1,Ybh1,Zbh1,1.1D0*BigM)
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
          adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry);
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
     export_data(5) = sqrt(resid_norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     data_headers(6) = "Ham Res2P2.0N"
     PsiTau = PsiRes
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,2.0D0*BigM,Xbh1,Ybh1,Zbh1,2.0D0*BigM)
     ! Next the L2 norm.  Note I don't need to Sync since I don't care
     ! about the ghost zones.
     dIntegral = 0.D0
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
          adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
     export_data(6) = sqrt(resid)

     data_headers(7) = "Ham Res2P2.0D"
     PsiTau = PsiNorm
     dIntegral = 0.D0
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,2.0D0*BigM,Xbh1,Ybh1,Zbh1,2.0D0*BigM)
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
          adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry);
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
     export_data(7) = sqrt(resid_norm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     data_headers(8) = "HmRsShExOBN"
     PsiTau = PsiRes
     call vol_integrand_one_spherical_shell(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          0.D0,0.D0,0.D0,0.0D0,0.8*Zglobal(cctk_gsh(3)))
     ! Next the L2 norm.  Note I don't need to Sync since I don't care
     ! about the ghost zones.
     dIntegral = 0.D0
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
          adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
     export_data(8) = sqrt(resid)

     data_headers(9) = "HmRsShExOBD"
     PsiTau = PsiNorm
     dIntegral = 0.D0
     call vol_integrand_one_spherical_shell(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          0.D0,0.D0,0.D0,0.0D0,0.8*Zglobal(cctk_gsh(3)))
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
          adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry);
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
     export_data(9) = sqrt(resid_norm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     data_headers(10) = "HmRsShExOBN"
     PsiTau = PsiRes
     call vol_integrand_one_spherical_shell(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          0.D0,0.D0,0.D0,0.0D0,0.5*Zglobal(cctk_gsh(3)))
     ! Next the L2 norm.  Note I don't need to Sync since I don't care
     ! about the ghost zones.
     dIntegral = 0.D0
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
          adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
     export_data(10) = sqrt(resid)

     data_headers(11) = "HmRsShExOBD"
     PsiTau = PsiNorm
     dIntegral = 0.D0
     call vol_integrand_one_spherical_shell(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          0.D0,0.D0,0.D0,0.0D0,0.5*Zglobal(cctk_gsh(3)))
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
          adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry);
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
     export_data(11) = sqrt(resid_norm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
     !************ END HAMILTONIAN CONSTRAINT ************!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     filename = 'bhns-mom.con'
     num_cols = 17

     data_headers(1) = "Time"
     export_data(1) = CCTK_TIME

     !************ MOMENTUM CONSTRAINT ************!
     data_headers(2) = "Momx ResN"
     data_headers(3) = "Momy ResN"
     data_headers(4) = "Momz ResN"
     data_headers(5) = "Mom ResD"

     call norm2_dv(ext,MRsx,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_x,CCTK_VARIABLE_REAL)
     export_data(2) = sqrt(resid_x)

     call norm2_dv(ext,MRsy,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_y,CCTK_VARIABLE_REAL)
     export_data(3) = sqrt(resid_y)

     call norm2_dv(ext,MRsz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_z,CCTK_VARIABLE_REAL)
     export_data(4) = sqrt(resid_z)

     call norm2_dv(ext,MNorm,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_mom_norm,CCTK_VARIABLE_REAL)
     export_data(5) = sqrt(resid_mom_norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     data_headers(6) = "Momx Res1.1N"
     data_headers(7) = "Momy Res1.1N"
     data_headers(8) = "Momz Res1.1N"
     data_headers(9) = "Mom Res1.1D"

     PsiTau = MRsx
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,1.1D0*BigM,Xbh1,Ybh1,Zbh1,1.1D0*BigM)
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_x,CCTK_VARIABLE_REAL)
     export_data(6) = sqrt(resid_x)

     PsiTau = MRsy
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,1.1D0*BigM,Xbh1,Ybh1,Zbh1,1.1D0*BigM)
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_y,CCTK_VARIABLE_REAL)
     export_data(7) = sqrt(resid_y)

     PsiTau = MRsz
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,1.1D0*BigM,Xbh1,Ybh1,Zbh1,1.1D0*BigM)
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_z,CCTK_VARIABLE_REAL)
     export_data(8) = sqrt(resid_z)

     PsiTau = MNorm
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,1.1D0*BigM,Xbh1,Ybh1,Zbh1,1.1D0*BigM)
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_mom_norm,CCTK_VARIABLE_REAL)
     export_data(9) = sqrt(resid_mom_norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     data_headers(10) = "Momx Res2.0N"
     data_headers(11) = "Momy Res2.0N"
     data_headers(12) = "Momz Res2.0N"
     data_headers(13) = "Mom Res2.0D"

     PsiTau = MRsx
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,2.0D0*BigM,Xbh1,Ybh1,Zbh1,2.0D0*BigM)
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_x,CCTK_VARIABLE_REAL)
     export_data(10) = sqrt(resid_x)

     PsiTau = MRsy
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,2.0D0*BigM,Xbh1,Ybh1,Zbh1,2.0D0*BigM)
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_y,CCTK_VARIABLE_REAL)
     export_data(11) = sqrt(resid_y)

     PsiTau = MRsz
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,2.0D0*BigM,Xbh1,Ybh1,Zbh1,2.0D0*BigM)
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_z,CCTK_VARIABLE_REAL)
     export_data(12) = sqrt(resid_z)

     PsiTau = MNorm
     call vol_integrand_excise_two_punctures(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,2.0D0*BigM,Xbh1,Ybh1,Zbh1,2.0D0*BigM)
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_mom_norm,CCTK_VARIABLE_REAL)
     export_data(13) = sqrt(resid_mom_norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     data_headers(14) = "MxRsShExOBN"
     data_headers(15) = "MyRsShExOBN"
     data_headers(16) = "MzRsShExOBN"
     data_headers(17) = "MomRsShExOBD"

     PsiTau = MRsx
     call vol_integrand_one_spherical_shell(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,0.0D0,0.9*Zglobal(cctk_gsh(3)))
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_x,CCTK_VARIABLE_REAL)
     export_data(14) = sqrt(resid_x)

     PsiTau = MRsy
     call vol_integrand_one_spherical_shell(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,0.0D0,0.9*Zglobal(cctk_gsh(3)))
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_y,CCTK_VARIABLE_REAL)
     export_data(15) = sqrt(resid_y)

     PsiTau = MRsz
     call vol_integrand_one_spherical_shell(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,0.0D0,0.9*Zglobal(cctk_gsh(3)))
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_z,CCTK_VARIABLE_REAL)
     export_data(16) = sqrt(resid_z)

     PsiTau = MNorm
     call vol_integrand_one_spherical_shell(ext,X,Y,Z, &
          PsiTau,Symmetry, &
          Xbh1,Ybh1,Zbh1,0.0D0,0.9*Zglobal(cctk_gsh(3)))
     call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_mom_norm,CCTK_VARIABLE_REAL)
     export_data(17) = sqrt(resid_mom_norm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
     !************ END MOMENTUM CONSTRAINT ************!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !************ GAMMA CONSTRAINT ************!
     filename = 'bhns-gam.con'
     num_cols = 4

     data_headers(1) = "Time"
     export_data(1) = CCTK_TIME

     data_headers(2) = "Gamx Res"
     data_headers(3) = "Gamy Res"
     data_headers(4) = "Gamz Res"

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
     export_data(2) = sqrt(gxres)
     call norm2_dv(ext,gcony,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gyres,CCTK_VARIABLE_REAL)
     export_data(3) = sqrt(gyres)
     call norm2_dv(ext,gconz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
          adjjmax,adjkmin,adjkmax,Symmetry)
     dIntegral = dIntegral*multfactor
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,gzres,CCTK_VARIABLE_REAL)
     export_data(4) = sqrt(gzres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !************ END GAMMA CONSTRAINT ************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do jj=1,moncrief_gw_num_radii

        moncrief_gw_extr_radius = moncrief_radius_GW(jj)        

        ! Gravitational wave extraction!
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

        nmodes_moncrief = ubound(momentsr,1)

        call gw_wave_flux(cctkGH,cctk_nghostzones,moncrief_gw_extr_radius,Mass_GW,cctk_lsh,cctk_gsh, &
             dT,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal, &
             phi,gxx,gxy,gxz,gyy,gyz,gzz,momentsr,momentsi, &
             momentsr_old,momentsi_old,odd_momentsr,odd_momentsi, &
             odd_momentsr_old,odd_momentsi_old,int_momentsr,int_momentsi, &
             hplus,hcross,theta_GW,phi_GW, &
             PhysicalRadius,Radiusderivative)

        if(moncrief_gw_num_radii.eq.1) then
           filename = 'bhns-evenmodes.wave'
        else
           write(filename,41)jj
41         FORMAT('bhns-evenmodes.wave.',I1)
        endif

        num_cols = 2*nmodes_moncrief + 1

        data_headers(1) = "# Time"
        export_data(1) = CCTK_TIME

        do ii=1,nmodes_moncrief
           data_headers(2*ii) = " " 
           data_headers(2*ii+1) = " " 
           export_data(2*ii) = momentsr(ii)
           export_data(2*ii+1) = momentsi(ii)
        end do

        if(CCTK_MyProc(CCTKGH) .eq. 0) then
           call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
        end if
        !--------------------------------------------------------------------------------!

        if(moncrief_gw_num_radii.eq.1) then
           filename = 'bhns-oddmodes.wave'
        else
           write(filename,51)jj
51         FORMAT('bhns-oddmodes.wave.',I1)
        endif
        num_cols = 2*nmodes_moncrief + 1

        data_headers(1) = "# Time"
        export_data(1) = CCTK_TIME

        do ii=1,nmodes_moncrief
           data_headers(2*ii) = " "
           data_headers(2*ii+1) = " "
           export_data(2*ii) = odd_momentsr(ii)
           export_data(2*ii+1) = odd_momentsi(ii)
        end do

        if(CCTK_MyProc(CCTKGH) .eq. 0) then
           call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
        end if

     end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !--------------------------------------------------------------------------------!
     filename = 'bhns-lapmin.mon'
     num_cols = 5
     data_headers(1) = "# Time"
     export_data(1) = CCTK_TIME
     data_headers(2) = "x_coord"
     data_headers(3) = "y_coord"
     data_headers(4) = "z_coord"
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(index,"lapse::lapm1")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(5),1,index)
     data_headers(5) = "lapmin"
     call find_gfvalue_coords(cctkGH,ext,index,lapm1,X,Y,Z,export_data(5),export_data(2),export_data(3),export_data(4))
     if(CCTK_MyProc(cctkGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     filename = 'bhns-rho_bmax.mon'
     num_cols = 5
     data_headers(1) = "# Time"
     export_data(1) = CCTK_TIME
     data_headers(2) = "x_coord"
     data_headers(3) = "y_coord"
     data_headers(4) = "z_coord"
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(index,"mhd_evolve::rho_b")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(5),1,index)
     data_headers(5) = "rho_b_max"
     call find_gfvalue_coords(cctkGH,ext,index,rho_b,X,Y,Z,export_data(5),export_data(2),export_data(3),export_data(4))
     if(CCTK_MyProc(cctkGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     ! Make Courant factor smaller when BH, NS get too close.
     export_data(2) = export_data(2) - Xbh1
     export_data(3) = export_data(3) - Ybh1
     export_data(4) = export_data(4) - Zbh1
     if(sqrt(export_data(2)**2 + export_data(3)**2 + export_data(4)**2) .lt. 0.05D0 &
          .and. CCTK_DELTA_TIME/CCTK_DELTA_SPACE(1) .gt. 0.35D0) then
        cctk_delta_time = 0.875D0 * cctk_delta_time
        write(*,*) "RESETTING COURANT FACTOR!!!!!!!!!"
     end if
     !--------------------------------------------------------------------------------!

  end if

  !Need to reset mask.  Otherwise AH finder will basically turn off.
  emask = 1.D0

end subroutine bhns_diagnostics

subroutine bhns_compute_neg_one_minus_u_0(ext,neg_one_minus_u_0, &
     u0,vx,vy,vz, &
     lapse,shiftx,shifty,shiftz, &
     phi,gxx,gxy,gxz,gyy,gyz,gzz)
  implicit none
  integer, dimension(3)                          :: ext
  real*8, dimension(ext(1),ext(2),ext(3))        :: u0,vx,vy,vz
  real*8, dimension(ext(1),ext(2),ext(3))        :: lapse,shiftx,shifty,shiftz
  real*8, dimension(ext(1),ext(2),ext(3))        :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))        :: neg_one_minus_u_0

  ! u_0 = u^0 * (-alpha^2 + gamma_ij beta^i (beta^j + v^j))
  neg_one_minus_u_0 = -1.D0 - u0*(-(lapse+1.D0)**2 + &
       exp(4.D0*phi) * ( &
       gxx*shiftx*(shiftx+vx) + &
       2.D0*gxy*shiftx*(shifty+vy) + &
       2.D0*gxz*shiftx*(shiftz+vz) + &
       gyy*shifty*(shifty+vy) + &
       2.D0*gyz*shifty*(shiftz+vz) + &
       gzz*shiftz*(shiftz+vz)))

end subroutine bhns_compute_neg_one_minus_u_0
