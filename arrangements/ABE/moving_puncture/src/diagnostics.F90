#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------
! Moving puncture diagnostics
!-----------------------------
subroutine movpunc_diagnostics(CCTK_ARGUMENTS)

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

  if(mod(cctk_iteration,diags_every).eq.0) then
  
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

  !  call setup_global_coord_arrays(CCTK_PASS_FTOF)

  !header_flag == 1 -> write file headers
  if(CCTK_TIME == 0.D0) then
     header_flag = 1
  else
     header_flag = 0
  end if

  do ii=1,nsurf

     surf_radius=surfvec(ii)
     rbr=rbrvec(ii)
     drbdr=drbrvec(ii)
     ddrbddr=ddrbrvec(ii)

     !--------------------------------------------------------------------------------!
     if(nsurf.eq.1) then
        filename = 'movpunc.jon'
     else
        write(filename,31)ii
31      FORMAT('movpunc.jon.',I1)
     endif

     num_cols = 6

     data_headers(1) = "Time"
     export_data(1) = CCTK_TIME

     data_headers(2) = "Mass_sur"
     call surf_Mass2_fisheye(CCTK_PASS_FTOF)
     export_data(2) = out_surf_int_sum

     data_headers(3) = "Komar Mass"
     call Komar_Mass(CCTK_PASS_FTOF)
     export_data(3) = out_surf_int_sum

     data_headers(4) = "P_x"
     data_headers(5) = "P_y"
     data_headers(6) = "P_z"
     call P_surf_integral(CCTK_PASS_FTOF)
     export_data(4) = out_surf_px
     export_data(5) = out_surf_py
     export_data(6) = out_surf_pz

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
  enddo
  !--------------------------------------------------------------------------------!

  filename = 'movpunc.con'
  num_cols = 14

  data_headers(1) = "Time"
  export_data(1) = CCTK_TIME

  data_headers(2) = "Ham_Res"
  data_headers(3) = "Ham_norm"
  data_headers(4) = "Ham_Resp"
  data_headers(5) = "Ham_normp"
  !************ HAMILTONIAN CONSTRAINT ************!
  !  call pre_con(ext,phi,Psi,KK, &
       !       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
  !       Axx,Axy,Axz,Ayy,Ayz,Azz,trK)
  !  call ham_residual(ext,X,Y,Z, &
       !       Psi,PsiRes,PsiTau,KK,trRtilde,rho, &
  !       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       !       Gammax,Gammay,Gammaz,Symmetry)
  !imin,jmin,kmin,sNx,sNz,PI_SYMM_Psi_gz)
  !  PsiTau = PsiRes
  !
  ! Next the L2 norm.  Note I don't need to Sync since I don't care
  ! about the ghost zones.
  !
  PsiTau = PsiRes
  dIntegral = 0.D0
  call norm2_dv(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
  resid = sqrt(resid)
  export_data(2) = resid

  !  call ham_norm(ext,X,Y,Z, &
       !       Psi,PsiNorm,PsiTau,KK,trRtilde,rho, &
  !       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       !       Gammax,Gammay,Gammaz,Symmetry) !, &
  !imin,jmin,kmin,sNx,sNz,PI_SYMM_Psi_gz)
  dIntegral = 0.D0
  call norm2_dv(ext,PsiNorm,X,Y,Z,dIntegral,adjimin, &
       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
  resid_norm = sqrt(resid_norm)
  export_data(3) = resid_norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! the puncture excluding version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Next the L2 norm.  Note I don't need to Sync since I don't care
  ! about the ghost zones.
  !
  PsiTau = PsiRes
  dIntegral = 0.D0
  call norm2_dv_punc(ext,PsiTau,X,Y,Z,dIntegral,adjimin, &
       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry,xbh1,ybh1,zbh1)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid,CCTK_VARIABLE_REAL)
  resid = sqrt(resid)
  export_data(4) =resid

  dIntegral = 0.D0
  call norm2_dv_punc(ext,PsiNorm,X,Y,Z,dIntegral,adjimin, &
       adjimax,adjjmin,adjjmax,adjkmin,adjkmax,Symmetry,xbh1,ybh1,zbh1)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_norm,CCTK_VARIABLE_REAL)
  resid_norm = sqrt(resid_norm)
  export_data(5) =resid_norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(CCTK_PASS_FTOF)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !write(*,*) "Resid = ",resid,resid_norm,resid/resid_norm
  !************ END HAMILTONIAN CONSTRAINT ************!

  data_headers(6) = "Momx Res"
  data_headers(7) = "Momy Res"
  data_headers(8) = "Momz Res"
  data_headers(9) = "Momx Resp"
  data_headers(10) = "Momy Resp"
  data_headers(11) = "Momz Resp"
  !************ MOMENTUM CONSTRAINT ************!
  call norm2_dv(ext,MRsx,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_x,CCTK_VARIABLE_REAL)
  resid_x = sqrt(resid_x)

  call norm2_dv(ext,MRsy,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_y,CCTK_VARIABLE_REAL)
  resid_y = sqrt(resid_y)

  call norm2_dv(ext,MRsz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_z,CCTK_VARIABLE_REAL)
  resid_z = sqrt(resid_z)

  call norm2_dv(ext,MNorm,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_mom_norm,CCTK_VARIABLE_REAL)
  resid_mom_norm = sqrt(resid_mom_norm)
  export_data(6) = resid_x/(resid_mom_norm+1.0d-10)
  export_data(7) = resid_y/(resid_mom_norm+1.0d-10)
  export_data(8) = resid_z/(resid_mom_norm+1.0d-10)

  !************ END MOMENTUM CONSTRAINT ************!
  !************ MOMENTUM CONSTRAINT w/excluded puncture************!
  call norm2_dv_punc(ext,MRsx,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_x,CCTK_VARIABLE_REAL)
  resid_x = sqrt(resid_x)

  call norm2_dv_punc(ext,MRsy,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_y,CCTK_VARIABLE_REAL)
  resid_y = sqrt(resid_y)

  call norm2_dv_punc(ext,MRsz,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_z,CCTK_VARIABLE_REAL)
  resid_z = sqrt(resid_z)

  call norm2_dv_punc(ext,MNorm,X,Y,Z,dIntegral,adjimin,adjimax,adjjmin, &
       adjjmax,adjkmin,adjkmax,Symmetry)
  dIntegral = dIntegral*multfactor
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,resid_mom_norm,CCTK_VARIABLE_REAL)
  resid_mom_norm = sqrt(resid_mom_norm)

  export_data(9)  = resid_x/(resid_mom_norm+1.0d-10)
  export_data(10) = resid_y/(resid_mom_norm+1.0d-10)
  export_data(11) = resid_z/(resid_mom_norm+1.0d-10)

  data_headers(12) = "xbh"
  data_headers(13) = "ybh"
  data_headers(14) = "zbh"
  export_data(12) = xbh1
  export_data(13) = ybh1
  export_data(14) = zbh1

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !--------------------------------------------------------------------------------!  
endif

end subroutine movpunc_diagnostics
