!--------------------------------------------------
! Compute B^i fluxes from induction equation: v1.5
!--------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine flux_induction_1D_v2(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS: 
  real*8                :: dX,dY,dZ
  integer               :: pow_axi,i,j,k
  integer, dimension(3) :: ext

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(reconstruct_Bitildes_instead_of_Bis==0) then
     call convert_Bil_Bir_to_Btildeil_Btildeir(cctkGH,ext,Bxl,Bxr,Byl,Byr,Bzl,Bzr,phi_f)
  end if

  if(reconstruct_Bitildes_instead_of_Bis==1) then
     call compute_cmax_cmin_hybrid_cpp(cctkGH,ext,v02r,v02l,cmax,cmin, &
          rho_br,rho_bl, &
          Pr,Pl, vxr,vxl,vyr,vyl,vzr,vzl, &
          Bxr,Bxl,Byr,Byl,Bzr,Bzl,lapm1_f, &
          shiftx_f,shifty_f,shiftz_f, phi_f, &
          gxx_f,gxy_f,gxz_f,gyy_f,gyz_f,gzz_f, &
          gupxx_f,gupyy_f,gupzz_f,m, &
          neos,ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
          temp1,temp2,temp3,temp4,temp5)
  end if

  if(m==1) then
     !  fxy = vx*By - vy*Bx,   fxz = vx*Bz - vz*Bx
     if(enable_disk_em_flux_induction==0) then
        pow_axi = 0
        call flux_induction_cpp(cctkGH,ext, X,fxy,vxr,vyr,Bxr,Byr, &
             vxl,vyl,Bxl,Byl,cmax,cmin, &
             pow_axi,Symmetry)
        pow_axi = 1
        call flux_induction_cpp(cctkGH,ext, X,fxz,vxr,vzr,Bxr,Bzr, &
             vxl,vzl,Bxl,Bzl,cmax,cmin, &
             pow_axi,Symmetry)
     else if(enable_disk_em_flux_induction==1) then
!!$        pow_axi = 0
!!$        call flux_induction_disk(ext, X,fxy,vxr,vyr,Bxr,Byr, &
!!$             vxl,vyl,Bxl,Byl,cmax,cmin, &
!!$             pow_axi,Symmetry,m)
!!$        pow_axi = 1
!!$        call flux_induction_disk(ext, X,fxz,vxr,vzr,Bxr,Bzr, &
!!$             vxl,vzl,Bxl,Bzl,cmax,cmin, &
!!$             pow_axi,Symmetry,m)
     end if

     if(excision_enable == 1) then
        call scalar_excision_bc(ext,X,Y,Z,fxy,Symmetry,excision_zone_gf);
        call scalar_excision_bc(ext,X,Y,Z,fxz,Symmetry,excision_zone_gf);
     end if

  else if(m==2) then
     !  fyx = vy*Bx - vx*By,   fyz = vy*Bz - vz*By
     ! No need to call flux_induction_disk here, even when you are
     ! evolving a disk. Zach says: I believe this is since we do disks in axisymmetry, so m!=2 always.
     pow_axi = 0
     call flux_induction_cpp(cctkGH,ext,X, fyx,vyr,vxr,Byr,Bxr, &
          vyl,vxl,Byl,Bxl,cmax,cmin, &
          pow_axi,Symmetry)
     call flux_induction_cpp(cctkGH,ext,X, fyz,vyr,vzr,Byr,Bzr, &
          vyl,vzl,Byl,Bzl,cmax,cmin, &
          pow_axi,Symmetry)                     
     if(excision_enable == 1) then
        call scalar_excision_bc(ext,X,Y,Z,fyx,Symmetry,excision_zone_gf);
        call scalar_excision_bc(ext,X,Y,Z,fyz,Symmetry,excision_zone_gf);
     end if

  else if(m==3) then     
     !  fzx = vz*Bx - vx*Bz,   fzy = vz*By - vy*Bz

     if (enable_disk_em_flux_induction==0) then
        ! next line contains bugfix by Branson for axisymmetric case.
        pow_axi = 2
        call flux_induction_cpp(cctkGH,ext, X,fzx,vzr,vxr,Bzr,Bxr, &
             vzl,vxl,Bzl,Bxl,cmax,cmin, &
             pow_axi,Symmetry)
        pow_axi = 0
        call flux_induction_cpp(cctkGH,ext, X,fzy,vzr,vyr,Bzr,Byr, &
             vzl,vyl,Bzl,Byl,cmax,cmin, &
             pow_axi,Symmetry)
     else if(enable_disk_em_flux_induction==1) then
!!$        pow_axi = 2
!!$        call flux_induction_disk(ext, X,fzx,vzr,vxr,Bzr,Bxr, &
!!$             vzl,vxl,Bzl,Bxl,cmax,cmin, &
!!$             pow_axi,Symmetry,m)
!!$        pow_axi = 0
!!$        call flux_induction_disk(ext, X,fzy,vzr,vyr,Bzr,Byr, &
!!$             vzl,vyl,Bzl,Byl,cmax,cmin, &
!!$             pow_axi,Symmetry,m)
     end if

     if(excision_enable == 1) then
        call scalar_excision_bc(ext,X,Y,Z,fzx,Symmetry,excision_zone_gf);
        call scalar_excision_bc(ext,X,Y,Z,fzy,Symmetry,excision_zone_gf);
     end if
  end if

end subroutine flux_induction_1D_v2
