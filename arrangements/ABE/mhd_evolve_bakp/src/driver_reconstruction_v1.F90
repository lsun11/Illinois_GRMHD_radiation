!--------------------------------------------------------------------------
! Reconstruction driver (slope limiter, metric & primitve facevals) : v1.0
!--------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine reconstruction_v1(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS: 
  real*8                :: dX,dY,dZ
  integer, dimension(3) :: ext

  CCTK_REAL :: reduction_value
  integer   :: dummy,handle,index,ierr
  integer   :: PPM_PLUS, PPM, CENO,MC,SPPM
  parameter(PPM_PLUS = 1, PPM = 2,CENO = 3,MC = 4, SPPM = 5)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ! compute \nabla P
  call compute_nabla_hydro(ext,dX,dY,dZ,vx,dvx_m,vy,dvy_m, &
       vz,dvz_m,rho_b,drho_b_m, &
       P,dP_m,m,Symmetry,X,Y,Z,Reconstruction)
  call compute_nabla_emfields(ext,dX,dY,dZ,Bx,By,Bz,dBx_m, &
       dBy_m,dBz_m,m,Symmetry,Sym_Bz,X,Y,Z, &
       Reconstruction)

  if(excision_enable == 1) then
     call scalar_excision_bc(ext,X,Y,Z, &
          drho_b_m,Symmetry,excision_zone_gf)
     call scalar_excision_bc(ext,X,Y,Z, &
          dP_m,Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z, &
          dvx_m,dvy_m,dvz_m, &
          Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z, &
          dBx_m,dBy_m,dBz_m, &
          Symmetry,excision_zone_gf)
  end if

  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_nablas')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_nablas')

  ! Compute the left and right states (P_L and P_R)
  call find_ur_ul_hydro(ext,X,Y,Z,vx,vxr,vxl,dvx_m, &
       vy,vyr,vyl,dvy_m, &
       vz,vzr,vzl,dvz_m, &
       rho_b,rho_br,rho_bl,drho_b_m, &
       P,Pr,Pl,dP_m, &
       m,Symmetry,Reconstruction)
  
  call find_ur_ul_emfields(ext,X,Y,Z,Bx,By,Bz,Bxr,Bxl, &
       Byr,Byl,Bzr,Bzl,dBx_m, &
       dBy_m,dBz_m,m,Symmetry,Sym_Bz, &
       Reconstruction)

!  write(*,*) "HELLO AFTER COMPUTE UR UL BZL = ",Bzl(2,2,2),Bzr(2,2,2),dBz_m(1,2,2),dBz_m(2,2,2),Bz(1,2,2),Bz(2,2,2)

  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')

  ! Average to find metric quantities at the interfaces  
  call face_avg_hydro(ext,X,Y,Z,lapm1,lapm1_f,shiftx,shiftx_f, &
       shifty,shifty_f,shiftz,shiftz_f, &
       gxx,gxx_f,gxy,gxy_f,gxz,gxz_f, &
       gyy,gyy_f,gyz,gyz_f,gzz,gzz_f, &
       phi,phi_f,gupxx,gupxx_f,gupyy, &
       gupyy_f,gupzz,gupzz_f,m,Symmetry)

  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_metric_facevals')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_metric_facevals')

  if(excision_enable == 1) then
     call scalar_excision_bc(ext,X,Y,Z, &
          lapm1_f,Symmetry,excision_zone_gf)
     call scalar_excision_bc(ext,X,Y,Z, &
          phi_f,Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z, &
          shiftx_f,shifty_f,shiftz_f, &
          Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z, &
          gxx_f,gxy_f,gxz_f, &
          Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z, &
          gyy_f,gyz_f,gzz_f, &
          Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z, &
          gupxx_f,gupyy_f,gupzz_f, &
          Symmetry,excision_zone_gf)
  end if

  if(Reconstruction.ne.MC) then
     write(*,*) "SORRY, CENO OR PPM RECONSTRUCTION NOT SUPPORTED WITH RECONSTRUCTION V1.  Please try reconstruction_v2 for PPM, or port over CENO!"
     stop
     !call mhd_advect_1D_PPM_CENO_predict(CCTK_PASS_FTOF)
  end if

  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')

  if(Reconstruction==SPPM) then
     call sppm_shift_hydro_ul(ext,Z,rho_br,rho_bl,Pr,Pl, &
          vxr,vxl,vyr,vyl,vzr,vzl,m,Symmetry)
     call sppm_shift_emfields_ul(ext,Z,Bxr,Bxl,Byr,Byl, &
          Bzr,Bzl,m,Symmetry,Sym_Bz)
     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')
  end if

end subroutine reconstruction_v1
