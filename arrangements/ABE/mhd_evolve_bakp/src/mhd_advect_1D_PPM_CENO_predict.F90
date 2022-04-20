
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!-----------------------------------------------------------------------------
! PPM, PPM+, and CENO Advection routines for hydro variables
!-----------------------------------------------------------------------------
subroutine mhd_advect_1D_PPM_CENO_predict(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS: 
  integer                                  :: i,j,k
  real*8                                   :: dX,dY,dZ
  integer, dimension(3)                    :: ext,U_syms1,U_syms2,U_syms3,U_syms4,U_syms5,U_syms6,U_syms7,U_syms8
  real*8                                   :: n_th

  integer :: dummy,handle
  integer :: index
  integer :: ierr
  CCTK_REAL reduction_value
  integer :: PPM_PLUS, PPM, CENO,MC,SPPM
  parameter(PPM_PLUS = 1, PPM = 2,CENO = 3,MC = 4, SPPM = 5)
  integer :: SCALAR,VECTORX,VECTORY,VECTORZ,BVECTORX,BVECTORY,BVECTORZ
  parameter(SCALAR = 1,VECTORX = 2,VECTORY = 3,VECTORZ = 4,BVECTORX = 5,BVECTORY = 6,BVECTORZ = 7)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(Reconstruction==PPM .or. Reconstruction==PPM_PLUS) then
     !=================================================
     !  Advance to PPM fourth-order reconstruction on faces!
     !=================================================
     !compute first and second centered derivatives
     ! modify to get rid of unneeded derivatives
     call find_centderivs(ext,X,Y,Z,rho_b,d0rho_b_m,d02rho_b_m, &
          P,d0P_m,d02P_m, &
          vx,d0vx_m,d02vx_m, &
          vy,d0vy_m,d02vy_m, &
          vz,d0vz_m,d02vz_m, &
          m,Symmetry)
     if(excision_enable == 1) then
	call scalar_excision_bc(ext,X,Y,Z, &
             d0rho_b_m,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             d02rho_b_m,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             d0P_m,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             d02P_m,Symmetry,excision_zone_gf)
	call vector_excision_bc(ext,X,Y,Z, &
             d0vx_m,d0vy_m,d0vz_m, &
             Symmetry,excision_zone_gf)
	call vector_excision_bc(ext,X,Y,Z, &
             d02vx_m,d02vy_m,d02vz_m, &
             Symmetry,excision_zone_gf)
     end if

     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_d0_quantities')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_d0_quantities')

     call hydro_forwardback(ext,X,Y,Z,d02rho_b_m,d2rho_b_bck,d2rho_b_fwd, &
          d02P_m,d2P_bck,d2P_fwd, &
          d02vx_m,d2vx_bck,d2vx_fwd, &
          d02vy_m,d2vy_bck,d2vy_fwd, &
          d02vz_m,d2vz_bck,d2vz_fwd, &
          m,Symmetry)

     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_hydro_fwdbck_quantities')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_hydro_fwdbck_quantities')

     ! find a_{j+1/2}
     if(adm_ppm_b .eq. 1) then 
        call mhdppm_find_face_vals(ext, X,Y,Z, &
             rho_b,drho_b_m,rho_br_ppm,rho_bl_ppm, &
             P,dP_m,Pr,Pl, &
             vx,dvx_m,vxr,vxl, &
             vy,dvy_m,vyr,vyl, &
             vz,dvz_m,vzr,vzl, &
             Bx,dBx_m,Bxr,Bxl, &
             By,dBy_m,Byr,Byl, &
             Bz,dBz_m,Bzr,Bzl, &
             m,Symmetry,Sym_Bz)
     else
        call ppm_find_face_vals(ext,X,Y,Z, &
             rho_b,drho_b_m,rho_br_ppm,rho_bl_ppm, &
             P,dP_m,Pr,Pl, &
             vx,dvx_m,vxr,vxl, &
             vy,dvy_m,vyr,vyl, &
             vz,dvz_m,vzr,vzl, &
             m,Symmetry)
     end if

     ! steepen rho
     n_th = 1.0/(gamma_th-1.0)
     call ppm_steepen_rho(ext,X,Y,Z,rho_b,d0rho_b_m,d02rho_b_m, &
          rho_br_ppm,rho_bl_ppm,rho_br,rho_bl, &
          P,n_th,rho_b_max,m,Symmetry)    
     if(excision_enable == 1) then
	call scalar_excision_bc(ext,X,Y,Z, &
             rho_br_ppm,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             rho_bl_ppm,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             Pr,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             Pl,Symmetry,excision_zone_gf)
	call vector_excision_bc(ext,X,Y,Z, &
             Bxr,Byr,Bzr, &
             Symmetry,excision_zone_gf)
	call vector_excision_bc(ext,X,Y,Z, &
             Bxl,Byl,Bzl, &
             Symmetry,excision_zone_gf)
	call vector_excision_bc(ext,X,Y,Z, &
             vxr,vyr,vzr, &
             Symmetry,excision_zone_gf)
	call vector_excision_bc(ext,X,Y,Z, &
             vxl,vyl,vzl, &
             Symmetry,excision_zone_gf)
     end if

     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')

     ! monotonize
     call ppm_monotonize(ext,X,Y,Z,rho_b_max, &
          rho_b,rho_br_ppm,rho_bl_ppm, &
          d02rho_b_m,d2rho_b_bck,d2rho_b_fwd, &
          P,Pr,Pl, &
          d02P_m,d2P_bck,d2P_fwd, &
          vx,vxr,vxl, &
          d02vx_m,d2vx_bck,d2vx_fwd, &
          vy,vyr,vyl, &
          d02vy_m,d2vy_bck,d2vy_fwd, &
          vz,vzr,vzl, &
          d02vz_m,d2vz_bck,d2vz_fwd, &
          m,Symmetry,Reconstruction)
     if(adm_ppm_b .eq. 1) then
        call mhdppm_monotonize(ext,X,Y,Z,Bx,Bxr,Bxl,By,Byr,Byl, &
             Bz,Bzr,Bzl,m,Symmetry,Sym_Bz)    
     end if
     ! find ftilde
     call ppm_ftilde(ext,X,Y,Z,dP_m,P,d0P_m,vx,vy,vz, &
          P_max,m,Symmetry)

     if(excision_enable == 1) then
        call scalar_excision_bc(ext,X,Y,Z, &
             rho_br_ppm,Symmetry,excision_zone_gf)
        call scalar_excision_bc(ext,X,Y,Z, &
             rho_bl_ppm,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             Pr,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             Pl,Symmetry,excision_zone_gf)
        if(adm_ppm_b .eq. 1) then
           call vector_excision_bc(ext,X,Y,Z, &
                Bxr,Byr,Bzr, &
                Symmetry,excision_zone_gf)
           call vector_excision_bc(ext,X,Y,Z, &
                Bxl,Byl,Bzl, &
                Symmetry,excision_zone_gf)
        end if
	call vector_excision_bc(ext,X,Y,Z, &
             vxr,vyr,vzr, &
             Symmetry,excision_zone_gf)
	call vector_excision_bc(ext,X,Y,Z, &
             vxl,vyl,vzl, &
             Symmetry,excision_zone_gf)
     end if

     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')

     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_nablas')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_nablas')
     if(adm_ppm_b .eq. 1) then 
        call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')
        call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')
     end if
     ! redefine ar, al
     if(adm_ppm_b .eq. 1) then 
        call mhdppm_shift(ext,X,Y,Z,rho_br_ppm,rho_bl_ppm, &
             Pr,Pl, &
             vxr,vxl,vyr,vyl,vzr,vzl, &
             Bxr,Bxl,Byr,Byl,Bzr,Bzl, &
             m,Symmetry)
     else
        call ppm_shift(ext,X,Y,Z,rho_br_ppm,rho_bl_ppm, &
             Pr,Pl,vxr,vxl, &
             vyr,vyl,vzr,vzl,m,Symmetry)       
     end if

     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
     if(adm_ppm_b .eq. 1) then 
        call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')
        call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')
     end if

     ! flatten
     if(adm_ppm_b .eq. 1) then 
        call mhdppm_flatten(ext,X,Y,Z,dP_m, &
             rho_b,rho_br_ppm,rho_bl_ppm, &
             P,Pr,Pl, &
             vx,vxr,vxl, &
             vy,vyr,vyl, &
             vz,vzr,vzl, &
             Bx,Bxr,Bxl, &
             By,Byr,Byl, &
             Bz,Bzr,Bzl, &
             m,Symmetry)
     else 
        call ppm_flatten(ext,X,Y,Z,dP_m, &
             rho_b,rho_br_ppm,rho_bl_ppm, &
             P,Pr,Pl, &
             vx,vxr,vxl, &
             vy,vyr,vyl, &
             vz,vzr,vzl, &
             m,Symmetry)  
     end if
     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')

     rho_br=rho_br_ppm
     rho_bl=rho_bl_ppm
  else if(Reconstruction.eq.CENO) then 
     !=================================================
     !  Advance to CENO third-order reconstruction on faces!
     !=================================================
     !compute first and second centered derivatives
     call find_centderivs(ext,X,Y,Z,rho_b,d0rho_b_m,d02rho_b_m, &
          P,d0P_m,d02P_m, &
          vx,d0vx_m,d02vx_m, &
          vy,d0vy_m,d02vy_m, &
          vz,d0vz_m,d02vz_m, &
          m,Symmetry)        
     if(excision_enable == 1) then
	call scalar_excision_bc(ext,X,Y,Z, &
             d0rho_b_m,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             d02rho_b_m,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             d0P_m,Symmetry,excision_zone_gf)
	call scalar_excision_bc(ext,X,Y,Z, &
             d02P_m,Symmetry,excision_zone_gf)
	call vector_excision_bc(ext,X,Y,Z, &
             d0vx_m,d0vy_m,d0vz_m, &
             Symmetry,excision_zone_gf)
	call vector_excision_bc(ext,X,Y,Z, &
             d02vx_m,d02vy_m,d02vz_m, &
             Symmetry,excision_zone_gf)
     end if

     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_d0_quantities')  
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_d0_quantities')
     ! Compute the third-order left and right face states (P_L and P_R)
     call find_ur_ul_ceno(ext, &
          rho_b,d0rho_b_m, &
          d02rho_b_m,rho_br,rho_bl, &
          P,d0P_m,d02P_m,Pr,Pl, &
          vx,d0vx_m,d02vx_m,vxr,vxl, &
          vy,d0vy_m,d02vy_m,vyr,vyl, &
          vz,d0vz_m,d02vz_m,vzr,vzl, &
          m,Symmetry,lapm1_f, &
          shiftx_f,shifty_f,shiftz_f,phi_f, &
          gxx_f,gxy_f,gxz_f, &
          gyy_f,gyz_f,gzz_f)
     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')
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


end subroutine mhd_advect_1D_PPM_CENO_predict
