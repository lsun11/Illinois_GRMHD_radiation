#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for bondi thorn
!-----------------------------------------------------------------------------
subroutine statpunc_bondi_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,Zmin
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20
  integer                                  :: vindex,found_horizon,horizon_number,foundflag
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  !Needed for GW extraction:
  real*8                                   :: Mass_GW,hplus,hcross,moncrief_gw_extr_radius
  real*8                                   :: psi4r,psi4i
  integer                                  :: jj,nmodes_moncrief
  real*8,dimension(2)                      :: gw_amplitudes
  !End GW extraction variables  

  real*8                                   :: old_surf_radius,SymmFactor
  real*8                                   :: rho_tiny,F_M0,vmag_surf,surf_area,ave_radius
  real*8                                   :: F_M0_BH
  real*8                                   :: umag_surf,u0_surf
  real*8                                   :: rho_b_surf,P_surf,surf
  real*8                                   :: out_surf_int_sum,out_surf_px,out_surf_py,out_surf_pz
  real*8, dimension(1,3)                   :: pointcoords
  CCTK_POINTER :: reduction_value_pointer
  integer :: header_flag,handle
  integer :: num_cols
  integer :: ierr,myproc_rank,ii
  real*8,dimension(4)                        :: drho_dtau_surf
  real *8                                  :: rho_pos,rho_old_pos,drhodt_pos

!if(1==0) then

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(MOD(cctk_iteration,out_every)==0) then

     ! Check parameter ntot_bondi
     if (cctk_iteration==0) then 
	if (ntot_bondi .ne. N_theta*N_phi) then 
	   write(*,*) 'In bondi_amr, ntot_bondi = ',ntot_bondi
	   write(*,*) 'N_theta, N_phi = ',N_theta, N_phi
	   write(*,*) 'ntot_bondi must be set to N_theta*N_phi'
	   write(*,*) 'Need to reset the parameter. Bye now...'
	   stop
	end if
     end if


     call set_SymmFactor(Symmetry,SymmFactor)

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1
     else
        header_flag = 0
     end if

     !--------------------------------------------------------------------------------!
     filename = 'bondi.mon'
     num_cols = 13
     export_data = 0.D0
     data_headers = 'unknowndude'

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Phicent'
     call CCTK_VarIndex (vindex, "bssn::phi")
     call interp_center_carp(cctkGH,vindex,export_data(2))

     data_headers(3) = 'Kcent'
     call CCTK_VarIndex (vindex, "bssn::trK")
     call interp_center_carp(cctkGH,vindex,export_data(3))

     data_headers(4) = 'lapsecent'
     call CCTK_VarIndex (vindex, "lapse::lapm1")
     call interp_center_carp(cctkGH,vindex,export_data(4))
     export_data(4) = export_data(4) + 1.D0

     data_headers(5) = 'psicent'
     !     Psi = exp(phi)
     call CCTK_VarIndex (vindex, "bssn::psi")
     call interp_center_carp(cctkGH,vindex,export_data(5))

     data_headers(6) = 'lapsemin'
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(vindex,"lapse::lapm1")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(6),1,vindex)
     export_data(6) = export_data(6) + 1.D0

     data_headers(7) = 'phimax'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"bssn::phi")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(7),1,vindex)

!!$     data_headers(8) = 'rhocent'
!!$     call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
!!$     call interp_center_carp(cctkGH,vindex,export_data(8))
!!$
!!$     data_headers(9) = 'rho_b_max'
!!$     call CCTK_ReductionHandle(handle,"maximum")
!!$     call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
!!$     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(9),1,vindex)


     data_headers(8) = 'Jsurfpunc1.1M'
     old_surf_radius = surf_radius
     surf_radius=1.1D0*bh_mass
     call J_surf_integral_offcenter(cctkGH,0.D0,0.D0,0.D0,export_data(8))
     surf_radius = old_surf_radius

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(9) = "J_1.1exVol"
     export_data(9) = J_ADM_VolInt*SymmFactor

     data_headers(10) = "J_Surf_+_Vol"
     export_data(10) = J_ADM_VolInt*SymmFactor + export_data(8)

     data_headers(11) = "M_ADMsurf_BH_1.1M"
     old_surf_radius = surf_radius
     surf_radius=1.1D0*bh_mass
     call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(11))
     surf_radius = old_surf_radius

     ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(12) = "M_ADM1.1exVol"
     export_data(12) = M_ADM_VolInt*SymmFactor

     data_headers(13) = "M_Surf_+_Vol"
     export_data(13) = M_ADM_VolInt*SymmFactor + export_data(11)


     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     do ii=1,nsurf

        surf_radius=surfvec(ii)
        rbr=rbrvec(ii)
        drbdr=drbrvec(ii)
        ddrbddr=ddrbrvec(ii)

        if(nsurf.eq.1) then
           filename = 'bondi.jon'
        else
           write(filename,21)ii
21         FORMAT("bondi.jon.",I1)
        endif

        num_cols = 13

        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME

        data_headers(2) = "fisheye radius"
        export_data(2) = surf_radius

        data_headers(3) = "phys radius"
        export_data(3) = surf_radius*rbr

        data_headers(4) = "Mass_sur"
        call M_surf_integral(cctkGH,export_data(4))

        data_headers(5) = "Ang_mom_surf"
        call J_surf_integral_offcenter(cctkGH,0.D0,0.D0,0.D0,export_data(5))

        data_headers(6) = "Komar Mass"
        call Komar_surf_integralv2(cctkGH,export_data(6))

        data_headers(7) = "P_x2"
        data_headers(8) = "P_y2"
        data_headers(9) = "P_z2"
        call P_surf_integralv2(cctkGH,export_data(7),export_data(8),export_data(9))
	data_headers(10) = "M_surf_mov"
	call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(10))

        data_headers(11) = "P_x"
        data_headers(12) = "P_y"
        data_headers(13) = "P_z"
        call P_surf_integral_offcenter(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(11),export_data(12),export_data(13))


        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end do

!!$     !--------------------------------------------------------------------------------!
!!$     filename = 'bondi.bon'
!!$     num_cols = 13
!!$     export_data = 0.D0
!!$
!!$     data_headers(1) = '# Time'
!!$     export_data(1) = CCTK_TIME
!!$
!!$     call find_drho_dtau(cctkGH,bondi_rhosurf1,drho_dtau_surf,0,N_rad_ray)
!!$     data_headers(2) = 'drho_dt 1 +z'
!!$     export_data(2) = drho_dtau_surf(1)
!!$     data_headers(3) = 'drho_dt 1 +x'
!!$     export_data(3) = drho_dtau_surf(2)
!!$     data_headers(4) = 'drho_dt 1 -x'
!!$     export_data(4) = drho_dtau_surf(3)
!!$     call find_drho_dtau(cctkGH,bondi_rhosurf2,drho_dtau_surf,0,N_rad_ray)
!!$     data_headers(5) = 'drho_dt 2 +z'
!!$     export_data(5) = drho_dtau_surf(1)
!!$     data_headers(6) = 'drho_dt 2 +x'
!!$     export_data(6) = drho_dtau_surf(2)
!!$     data_headers(7) = 'drho_dt 2 -x'
!!$     export_data(7) = drho_dtau_surf(3)
!!$     call find_drho_dtau(cctkGH,bondi_rhosurf3,drho_dtau_surf,0,N_rad_ray)
!!$     data_headers(8) = 'drho_dt 3 +z'
!!$     export_data(8) = drho_dtau_surf(1)
!!$     data_headers(9) = 'drho_dt 3 +x'
!!$     export_data(9) = drho_dtau_surf(2)
!!$     data_headers(10) = 'drho_dt 3 -x'
!!$     export_data(10) = drho_dtau_surf(3)
!!$     call find_drho_dtau(cctkGH,bondi_rhosurf4,drho_dtau_surf,0,N_rad_ray)
!!$     data_headers(11) = 'drho_dt 4 +z'
!!$     export_data(11) = drho_dtau_surf(1)
!!$     data_headers(12) = 'drho_dt 4 +x'
!!$     export_data(12) = drho_dtau_surf(2)
!!$     data_headers(13) = 'drho_dt 4 -x'
!!$     export_data(13) = drho_dtau_surf(3)
!!$    
!!$     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
!!$     
!!$     
!!$     !-------------------------------------------------------------------------------------------------------------!
!!$     filename = 'bondi.stu'
!!$     num_cols = 4
!!$     export_data = 0.D0
!!$     
!!$     data_headers(1) = '# Timee'
!!$     export_data(1) = CCTK_TIME
!!$    
!!$     call get_rho_info(cctkGH,12.d0,12.d0,12.d0,rho_pos,rho_old_pos,drhodt_pos)
!!$     
!!$     data_headers(2) = 'rho_pos'
!!$     export_data(2) = rho_pos
!!$     
!!$     data_headers(3) = 'rho_old_pos'
!!$     export_data(3) = rho_old_pos
!!$   
!!$     data_headers(4) = 'drhodt_pos'
!!$     export_data(4) = drhodt_pos
!!$
!!$     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

!!$     !--------------------------------------------------------------------------------!
!!$     if (1.eq.1) then
!!$        do ii=1,nsurfrho
!!$           rhosurf=rhosurfvec(ii)
!!$           
!!$           if(nsurf.eq.1) then
!!$              filename = 'bondi.fon'
!!$           else
!!$              write(filename,31)ii
!!$31            FORMAT("bondi.fon.",I1)
!!$           endif
!!$
!!$           num_cols = 5
!!$
!!$           data_headers(1) = '# Time'
!!$           export_data(1) = CCTK_TIME
!!$
!!$           data_headers(2) = "rhosurf"
!!$           export_data(2) = rhosurf
!!$           
!!$           surf_radius=outer_surf_radius
!!$           call surf_integral_fixed_rad(cctkGH,F_M0)
!!$
!!$           data_headers(3) = "outer_M0_flux"
!!$           export_data(3) = F_M0
!!$
!!$           ! M0 volume integral is computed in scheduler of diagnostics_mhd
!!$           data_headers(4) = "M0_between"
!!$           export_data(4) = M0_between_surf_VolInt*SymmFactor
!!$
!!$           ! M0 volume integral is computed in scheduler of diagnostics_mhd
!!$           data_headers(5) = "M0_AH_between"
!!$           export_data(5) = M0_between_AH_surf_VolInt*SymmFactor
!!$          
!!$           if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
!!$        end do
!!$     endif
!!$     !--------------------------------------------------------------------------------!

     filename = 'bondi.don'

     num_cols = 7
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     ! M0 volume integral is computed in scheduler of diagnostics_mhd
     data_headers(2) = "M0"
     export_data(2) = M0_VolInt*SymmFactor

     ! M0 inside AH volume integral is computed in scheduler of diagnostics_mhd
     data_headers(3) = "M0 inside AH"
     export_data(3) = M0_AH_VolInt*SymmFactor

     ! Eint (minternal_hybrid) volume integral is computed in scheduler of diagnostics_mhd
     data_headers(4) = 'Eint'
     export_data(4) = Minternal_VolInt*SymmFactor

     ! Kinetic energy (T) volume integral is computed in scheduler of diagnostics_mhd
     data_headers(5) = 'T'
     export_data(5) = T_VolInt*SymmFactor

     ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(6) = "M_ADM_Vol"
     export_data(6) = M_ADM_VolInt*SymmFactor
     data_headers(7) = "M_Komar Surf"
     call Komar_surf_integralv2(cctkGH,export_data(7))
   
     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
 
    !--------------------------------------------------------------------------------!
     filename = 'bondi.ron'

     num_cols = 10
     !num_cols=5
     export_data = 0.D0
 
     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     ! M0 flux across AH
     F_M0_BH = 0.d0
     found_horizon = 0
     horizon_number = 1
     
   
      if (cctk_iteration .gt. 0) then
         call M0_flux_through_nonstatic_BH(cctkGH,F_M0_BH,horizon_number, & 
             CCTK_TIME,t_last_bondi, xh_last_bondi,yh_last_bondi, & 
             zh_last_bondi,ah_radii_last_bondi,N_theta,N_phi, & 
             Symmetry,found_horizon)
      end if
     data_headers(2) = "F_M0_BH" 
     export_data(2) = F_M0_BH
     
     ! M0 flux across 1st rho surface
     F_M0_BH = 0.d0
     if (cctk_iteration .gt. 0) then 
        call rhosurf_integrals(cctkGH,F_M0_BH,bondi_rhosurf1,surf_area,ave_radius, & 
             CCTK_TIME,t_last_bondi, xh_last_bondi,yh_last_bondi, & 
             zh_last_bondi,surf1_radii_last_bondi,N_rad_ray,N_theta,N_phi, & 
             rhosurf_rmin,rhosurf_rmax,Symmetry,0)
     end if
     data_headers(3) = "F_M0_rhosurf1" 
     export_data(3) = F_M0_BH
          
     data_headers(4) = "ave_radius1"
     export_data(4) = ave_radius

     ! M0 flux across 2nd rho surface
     F_M0_BH = 0.d0
     if (cctk_iteration .gt. 0) then 
        call rhosurf_integrals(cctkGH,F_M0_BH,bondi_rhosurf2,surf_area,ave_radius, & 
             CCTK_TIME,t_last_bondi, xh_last_bondi,yh_last_bondi, & 
             zh_last_bondi,surf2_radii_last_bondi,N_rad_ray,N_theta,N_phi, & 
             rhosurf_rmin,rhosurf_rmax,Symmetry,0)
     end if
     data_headers(5) = "F_M0_rhosurf2" 
     export_data(5) = F_M0_BH
     
     data_headers(6) = "ave_radius2"
     export_data(6) = ave_radius
     
     ! M0 flux across 3rd rho surface
     F_M0_BH = 0.d0
     if (cctk_iteration .gt. 0) then 
        call rhosurf_integrals(cctkGH,F_M0_BH,bondi_rhosurf3,surf_area,ave_radius, & 
             CCTK_TIME,t_last_bondi, xh_last_bondi,yh_last_bondi, & 
             zh_last_bondi,surf3_radii_last_bondi,N_rad_ray,N_theta,N_phi, & 
             rhosurf_rmin,rhosurf_rmax,Symmetry,0)
     end if
     data_headers(7) = "F_M0_rhosurf3" 
     export_data(7) = F_M0_BH
   
     data_headers(8) = "ave_radius3"
     export_data(8) = ave_radius

     ! M0 flux across 4th rho surface
     F_M0_BH = 0.d0
     if (cctk_iteration .gt. 0) then 
        call rhosurf_integrals(cctkGH,F_M0_BH,bondi_rhosurf4,surf_area,ave_radius, & 
             CCTK_TIME,t_last_bondi, xh_last_bondi,yh_last_bondi, & 
             zh_last_bondi,surf4_radii_last_bondi,N_rad_ray,N_theta,N_phi, & 
             rhosurf_rmin,rhosurf_rmax,Symmetry,0)
        !reset stuff
        t_last_bondi = CCTK_TIME
        foundflag = HorizonLocalCoordinateOrigin(horizon_number,xh_last_bondi,yh_last_bondi,zh_last_bondi)
     end if
     data_headers(9) = "F_M0_rhosurf4" 
     export_data(9) = F_M0_BH
   
     data_headers(10) = "ave_radius4"
     export_data(10) = ave_radius


     
     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!

     filename = 'bondi.kon'
     num_cols = 8
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Xbh1'
     export_data(2) = bh_posn_x(1)

     data_headers(3) = 'Ybh1'
     export_data(3) = bh_posn_y(1)

     data_headers(4) = 'Zbh1'
     export_data(4) = bh_posn_z(1)

     data_headers(5) = 'X_CoM'
     export_data(5) = CoMx_VolInt*SymmFactor/M0_VolInt/SymmFactor
     data_headers(6) = 'Y_CoM'
     export_data(6) = CoMy_VolInt*SymmFactor/M0_VolInt/SymmFactor
     data_headers(7) = 'Z_CoM'
     if(Symmetry==EQUATORIAL) then
        export_data(7) = 0.D0
     else
        export_data(7) = CoMz_VolInt*SymmFactor/M0_VolInt/SymmFactor
     end if

     data_headers(8) = 'rho_b_CoM'
     call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
     pointcoords(1,1) = CoMx_VolInt*SymmFactor/M0_VolInt/SymmFactor
     pointcoords(1,2) = CoMy_VolInt*SymmFactor/M0_VolInt/SymmFactor
     pointcoords(1,3) = 0.D0
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(8))

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     !--------------------------------------------------------------------------------!
     !************ HAMILTONIAN CONSTRAINT ************!

     filename = "bondi-ham.con"
     num_cols = 7
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Haml ResN'
     export_data(2) = Ham_const_VolIntN*SymmFactor

     data_headers(3) = 'Haml ResD'
     export_data(3) = Ham_const_VolIntN*SymmFactor

     data_headers(4) = "Ham Res2P1.1N"
     export_data(4) = Ham_const_excised_VolIntN*SymmFactor
     write(*,*) "Ham_const_excised_VolIntN: ",Ham_const_excised_VolIntN
     data_headers(5) = "Ham Res2P1.1D"
     export_data(5) = Ham_const_excised_VolIntD*SymmFactor
     write(*,*) "Ham_const_excised_VolIntD: ",Ham_const_excised_VolIntD

     data_headers(6) = "HRes2P1.1NNR"
     export_data(6) = Ham_const_excised_innerregion_VolIntN*SymmFactor

     data_headers(7) = "HRes2P1.1DNR"
     export_data(7) = Ham_const_excised_innerregion_VolIntD*SymmFactor

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     !--------------------------------------------------------------------------------!
     !************ MOMENTUM CONSTRAINT ************!

     filename = 'bondi-mom.con'
     num_cols = 9
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = "Momx ResN"
     export_data(2) = momx_const_VolIntN*SymmFactor
     data_headers(3) = "Momy ResN"
     export_data(3) = momy_const_VolIntN*SymmFactor
     data_headers(4) = "Momz ResN"
     export_data(4) = momz_const_VolIntN*SymmFactor
     data_headers(5) = "Mom ResD"
     export_data(5) = mom_const_VolIntD*SymmFactor

     data_headers(6) = "Momx Res1.1N"
     export_data(6) = momx_const_excised_VolIntN*SymmFactor
     data_headers(7) = "Momy Res1.1N"
     export_data(7) = momy_const_excised_VolIntN*SymmFactor
     data_headers(8) = "Momz Res1.1N"
     export_data(8) = momz_const_excised_VolIntN*SymmFactor
     data_headers(9) = "Mom Res1.1D"
     export_data(9) = mom_const_excised_VolIntD*SymmFactor

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ GAMMA CONSTRAINT ************!
     filename = 'bondi-gam.con'
     num_cols = 4

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = "Gamx Res"
     export_data(2) = Gamx_const_VolInt*SymmFactor
     data_headers(3) = "Gamy Res"
     export_data(3) = Gamy_const_VolInt*SymmFactor
     data_headers(4) = "Gamz Res"
     export_data(4) = Gamz_const_VolInt*SymmFactor

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !!call CCTK_ReductionHandle(handle,"minimum")
     !!call CCTK_VarIndex(vindex,"grid::Z")
     !!call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Zmin,1,vindex)

     !!do jj=1,moncrief_gw_num_radii

     !!   moncrief_gw_extr_radius = moncrief_radius_GW(jj)        

     !!   ! Gravitational wave extraction!
     !!   Mass_GW = 0.0D0

     !!   ! Initialize several variables at t=0
     !!   if (cctk_iteration==0) then
     !!      momentsr_old = 0.d0
     !!      momentsi_old = 0.d0
     !!      odd_momentsr_old = 0.d0
     !!      odd_momentsi_old = 0.d0
     !!      int_momentsr = 0.d0
     !!      int_momentsi = 0.d0
     !!      momentsr = 0.D0
     !!      momentsi = 0.D0
     !!   end if

     !!   !gw_amplitudes = 0.D0
     !!   hplus = 0.d0
     !!   hcross = 0.d0

     !!   nmodes_moncrief = ubound(momentsr,1)

     !!   call gw_wave_flux(cctkGH,moncrief_gw_extr_radius,Mass_GW, &
     !!        dT,momentsr,momentsi, &
     !!        momentsr_old,momentsi_old,odd_momentsr,odd_momentsi, &
     !!        odd_momentsr_old,odd_momentsi_old,int_momentsr,int_momentsi, &
     !!        hplus,hcross,theta_GW,phi_GW, &
     !!        Zmin,dZ)

     !!   if(moncrief_gw_num_radii.eq.1) then
     !!      filename = 'bondi-evenmodes.wave'
     !!   else
     !!      write(filename,51)jj
!!51         FORMAT('bondi-evenmodes.wave.',I1)
     !!   endif

     !!   num_cols = 2*nmodes_moncrief + 1

     !!   data_headers(1) = "# Time"
     !!   export_data(1) = CCTK_TIME

     !!   do ii=1,nmodes_moncrief
     !!      data_headers(2*ii) = " " 
     !!      data_headers(2*ii+1) = " " 
     !!      export_data(2*ii) = momentsr(ii)
     !!      export_data(2*ii+1) = momentsi(ii)
     !!   end do

     !!   if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     !!   !--------------------------------------------------------------------------------!

     !!   if(moncrief_gw_num_radii.eq.1) then
     !!      filename = 'bondi-oddmodes.wave'
     !!   else
     !!      write(filename,61)jj
!!61         FORMAT('bondi-oddmodes.wave.',I1)
     !!   endif
     !!   num_cols = 2*nmodes_moncrief + 1

     !!   data_headers(1) = "# Time"
     !!   export_data(1) = CCTK_TIME

     !!   do ii=1,nmodes_moncrief
     !!      data_headers(2*ii) = " "
     !!      data_headers(2*ii+1) = " "
     !!      export_data(2*ii) = odd_momentsr(ii)
     !!      export_data(2*ii+1) = odd_momentsi(ii)
     !!   end do

     !!   if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !!end do

     !--------------------------------------------------------------------------------!
     !!filename = 'bondi.lum'
     !!num_cols = 4
     !!export_data = 0.d0
     !!data_headers = 'unknowndud'

     !!data_headers(1) = '# Time'
     !!export_data(1) = CCTK_TIME
     !!
     !!data_headers(2) = 'brem_qei'
     !!export_data(2) = brem_qei_VolInt*SymmFactor
     !!
     !!data_headers(3) = 'brem_qee'
     !!export_data(3) = brem_qee_VolInt*SymmFactor
     !!
     !!data_headers(4) = 'synchrotron'
     !!export_data(4) = synch_VolInt*SymmFactor
     !!
     !!if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     !--------------------------------------------------------------------------------
    
  end if
!end if
end subroutine statpunc_bondi_diagnostics
