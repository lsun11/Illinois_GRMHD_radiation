#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!------------------------------------
! Do TwoPuncturesAEI diagnostics
!------------------------------------
subroutine TwoPuncturesAEI_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,ext_noghost
  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,Zmin,BigM
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20
  integer                                  :: vindex
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  !Needed for GW extraction:
  real*8                                   :: Mass_GW,hplus,hcross,moncrief_gw_extr_radius
  real*8                                   :: psi4r,psi4i
  integer                                  :: jj,nmodes_moncrief
  real*8,dimension(2)                      :: gw_amplitudes
  !End GW extraction variables  
  real*8                                   :: old_surf_radius,SymmFactor
  CCTK_POINTER :: reduction_value_pointer
  integer :: header_flag,handle
  integer :: index,num_cols
  integer :: ierr,myproc_rank,ii

  ! Variables used for BH local diagnotsics:
  integer :: found_horizon, horizon_number
  integer :: N_theta_bh_diagnostics, N_phi_bh_diagnostics
  real*8  :: Mirr, M_BH, J_BH, xi_err

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(trA_detg_enforce.ne.2) then
     !First perform sanity check.  Already done if trA_detg_enforce==2!
     !     call sanitycheck_restore_Aij(cctkGH,cctk_lsh,gxx,gxy,gxz,gyy,gyz,gzz, &
     !          Axx,Axy,Axz,Ayy,Ayz,Azz)
  end if

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(1),1,index)

  if(MOD(cctk_iteration,out_every)==0) then

     write(*,*) "Setting BigM to 1.0!"
     BigM = 1.D0

     call set_SymmFactor(Symmetry,SymmFactor)

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1
     else
        header_flag = 0
     end if

     !--------------------------------------------------------------------------------!

     filename = 'binary_bh.mon'
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
     Psi = exp(phi)
     call CCTK_VarIndex (vindex, "bssn::psi")
     call interp_center_carp(cctkGH,vindex,export_data(5))

     data_headers(6) = 'lapsemin'
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(index,"lapse::lapm1")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(6),1,index)
     export_data(6) = export_data(6) + 1.D0

     data_headers(7) = 'phimax'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(index,"bssn::phi")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(7),1,index)


     data_headers(8) = 'Jsurf_BH1_1.1M'
     old_surf_radius = surf_radius
     surf_radius=1.1D0*BigM
     call J_surf_integral_offcenter(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(8))
     surf_radius = old_surf_radius

     data_headers(9) = 'Jsurf_BH2_1.1M'
     old_surf_radius = surf_radius
     surf_radius=1.1D0*BigM
     call J_surf_integral_offcenter(cctkGH,bh_posn_x(2),bh_posn_y(2),bh_posn_z(2),export_data(9))
     surf_radius = old_surf_radius

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(10) = "J_1.1exVol"
     export_data(10) = J_ADM_VolInt*SymmFactor

     data_headers(11) = "M_ADMsurf_BH1_1.1M"
     old_surf_radius = surf_radius
     surf_radius=1.1D0*BigM
     call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(11))
     surf_radius = old_surf_radius

     data_headers(12) = "M_ADMsurf_BH2_1.1M"
     old_surf_radius = surf_radius
     surf_radius=1.1D0*BigM
     call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(2),bh_posn_y(2),bh_posn_z(2),export_data(12))
     surf_radius = old_surf_radius

     ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(13) = "M_ADM1.1exVol"
     export_data(13) = M_ADM_VolInt*SymmFactor

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!

     filename = 'binary_bh.kon'
     num_cols = 7
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Xbh1'
     export_data(2) = bh_posn_x(1)

     data_headers(3) = 'Ybh1'
     export_data(3) = bh_posn_y(1)

     data_headers(4) = 'Zbh1'
     export_data(4) = bh_posn_z(1)

     data_headers(5) = 'Xbh2'
     export_data(5) = bh_posn_x(2)

     data_headers(6) = 'Ybh2'
     export_data(6) = bh_posn_y(2)

     data_headers(7) = 'Zbh2'
     export_data(7) = bh_posn_z(2)

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     ! BH local diagnostics, using isolated horizon formalism
     !
     filename = 'binary_bh_IHF_diagnostics.mon'
     num_cols = 13
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     N_theta_bh_diagnostics = 50
     N_phi_bh_diagnostics = 50

     Mirr = 0.d0
     M_BH = 0.d0
     J_BH = 0.d0
     xi_err = 0.d0

     horizon_number = 1
     call BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, &
                       horizon_number,N_theta_bh_diagnostics, &
                       N_phi_bh_diagnostics,Symmetry,found_horizon)

     data_headers(2) = 'Mirr_bh1'
     export_data(2) = Mirr

     data_headers(3) = 'M_bh1'
     export_data(3) = M_BH

     data_headers(4) = 'J_bh1'
     export_data(4) = J_BH

     data_headers(5) = 'xi_err1'
     export_data(5) = xi_err

     horizon_number = 2
     call BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, &
                       horizon_number,N_theta_bh_diagnostics, &
                       N_phi_bh_diagnostics,Symmetry,found_horizon)

     data_headers(6) = 'Mirr_bh2'
     export_data(6) = Mirr

     data_headers(7) = 'M_bh2'
     export_data(7) = M_BH

     data_headers(8) = 'J_bh2'
     export_data(8) = J_BH

     data_headers(9) = 'xi_err2'
     export_data(9) = xi_err

     ! Common horizon?
     horizon_number = 3
     call BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, &
                       horizon_number,N_theta_bh_diagnostics, &
                     N_phi_bh_diagnostics,Symmetry,found_horizon)

     data_headers(10) = 'Mirr_bh3'
     export_data(10) = Mirr

     data_headers(11) = 'M_bh3'
     export_data(11) = M_BH

     data_headers(12) = 'J_bh3'
     export_data(12) = J_BH

     data_headers(13) = 'xi_err3'
     export_data(13) = xi_err

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ HAMILTONIAN CONSTRAINT ************!

     filename = 'binary_bh-ham.con'
     num_cols = 5
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Haml ResN'
     export_data(2) = Ham_const_VolIntN*SymmFactor

     data_headers(3) = 'Haml ResD'
     export_data(3) = Ham_const_VolIntD*SymmFactor

     data_headers(4) = "Ham Res2P1.1N"
     export_data(4) = Ham_const_excised_VolIntN*SymmFactor

     data_headers(5) = "Ham Res2P1.1D"
     export_data(5) = Ham_const_excised_VolIntD*SymmFactor

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ MOMENTUM CONSTRAINT ************!

     filename = 'binary_bh-mom.con'
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
     filename = 'binary_bh-gam.con'
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
     do ii=1,nsurf

        surf_radius=surfvec(ii)
        rbr=rbrvec(ii)
        drbdr=drbrvec(ii)
        ddrbddr=ddrbrvec(ii)

        if(nsurf.eq.1) then
           filename = 'binary_bh.jon'
        else
           write(filename,31)(ii-1)
31         FORMAT('binary_bh.jon.',I1)
        endif

        num_cols = 9

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

        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end do

     !--------------------------------------------------------------------------------!
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(index,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Zmin,1,index)

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

        call gw_wave_flux(cctkGH,moncrief_gw_extr_radius,Mass_GW, &
             dT,momentsr,momentsi, &
             momentsr_old,momentsi_old,odd_momentsr,odd_momentsi, &
             odd_momentsr_old,odd_momentsi_old,int_momentsr,int_momentsi, &
             hplus,hcross,theta_GW,phi_GW, &
             Zmin,dZ)

        if(moncrief_gw_num_radii.eq.1) then
           filename = 'binary_bh-evenmodes.wave'
        else
           write(filename,41)jj
41         FORMAT('binary_bh-evenmodes.wave.',I1)
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

        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
        !--------------------------------------------------------------------------------!

        if(moncrief_gw_num_radii.eq.1) then
           filename = 'binary_bh-oddmodes.wave'
        else
           write(filename,51)jj
51         FORMAT('binary_bh-oddmodes.wave.',I1)
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

        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !--------------------------------------------------------------------------------!
     if(1==0) then
        filename = 'binary_bh-lapmin.mon'
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
     end if
  end if

  !  PsiTau = (shiftx*lapsex + shifty*lapsey + shiftz*lapsez)/(2.D0 * (1.D0 + lapm1)) - trK

end subroutine TwoPuncturesAEI_diagnostics
