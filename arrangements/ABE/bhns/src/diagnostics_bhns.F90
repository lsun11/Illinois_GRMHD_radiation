#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for magnetar thorn
!-----------------------------------------------------------------------------
subroutine magnetar_diagnostics_bhns(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,Zmin
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

  real*8                                   :: old_surf_radius,SymmFactor,magnetar_admmass
  real*8                                   :: rho_tiny,F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em,RestMass
  CCTK_POINTER :: reduction_value_pointer
  integer :: header_flag,handle
  integer :: index,num_cols
  integer :: ierr,myproc_rank,ii

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(cctk_iteration.ge.0 .and. 1==0) then
     !Set position in CarpetRegrid2
     magnetar_admmass = 0.278798248452541
     position_x(1) = 1.79*magnetar_admmass*sin(0.168d0/magnetar_admmass*cctk_time)
     position_y(1) = 1.43*magnetar_admmass*(1.D0 - cos(0.168d0/magnetar_admmass*cctk_time))
     position_z(1) = 0.D0
     write(*,*) "NEW AMR box POSITION: ",position_x(1),position_y(1)
  end if


  if(MOD(cctk_iteration,out_every)==0) then
     call set_SymmFactor(Symmetry,SymmFactor)

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1
     else
        header_flag = 0
     end if

     !--------------------------------------------------------------------------------!

     filename = 'magnetar.mon'
     num_cols = 7
     export_data = 0.D0

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

     data_headers(7) = 'rhocent'
     call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
     call interp_center_carp(cctkGH,vindex,export_data(7))

     ! Sadly, we can't set local gridfunctions (like Pr here) inside GLOBAL routines like this one.
     !     data_headers(8) = 'Omega_c'
     !     Pr = (vy*X - vx*Y) / (X*X + Y*Y)
     !     call CCTK_VarIndex (vindex, "mhd_evolve::Pr")
     !     call interp_center_carp(cctkGH,vindex,export_data(8))

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!

     filename = 'magnetar.don'

     if(em_evolve_enable==1) then
        num_cols = 15
     else
        num_cols = 6
     end if
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     ! M0 volume integral is computed in scheduler of diagnostics_mhd
     data_headers(2) = "M0"
     export_data(2) = M0_VolInt*SymmFactor

     ! Eint (minternal_hybrid) volume integral is computed in scheduler of diagnostics_mhd
     data_headers(3) = 'Eint'
     export_data(3) = Minternal_VolInt*SymmFactor

     ! Kinetic energy (T) volume integral is computed in scheduler of diagnostics_mhd
     data_headers(4) = 'T'
     export_data(4) = T_VolInt*SymmFactor

     ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(5) = "M_ADM_Vol"
     export_data(5) = M_ADM_VolInt*SymmFactor

     data_headers(6) = "M_Komar Surf"
     call Komar_surf_integralv2(cctkGH,export_data(6))

     if(em_evolve_enable==1) then
        ! EM energy volume integral is computed in scheduler of diagnostics_mhd
        data_headers(7) = 'Eem'
        export_data(7) = em_energy_VolInt*SymmFactor

        ! (Alt.) EM energy volume integral is computed in scheduler of diagnostics_mhd
        data_headers(8) = 'Eem2'
        export_data(8) = em_energy2_VolInt*SymmFactor

        ! b_phi volume integral is computed in scheduler of diagnostics_mhd
        data_headers(9) = 'int B_phi'
        export_data(9) = b_phi_VolInt*SymmFactor

        data_headers(10) = 'b2_max'
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex (vindex, "mhd_evolve::temp8")
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(10),1,vindex)

        data_headers(11) = 'b2_max_x'
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex (vindex, "mhd_evolve::temp9")
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(11),1,vindex)

        data_headers(12) = 'b2_max_y'
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex (vindex, "mhd_evolve::temp10")
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(12),1,vindex)

        data_headers(13) = 'b2_max_z'
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex (vindex, "mhd_evolve::temp11")
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(13),1,vindex)

        data_headers(14) = "Monopole"
        export_data(14) = monopole_VolInt*SymmFactor

        data_headers(15) = "Monop.t=0"
        export_data(15) = initial_monopole_value*SymmFactor

     end if

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !--------------------------------------------------------------------------------!
     filename = 'magnetar.jon'
     num_cols = 11

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(2) = "AngMomVol"
     export_data(2) = J_ADM_VolInt*SymmFactor

     data_headers(3) = "J_surf"
     call J_surf_integral_offcenter(cctkGH,0.D0,0.D0,0.D0,export_data(3))

     data_headers(4) = "Mass_sur"
     call M_surf_integral(cctkGH,export_data(4))

     data_headers(5) = "Komar Mass(NULL)"
     export_data(5) = 0.D0
     ! Sadly, we can't set local gridfunctions (like lapsex,lapsey,lapsez here) inside GLOBAL routines like this one.
     !     call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry)
     !     call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)
     !     call Komar_surf_integralv2(cctkGH,export_data(5))

     data_headers(6) = "out_M0_flux"
     data_headers(7) = "out_fluid_Eng_flux"
     data_headers(8) = "out_EM_Energy_flux"
     data_headers(9) = "out_fluid_J_flux"
     data_headers(10) = "out_EM_J_flux"
     if(excision_enable==0) then
        call mhd_surface_fluxes(cctkGH,export_data(6),export_data(7),export_data(8), &
             export_data(9),export_data(10))
     else
        write(*,*) "OUT FLUXES NOT SUPPORTED ACROSS HORIZON (YET)"
        !FIXME:
        !        call ahfluxes(CCTK_PASS_FTOF)
     end if

     !FIXME: not sure what this was...
     data_headers(11) = "rho_isoflux"
     !  call rhoflux(CCTK_PASS_FTOF)
     export_data(11) = 0.D0

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !--------------------------------------------------------------------------------!
     !************ HAMILTONIAN CONSTRAINT ************!

     filename = 'magnetar-ham.con'
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

     filename = 'magnetar-mom.con'
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
     filename = 'magnetar-gam.con'
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
     filename = 'magnetar.ron'
     num_cols = 5

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = "Gammaxmax"
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(index,"bssn::Gammax")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(2),1,index)

     data_headers(3) = "rho_bmax"
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(index,"mhd_evolve::rho_b")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(3),1,index)

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(4) = "AngMomVol"
     export_data(4) = J_ADM_VolInt*SymmFactor

     data_headers(5) = "lapse_min"
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(index,"lapse::lapm1")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(5),1,index)
     export_data(5) = export_data(5) + 1.D0

     ! Sadly, we can't set local gridfunctions inside GLOBAL routines like this one, since ext is not defined!
     !     data_headers(6) = "Rp"
     !     call CCTK_ReductionHandle(handle,"maximum")
     !     call CCTK_VarIndex(index,"mhd_evolve::rho_b")
     !     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,rho_tiny,1,index)
     !     rho_tiny = rho_tiny * 1.D-3
     !     call find_equatorial_radius(ext,cctkgh,X,Y,Z,dX,dY,dZ,rho_tiny,rho_star,phi,gxx,gxy,gyy,Symmetry,export_data(6),CCTK_MyProc(CCTKGH))

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !--------------------------------------------------------------------------------!
     filename = 'magnetar.xon'
     num_cols = 5

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = "Xcen"
     data_headers(3) = "Ycen"
     data_headers(4) = "Zcen"

     ! M0 volume integral is computed in scheduler of diagnostics_mhd
     RestMass = M0_VolInt*SymmFactor

     ! CoM volume integrals are computed in scheduler of diagnostics_mhd
     export_data(2) = SymmFactor*CoMx_VolInt/RestMass
     export_data(3) = SymmFactor*CoMx_VolInt/RestMass
     export_data(4) = SymmFactor*CoMx_VolInt/RestMass

     data_headers(5) = "Rcen"
     export_data(5) = sqrt(export_data(2)*export_data(2) + export_data(3)*export_data(3) + export_data(4)*export_data(4))

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

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
           filename = 'magnetar-evenmodes.wave'
        else
           write(filename,41)jj
41         FORMAT('magnetar-evenmodes.wave.',I1)
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
           filename = 'magnetar-oddmodes.wave'
        else
           write(filename,51)jj
51         FORMAT('magnetar-oddmodes.wave.',I1)
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
     !--------------------------------------------------------------------------------!
  end if
  write(*,*) "Finished with magnetar diagnostics in bhns"	
end subroutine magnetar_diagnostics_bhns
