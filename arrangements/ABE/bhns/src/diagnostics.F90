#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for bhns thorn
!-----------------------------------------------------------------------------
subroutine bhns_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i,j,k


  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,Zmin,dt_diag
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20

  integer                                  :: vindex,N_points
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  !Needed for convenient constraint output
  real*8 :: constraint_ratio
  !End of convenient constraint output variables

  !Needed for GW extraction:
  real*8                                   :: Mass_GW,hplus,hcross,moncrief_gw_extr_radius
  real*8                                   :: psi4r,psi4i,r_local
  integer                                  :: jj,nmodes_moncrief
  real*8,dimension(2)                      :: gw_amplitudes
  !End GW extraction variables  

  ! Variables used for BH local diagnotsics:
  integer :: found_horizon, horizon_number
  integer :: N_theta_bh_diagnostics, N_phi_bh_diagnostics
  real*8  :: Mirr, M_BH, J_BH, xi_err
  ! End BH local diagnotsics variables

  real*8                                   :: old_surf_radius,SymmFactor,horizdirn_x,horizdirn_y,horizdirn_z,out_radius
  real*8                                   :: rho_tiny,M0dot,F_E_fluid,F_E_em,F_J_fluid,F_J_em,RestMass
  real*8, dimension(1,3)                   :: pointcoords
  CCTK_POINTER :: reduction_value_pointer
  integer :: header_flag,handle
  integer :: num_cols
  integer :: ierr,myproc_rank,ii,narray

  !center of mass
  real*8,dimension(2,3)                    :: points_CM
  real*8  :: max_horiz_radius

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  do k=1, cctk_lsh(3)
     do j=1, cctk_lsh(2)
        do i=1, cctk_lsh(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1Inside bhns: diagnostics.F90, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2Inside bhns: diagnostics.F90, phi(i,j,k) is ", phi(i,j,k)
end if

    end do
    end do
    end do


!if(0==1) then

  if(mod(cctk_iteration,out_every)==0) then
     call set_SymmFactor(Symmetry,SymmFactor)

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1

        if (ntot_bhns .ne. N_theta*N_phi) then
           write(*,*) 'In bhns, ntot_bhns = ',ntot_bhns
           write(*,*) 'N_theta, N_phi = ',N_theta, N_phi
           write(*,*) 'ntot_bhns must be set to N_theta*N_phi'
           write(*,*) 'Need to reset the parameter. Bye now...'
           stop
        end if

     else
        header_flag = 0
     end if

     !-----------------------------------------!
     ! choose what planes we want 2D output
     !if((use_new_bhns_initial_data.ne.3).or.(abs(bh_posn_y(1))>0.d0)) then
     if((use_new_bhns_initial_data.eq.0).or.(use_new_bhns_initial_data.eq.1) &
           .or.(use_new_bhns_initial_data.eq.2).or.(use_new_bhns_initial_data.eq.6).or.(abs(bh_posn_y(1))>0.d0)) then
        out_xzplane_y = bh_posn_y(1)
        out_yzplane_x = bh_posn_x(1)
        out_xyplane_z = bh_posn_z(1)
        print *,"out_xzplane_y =",out_xzplane_y
        print *,"out_yzplane_x =",out_yzplane_x
        print *,"out_xyplane_z =",out_xyplane_z
     end if
     !--------------------------------------------------------------------------------!
     filename = 'bhns.mon'
     num_cols = 19
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

     data_headers(8) = 'rhocent'
     call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
     call interp_center_carp(cctkGH,vindex,export_data(8))

     data_headers(9) = 'rho_b_max'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(9),1,vindex)

     export_data(9)  = bhns_max_rho_b
     print *, "bhns_max_rho_b=", bhns_max_rho_b


     data_headers(10) = 'Jsurf_BH_1.1'
     surf_radius=BH_Vol_Excise_Radius
     old_surf_radius = surf_radius
!     surf_radius=1.1D0
     surf_radius=BH_Vol_Excise_Radius 
     call J_surf_integral_offcenter(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(10))
     surf_radius = old_surf_radius

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(11) = "J_1.1exVol"
     export_data(11) = J_ADM_VolInt*SymmFactor
!     print *,"J_ADM_VolInt=",J_ADM_VolInt,"SymmFactor=",SymmFactor

     data_headers(12) = "M_ADMsurf_BH_1.1"
     surf_radius=BH_Vol_Excise_Radius
     old_surf_radius = surf_radius
!     surf_radius=1.1D0
     surf_radius=BH_Vol_Excise_Radius 
     call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(12))
     surf_radius = old_surf_radius

     ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(13) = "M_ADM1.1exVol"
     export_data(13) = M_ADM_VolInt*SymmFactor

     data_headers(14) = 'rhostar_min'
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex (vindex, "mhd_evolve::rho_star")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(14),1,vindex)

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(15) = "J_1.1exVolr"
     export_data(15) = J_ADM_VolInt_inner*SymmFactor

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(16) = "J_1.1exVolr/2"
     export_data(16) = J_ADM_VolInt_inner2*SymmFactor

     !     data_headers(8) = 'Omega_c'
     !     Pr = (vy*X - vx*Y) / (X*X + Y*Y)
     !     call CCTK_VarIndex (vindex, "mhd_evolve::Pr")
     !     call interp_center_carp(cctkGH,vindex,export_data(8))

        data_headers(17) = 'rho_b_max_x'
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex (vindex, "mhd_evolve::temp12")
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(17),1,vindex)

        data_headers(18) = 'rho_b_max_y'
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex (vindex, "mhd_evolve::temp13")
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(18),1,vindex)
        
        data_headers(19) = 'rho_b_max_z'
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex (vindex, "mhd_evolve::temp14")
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(19),1,vindex)
        
!        print *, "here=", rho_b_max_x,rho_b_max_y,rho_b_max_z


     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     do ii=1,nsurf

        surf_radius=surfvec(ii)
        rbr=rbrvec(ii)
        drbdr=drbrvec(ii)
        ddrbddr=ddrbrvec(ii)

        if(nsurf.eq.1) then
           filename = 'bhns.jon'
        else
           write(filename,31)(ii-1)
31         FORMAT("bhns.jon.",I1)
        endif

        num_cols = 21

        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME

        data_headers(2) = "fisheye radius"
        export_data(2) = surf_radius

        data_headers(3) = "phys radius"
        export_data(3) = surf_radius*rbr

        data_headers(4) = "Mass_sur"
        call M_surf_integral(cctkGH,export_data(4))

        data_headers(5) = "Ang_mom_surf"
!        call J_surf_integral_offcenter(cctkGH,0.D0,0.D0,0.D0,export_data(5))
        call J_surf_integral_offcenter(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(5))

        data_headers(6) = "Komar Mass"
        call Komar_surf_integralv2(cctkGH,export_data(6))

        data_headers(7) = "P_x2"
        data_headers(8) = "P_y2"
        data_headers(9) = "P_z2"
        call P_surf_integralv2(cctkGH,export_data(7),export_data(8),export_data(9))

        !call mhd_surface_fluxes(cctkGH,F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em)
        data_headers(10) = "F_M0"
        data_headers(11) = "F_E_fluid"
        data_headers(12) = "F_E_em"
        data_headers(13) = "F_J_fluid"
        data_headers(14) = "F_J_em"
        call mhd_surface_fluxes(cctkGH,export_data(10),export_data(11),export_data(12),export_data(13),export_data(14))

	data_headers(15) = "F_E_rad"
        data_headers(16) = "F_J_rad"
        data_headers(17) = "F_E_rad_nue"
        data_headers(18) = "F_J_rad_nue"
        data_headers(19) = "F_E_rad_nux"
        data_headers(20) = "F_J_rad_nux"
	call rad_surface_fluxes(cctkGH,export_data(15),export_data(16),export_data(17),export_data(18),export_data(19),export_data(20),ii,cctk_time)
        
        data_headers(21) =  "F_Mesc"
        call mhd_surface_fluxes_2(cctkGH, export_data(21))
	if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end do

  
     !--------------------------------------------------------------------------------!
if(0==1) then   
     filename = 'bhns_spin.mon'
        
     N_theta_bh_diagnostics =  N_theta
     N_phi_bh_diagnostics   =  N_phi

     print *, "N_theta_bh_diagnostics=", N_theta
     print *, "N_phi_bh_diagnostics=", N_phi

     if (track_bhns.eq.1) then
        
        num_cols =9
        export_data = 0.D0

        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME

        if(CCTK_ITERATION.eq. 0) then
           points_CM(1,1) =  initial_ns_coord_x
           points_CM(1,2) =  initial_ns_coord_y
           
           points_CM(2,1) =  initial_ns2_coord_x
           points_CM(2,2) =  initial_ns2_coord_y
        else
           points_CM(1,1) = Box1X_VolInt1/Box1denom_VolInt1
           points_CM(1,2) = Box1Y_VolInt1/Box1denom_VolInt1
           
           points_CM(2,1) = Box1X_VolInt2/Box1denom_VolInt2
           points_CM(2,2) = Box1Y_VolInt2/Box1denom_VolInt2

        end if
        points_CM(1,3) = 0.D0
        points_CM(2,3) = 0.D0


        print *,"points_CM(1,1)=",points_CM(1,1)
        print *,"points_CM(1,2)=",points_CM(1,2)
        print *,"points_CM(1,3)=",points_CM(1,3)

        print *,"points_CM(2,1)=",points_CM(2,1)
        print *,"points_CM(2,2)=",points_CM(2,2)
        print *,"points_CM(2,3)=",points_CM(2,3)

        ! angular momentum star1
        data_headers(2) = "J_star1"
        r_local=bhns_R_NS
        call J_surf_integral_spin(cctkGH,N_theta_bh_diagnostics,N_phi_bh_diagnostics,points_CM(1,1),points_CM(1,2),points_CM(1,3),r_local,export_data(2),Symmetry)

        ! angular momentum star2
        data_headers(3) = "J_star2"
        r_local=bhns_R_NS
        call J_surf_integral_spin(cctkGH,N_theta_bh_diagnostics,N_phi_bh_diagnostics,points_CM(2,1),points_CM(2,2),points_CM(2,3),r_local,export_data(3),Symmetry)


        ! angular momentum star1
        data_headers(4) = "J_star1 R at 1.1"
        r_local=1.1*bhns_R_NS
        call J_surf_integral_spin(cctkGH,N_theta_bh_diagnostics,N_phi_bh_diagnostics,points_CM(1,1),points_CM(1,2),points_CM(1,3),r_local,export_data(4),Symmetry)

        ! angular momentum star2
        data_headers(5) = "J_star2 R at 1.1"
        r_local=1.1*bhns_R_NS
        call J_surf_integral_spin(cctkGH,N_theta_bh_diagnostics,N_phi_bh_diagnostics,points_CM(2,1),points_CM(2,2),points_CM(2,3),r_local,export_data(5),Symmetry)


        ! circulation star1
        data_headers(6) = "C star1"
        r_local = bhns_R_NS
        call circulation(cctkGH,points_CM(1,1),points_CM(1,2),points_CM(1,3),r_local,export_data(6))

        ! circulation star1
        data_headers(7) = "C star1 along y"
        r_local = bhns_R_NS2
        call circulation(cctkGH,points_CM(1,1),points_CM(1,2),points_CM(1,3),r_local,export_data(7))

        ! circulation star1
        data_headers(8) = "C star1 along y X 0.9"
        r_local = 0.9*bhns_R_NS2 
        call circulation(cctkGH,points_CM(1,1),points_CM(1,2),points_CM(1,3),r_local,export_data(8))

        ! circulation star1
        data_headers(9) = "C star1 along y X  0.85"
        r_local = 0.85*bhns_R_NS2
        call circulation(cctkGH,points_CM(1,1),points_CM(1,2),points_CM(1,3),r_local,export_data(9))


        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

!        r_local = bhns_R_NS
!        call circulation_trace_particles(cctkGH,iteration,particle_tracer_substep_every,CCTK_DELTA_TIME,x_cm,y_cm,z_cm,r_local,output_integral)

     else 

        if (Symmetry==NO_SYMM)  then

        num_cols =11
        export_data = 0.D0

        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME

        !     Get horizon radius in direction where it is likely to be MAXIMIZED:
        horizdirn_x = 0.D0
        horizdirn_y = 100000.D0
        horizdirn_z = 0.D0
        call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,max_horiz_radius)
        
        ! radius
        data_headers(2) = 'BH_radius'
        export_data(2)  =  max_horiz_radius

        ! angular momentum components
        data_headers(3) = "Jx at r_BH"
        data_headers(4) = "Jy at r_BH"
        data_headers(5) = "Jz at r_BH"
        r_local=max_horiz_radius
        call J_surf_integral_spin(cctkGH,N_theta_bh_diagnostics,N_phi_bh_diagnostics,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),r_local,export_data(3),export_data(4),export_data(5),Symmetry)

        ! angular momentum components
        data_headers(6) = "Jx at r = 100M"
        data_headers(7) = "Jy at r = 100M"
        data_headers(8) = "Jz at r = 100M"
        r_local= 100.d0*M_ADM   

        call J_surf_integral_spin(cctkGH,N_theta_bh_diagnostics,N_phi_bh_diagnostics,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),r_local,export_data(6),export_data(7),export_data(8),Symmetry)


        ! angular momentum components
        data_headers(9)  = "Jx at r = 200M"
        data_headers(10) = "Jy at r = 200M"
        data_headers(11) = "Jz at r = 200M"
        r_local= 200.d0*M_ADM   

        call J_surf_integral_spin(cctkGH,N_theta_bh_diagnostics,N_phi_bh_diagnostics,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),r_local,export_data(9),export_data(10),export_data(11),Symmetry)


        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     else
        print *,"For equatorial symmetry use Isolate horizon formalism"
     end if
  end if
     !--------------------------------------------------------------------------------!
end if
     filename = 'bhns.don'

     num_cols = 14
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     ! M0 volume integral is computed in scheduler of diagnostics_mhd
     data_headers(2) = "M0"
     export_data(2) = M0_VolInt*SymmFactor

     ! M0 inside AH volume integral is computed in scheduler of diagnostics_mhd
     data_headers(3) = "M0 inside AH"
     export_data(3) = M0_AH_VolInt*SymmFactor

     ! Mass of escaping matter for r> 30M integral is computed in scheduler of diagnostics_mhd
     data_headers(4) = "M0_esc(r>30M)"
     export_data(4) = M0_escape30M*SymmFactor

     ! Mass of escaping matter for r> 50M integral is computed in scheduler of diagnostics_mhd
     data_headers(5) = "M0_esc(r>50M)"
     export_data(5) = M0_escape50M*SymmFactor

     ! Mass of escaping matter for r> 70M integral is computed in scheduler of diagnostics_mhd
     data_headers(6) = "M0_esc(r>70M)"
     export_data(6) = M0_escape70M*SymmFactor

     ! Mass of escaping matter for r> 100M integral is computed in scheduler of diagnostics_mhd
     data_headers(7) = "M0_esc(r>100M)"
     export_data(7) = M0_escape100M*SymmFactor

     ! Eint (minternal_hybrid) volume integral is computed in scheduler of diagnostics_mhd
     data_headers(8) = 'Eint'
     export_data(8) = Minternal_VolInt*SymmFactor

     ! Kinetic energy (T) volume integral is computed in scheduler of diagnostics_mhd
     data_headers(9) = 'T'
     export_data(9) = T_VolInt*SymmFactor

     ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(10) = "M_ADM_Vol"
     export_data(10) = M_ADM_VolInt*SymmFactor

     data_headers(11) = "M_Komar Surf"
     call Komar_surf_integralv2(cctkGH,export_data(11))

     ! M0 flux across AH
     found_horizon = 0
     horizon_number = 1
     M0dot = 0.d0
     dt_diag = CCTK_TIME - t_last
     if (cctk_iteration .gt. 0) then
        call M0_flux_through_nonstatic_BH(cctkGH,M0dot,horizon_number, &
             CCTK_TIME,t_last, xh_last,yh_last, &
             zh_last,ah_radii_last,N_theta,N_phi, &
             Symmetry,found_horizon)
     end if
     data_headers(12) = "M0dot_BH"
     export_data(12) = M0dot

     ! Total amount of M0 flowing into BH
     if (found_horizon==1) then
        ! Note: the minus sign is to make int_M0dot positive
        int_M0dot = int_M0dot - 0.5d0*(M0dot+M0dot_last)*dt_diag
        M0dot_last = M0dot
     end if
     data_headers(13) = "int_M0dot"
     export_data(13) = int_M0dot

     ! fluid energy computed from Integrate d^3x sqrt(-g) Tfluid_0^0, computed in scheduler of diagnostics_mhd
     data_headers(14) = "Tfluid0_0_VolInt"
     export_data(14) = Tfluid0_0_VolInt*SymmFactor

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!

     !if((use_new_bhns_initial_data.ne.3)) then
     if((use_new_bhns_initial_data.eq.0).or.(use_new_bhns_initial_data.eq.1) &
           .or.(use_new_bhns_initial_data.eq.2)) then
        filename = 'bhns.kon'
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
        export_data(5) = CoMx_VolInt/CoM_VolInt_denominator
        data_headers(6) = 'Y_CoM'
        export_data(6) = CoMy_VolInt/CoM_VolInt_denominator
        data_headers(7) = 'Z_CoM'
        if(Symmetry==EQUATORIAL) then
           export_data(7) = 0.D0
        else
           export_data(7) = CoMz_VolInt/CoM_VolInt_denominator
        end if
        
        data_headers(8) = 'rho_b_CoM'
        call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
        pointcoords(1,1) = CoMx_VolInt/CoM_VolInt_denominator
        pointcoords(1,2) = CoMy_VolInt/CoM_VolInt_denominator
        pointcoords(1,3) = 0.D0
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(8))
        
        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
        filename = 'bhns.pon'
        num_cols = 8
        export_data = 0.D0
        
        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME
        
        data_headers(2) = 'phicent'
        call CCTK_VarIndex(vindex,"bssn::phi")
        
        pointcoords(1,1) = bh_posn_x(1)
        pointcoords(1,2) = bh_posn_y(1)
        pointcoords(1,3) = 0.D0
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(2))
        
        data_headers(3) = 'phixplus'
        call CCTK_VarIndex(vindex,"bssn::phi")
        horizdirn_x = 100000.D0
        horizdirn_y = 0.D0
        horizdirn_z = 0.D0
        call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
        write(*,*) "inside bhns diagnostics: out_radius = ",out_radius
        pointcoords(1,1) = bh_posn_x(1)+out_radius
        pointcoords(1,2) = bh_posn_y(1)
        pointcoords(1,3) = 0.D0
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(3))
        
        data_headers(4) = 'phiyplus'
        call CCTK_VarIndex(vindex,"bssn::phi")
        horizdirn_x = 0.D0
        horizdirn_y = 100000.D0
        horizdirn_z = 0.D0
        call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
        pointcoords(1,1) = bh_posn_x(1)
        pointcoords(1,2) = bh_posn_y(1)+out_radius
        pointcoords(1,3) = 0.D0
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(4))
        
        data_headers(5) = 'phixminus'
        call CCTK_VarIndex(vindex,"bssn::phi")
        horizdirn_x = -100000.D0
        horizdirn_y = 0.D0
        horizdirn_z = 0.D0
        call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
        pointcoords(1,1) = bh_posn_x(1)-out_radius
        pointcoords(1,2) = bh_posn_y(1)
        pointcoords(1,3) = 0.D0
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(5))
        
        data_headers(6) = 'phiyminus'
        call CCTK_VarIndex(vindex,"bssn::phi")
        horizdirn_x = 0.D0
        horizdirn_y = -100000.D0
        horizdirn_z = 0.D0
        call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
        pointcoords(1,1) = bh_posn_x(1)
        pointcoords(1,2) = bh_posn_y(1)-out_radius
        pointcoords(1,3) = 0.D0
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(6))
        
        data_headers(7) = 'phizplus'
        call CCTK_VarIndex(vindex,"bssn::phi")
        horizdirn_x = 0.D0
        horizdirn_y = 0.D0
        horizdirn_z = 100000.D0
        call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
        pointcoords(1,1) = bh_posn_x(1)
        pointcoords(1,2) = bh_posn_y(1)
        pointcoords(1,3) = out_radius
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(7))
        
        data_headers(8) = 'phixminus0.5'
        call CCTK_VarIndex(vindex,"bssn::phi")
        horizdirn_x = -100000.D0
        horizdirn_y = 0.D0
        horizdirn_z = 0.D0
        call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,out_radius);
        pointcoords(1,1) = bh_posn_x(1)-out_radius*0.5D0
        pointcoords(1,2) = bh_posn_y(1)
        pointcoords(1,3) = 0.D0
        call interp_driver_carp(cctkGH,1,pointcoords,vindex,export_data(8))
        
        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)     

     end if
     !--------------------------------------------------------------------------------!
     !Box1X_VolInt,Box1Y_VolInt,Box1Z_VolInt,Box1denom_VolInt
     filename = 'bhns.xon'
     !if(use_new_bhns_initial_data.ne.3) then
     if((use_new_bhns_initial_data.eq.0).or.(use_new_bhns_initial_data.eq.1) &
           .or.(use_new_bhns_initial_data.eq.2)) then
        num_cols = 3
        export_data = 0.D0
        
        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME
        
        data_headers(2) = 'Box1x'
        export_data(2) = position_x(2) 
        data_headers(3) = 'Box1y'
        export_data(3) = position_y(2)
     else
        num_cols = 7
        export_data = 0.D0
        
        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME
        
        data_headers(2) = 'Box1x'
        export_data(2) = position_x(1) 
        data_headers(3) = 'Box1y'
        export_data(3) = position_y(1)
        
        data_headers(4) = 'Box2x'
        export_data(4) = position_x(2)
        data_headers(5) = 'Box2y'
        export_data(5) = position_y(2)
        
        data_headers(6) = 'M0_star_1'
        export_data(6) = mass_star1*SymmFactor
        
        data_headers(7) = 'M0_star_2'
        export_data(7) = mass_star2*SymmFactor
     end if

!     data_headers(4) = 'Box2x'
!     export_data(4) = position_x(3)
!     data_headers(5) = 'Box2y'
!     export_data(5) = position_y(3)

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)


     !--------------------------------------------------------------------------------!
     !************ Density mode decomposition ************!

     filename = "bhns-dens_mode.con"
     num_cols = 14
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'm=0'
     export_data(2) = density_modes_r0*SymmFactor

     data_headers(3) = 'm=1 real'
     export_data(3) = density_modes_r1*SymmFactor

     data_headers(4) = "m=1 im"
     export_data(4) = density_modes_i1*SymmFactor

     data_headers(5) = "m=2 real"
     export_data(5) = density_modes_r2*SymmFactor

     data_headers(6) = "m=2 im"
     export_data(6) = density_modes_i2*SymmFactor

     data_headers(7) = "m=3 real"
     export_data(7) = density_modes_r3*SymmFactor

     data_headers(8) = "m=3 im"
     export_data(8) = density_modes_i3*SymmFactor

     data_headers(9) = "m=4 real"
     export_data(9) = density_modes_r4*SymmFactor

     data_headers(10) = "m=4 im"
     export_data(10) = density_modes_i4*SymmFactor

     data_headers(11) = "m=5 real"
     export_data(11) = density_modes_r5*SymmFactor

     data_headers(12) = "m=5 im"
     export_data(12) = density_modes_i5*SymmFactor

     data_headers(13) = "m=6 real"
     export_data(13) = density_modes_r6*SymmFactor

     data_headers(14) = "m=6 im"
     export_data(14) = density_modes_i6*SymmFactor
     

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     !--------------------------------------------------------------------------------!

     !--------------------------------------------------------------------------------!
     !************ HAMILTONIAN CONSTRAINT ************!

     filename = "bhns-ham.con"
     num_cols = 11
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

     data_headers(6) = "HR1.1-NoD"
     export_data(6) = Ham_const_excised_VolIntN/Ham_const_excised_VolIntD

     data_headers(7) = "HRes_outAH_N"
     !export_data(7) = Ham_const_outside_AH_mask_VolIntN*SymmFactor

     data_headers(8) = "HRes_outAH_D"
     !export_data(8) = Ham_const_outside_AH_mask_VolIntD*SymmFactor

     data_headers(9) = "HRes_outAHNoD"
     !export_data(9) = Ham_const_outside_AH_mask_VolIntN/Ham_const_outside_AH_mask_VolIntD

     data_headers(10) = "HRes2P1.1NNR"
     export_data(10) = Ham_const_excised_innerregion_VolIntN*SymmFactor

     data_headers(11) = "HRes2P1.1DNR"
     export_data(11) = Ham_const_excised_innerregion_VolIntD*SymmFactor

     

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     !--------------------------------------------------------------------------------!
     !************ MOMENTUM CONSTRAINT ************!

     filename = 'bhns-mom.con'
     num_cols = 10
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

     data_headers(10) = "MRx 1.1NoD"
     export_data(10) = momx_const_excised_VolIntN/mom_const_excised_VolIntD

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ GAMMA CONSTRAINT ************!
     filename = 'bhns-gam.con'
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
!!$     call CCTK_ReductionHandle(handle,"minimum")
!!$     call CCTK_VarIndex(vindex,"grid::Z")
!!$     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Zmin,1,vindex)
!!$
!!$     do jj=1,moncrief_gw_num_radii
!!$
!!$        moncrief_gw_extr_radius = moncrief_radius_GW(jj)        
!!$
!!$        ! Gravitational wave extraction!
!!$        Mass_GW = 0.0D0
!!$
!!$        ! Initialize several variables at t=0
!!$        if (cctk_iteration==0) then
!!$           momentsr_old = 0.d0
!!$           momentsi_old = 0.d0
!!$           odd_momentsr_old = 0.d0
!!$           odd_momentsi_old = 0.d0
!!$           int_momentsr = 0.d0
!!$           int_momentsi = 0.d0
!!$           momentsr = 0.D0
!!$           momentsi = 0.D0
!!$        end if
!!$
!!$        !gw_amplitudes = 0.D0
!!$        hplus = 0.d0
!!$        hcross = 0.d0
!!$
!!$        nmodes_moncrief = ubound(momentsr,1)
!!$
!!$        call gw_wave_flux(cctkGH,moncrief_gw_extr_radius,Mass_GW, &
!!$             dT,momentsr,momentsi, &
!!$             momentsr_old,momentsi_old,odd_momentsr,odd_momentsi, &
!!$             odd_momentsr_old,odd_momentsi_old,int_momentsr,int_momentsi, &
!!$             hplus,hcross,theta_GW,phi_GW, &
!!$             Zmin,dZ)
!!$
!!$        if(moncrief_gw_num_radii.eq.1) then
!!$           filename = 'bhns-evenmodes.wave'
!!$        else
!!$           write(filename,41)jj
!!$41         FORMAT('bhns-evenmodes.wave.',I1)
!!$        endif
!!$
!!$        num_cols = 2*nmodes_moncrief + 1
!!$
!!$        data_headers(1) = "# Time"
!!$        export_data(1) = CCTK_TIME
!!$
!!$        do ii=1,nmodes_moncrief
!!$           data_headers(2*ii) = " " 
!!$           data_headers(2*ii+1) = " " 
!!$           export_data(2*ii) = momentsr(ii)
!!$           export_data(2*ii+1) = momentsi(ii)
!!$        end do
!!$
!!$        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
!!$        !--------------------------------------------------------------------------------!
!!$
!!$        if(moncrief_gw_num_radii.eq.1) then
!!$           filename = 'bhns-oddmodes.wave'
!!$        else
!!$           write(filename,51)jj
!!$51         FORMAT('bhns-oddmodes.wave.',I1)
!!$        endif
!!$        num_cols = 2*nmodes_moncrief + 1
!!$
!!$        data_headers(1) = "# Time"
!!$        export_data(1) = CCTK_TIME
!!$
!!$        do ii=1,nmodes_moncrief
!!$           data_headers(2*ii) = " "
!!$           data_headers(2*ii+1) = " "
!!$           export_data(2*ii) = odd_momentsr(ii)
!!$           export_data(2*ii+1) = odd_momentsi(ii)
!!$        end do
!!$
!!$        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
!!$
!!$     end do
!!$

     !--------------------------------------------------------------------------------!
     ! BH local diagnostics
     ! ihf = Isolated horizon formalism
     if (Symmetry==EQUATORIAL)  then 
     filename = 'bhns-ihf-diags.mon'
     num_cols = 5
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     N_theta_bh_diagnostics =  50
     N_phi_bh_diagnostics   =  50

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

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
!  end if
     !--------------------------------------------------------------------------------!
     !************ EM DIAGNOSTICS  ************!
     if(em_evolve_enable==1) then

        filename = 'bhns-emf.con'
        num_cols = 16
        export_data = 0.D0

        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME

        ! EM energy volume integral is computed in scheduler of diagnostics_mhd
        data_headers(2) = 'Eem'
        export_data(2) = em_energy_VolInt*SymmFactor

        ! EM energy volume integral is computed in scheduler of diagnostics_mhd
        data_headers(3) = 'Eem_outsideAH'
        export_data(3) = em_energy_outsideBH_VolInt*SymmFactor

        ! (Alt.) EM energy volume integral is computed in scheduler of diagnostics_mhd <-- CRAPPY, don't use
        !data_headers(3) = 'Eem2'
        !export_data(3) = em_energy2_VolInt*SymmFactor

        ! b_phi volume integral is computed in scheduler of diagnostics_mhd
        data_headers(4) = 'int B_phi'
        export_data(4) = b_phi_VolInt*SymmFactor

        data_headers(5) = 'Bx_max'
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex (vindex, "mhd_evolve::Bx")
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(5),1,vindex)

        data_headers(6) = 'By_max'
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex (vindex, "mhd_evolve::By")
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(6),1,vindex)

        data_headers(7) = "Monopole"
        export_data(7) = monopole_VolInt*SymmFactor

        data_headers(8) = "Monop.t=0"
        export_data(8) = initial_monopole_value*SymmFactor

        data_headers(9) = "Mpole_outAH"
        export_data(9) = monopole_outsideBH_VolInt*SymmFactor

        data_headers(10) = 'b2_max'
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex (vindex, "mhd_evolve::smallb2")
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

        ! EM energy volume integral is computed in scheduler of diagnostics_mhd
        data_headers(14) = 'Eemoutside_r1'
        export_data(14) = em_energy_outsideradius1_VolInt*SymmFactor

        ! EM energy volume integral is computed in scheduler of diagnostics_mhd
        data_headers(15) = 'Eemoutside_r2'
        export_data(15) = em_energy_outsideradius2_VolInt*SymmFactor

        data_headers(16) = 'b2_maxv2'
        export_data(16)  = bhns_max_b2

        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     end if


     !************ RAD DIAGNOSTICS  ************!
        filename = 'bhns-rad.con'
        num_cols = 5
        export_data = 0.D0

        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME

        ! Rad energy volume integral is computed in scheduler of diagnostics_mhd
        data_headers(2) = 'Erad'
        export_data(2) = rad_energy_VolInt*SymmFactor

        data_headers(3) = 'Erad_nue'
        export_data(3) = rad_energy_nue_VolInt*SymmFactor

        data_headers(4) = 'Erad_nux'
        export_data(4) = rad_energy_nux_VolInt*SymmFactor

	data_headers(5) = 'T_fluid_max'
     	call CCTK_ReductionHandle(handle,"maximum")
     	call CCTK_VarIndex (vindex, "mhd_evolve::T_fluid")
     	call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(3),1,vindex)



	if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !  alpha diagnostic
     if(alpha_diagnostic.eq.1) then
        if(mod(cctk_iteration,out_every)==0) then    
        
           call bhns_surface_dens_1d(cctkGH,CCTK_TIME,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1), &
                 surfxnum,surfphinum,surfznum,surfzmax,surfxmin,surfxmax,Symmetry,CCTK_MyProc(CCTKGH))
        
           call bhns_compute_alpha_1d(cctkGH,CCTK_TIME,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1), &
                surfxnum,surfphinum,surfznum,surfzmax,surfxmin,surfxmax,Symmetry,CCTK_MyProc(CCTKGH),alpha_rho_cut_off)
        end if
     end if
     !--------------------------------------------------------------------------------!
  end if

end if
!end if
end subroutine bhns_diagnostics
